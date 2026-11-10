# ── Least-squares knot-curve fits ──────────────────────────────────────────────────────
#
# `Fit.LeastSquares` fits a knot curve to quotes that no curve need reprice: more quotes than knots,
# repeated maturities, the bonds of many issuers. The knot rates z minimize the weighted squared
# pricing residuals
#
#     ϕ(z) = Σᵢ wᵢ rᵢ(z)²,   rᵢ(z) = present_value(curve(z), qᵢ.instrument) - qᵢ.price,
#
# by Levenberg-Marquardt on the residuals' Jacobian (one dual layer of the knots' width, where a loss
# fit's Newton steps need the loss's Hessian). When every quote's cashflows are fixed, the prices are the
# cashflow matrix times the curve's discount factors at the distinct cashflow times, r(z) = A d(z) - p:
# the curve is evaluated once per distinct time rather than once per cashflow, and the Jacobian is
# A ∂d/∂z. The caller's dual numbers then enter only through A, p and the weights, so the fit is
# differentiated by the implicit function theorem on its optimality condition
#
#     g(z) = (∂d/∂z)ᵀ Aᵀ W r(z) = 0:   dz = -H⁻¹ ∂g,   H = (A ∂d/∂z)ᵀ W (A ∂d/∂z) + Σⱼ cⱼ ∇²dⱼ(z),
#
# where H is the Hessian of ϕ/2 and c = Aᵀ W r(z), where ∂g (the caller's partials of g) is plain arithmetic in their dual numbers
# at the fitted, primal z. The second term of H, the residuals' curvature, does not vanish when the
# residuals do not: without it the derivatives are those of the Gauss-Newton approximation.

# The weights of a least-squares fit of `n` quotes: one finite, non-negative number per quote.
__least_squares_weights(::Nothing, n) = ones(n)
function __least_squares_weights(w::AbstractVector, n)
    all(x -> isfinite(x), __primal.(w)) || throw(
        ArgumentError("Fit.LeastSquares: weights must be finite; got $(w)")
    )
    return w
end

# The quotes' cashflows as a sparse matrix over their distinct times: the price of quote i on a curve is
# Σⱼ A[i, j] discount(curve, times[j]). Amounts keep their type (dual numbers included).
function __cashflow_matrix(quotes)
    cfs = map(q -> collect(q.instrument), quotes)
    # A time of -0.0 is the instant 0.0; adding zero makes it so for `unique` and the lookup.
    timepoint(cf) = (__no_dual_time(cf.time); float(cf.time) + 0.0)
    times = sort!(unique!([timepoint(cf) for c in cfs for cf in c]))
    amounts = [cf.amount for c in cfs for cf in c]
    rows = [i for (i, c) in enumerate(cfs) for _ in c]
    cols = [searchsortedfirst(times, timepoint(cf)) for c in cfs for cf in c]
    values = convert(Vector{__float_eltype(amounts)}, amounts)
    A = SparseArrays.sparse(rows, cols, values, length(quotes), length(times))
    # Take absolute amounts before coalescing payments at the same time, matching __quote_size.
    gross = any(x -> __primal(x) < 0, values) ?
        SparseArrays.sparse(rows, cols, abs.(values), length(quotes), length(times)) : nothing
    return A, times, gross
end

# The residuals of a least-squares fit of fixed-cashflow quotes, on primal (or the solve's own dual) knot
# rates `z`: their cashflow matrix `A` (primal), its distinct `times`, and the quotes' `prices`. `gross` is
# the matrix of the cashflows' absolute amounts when some quote has cashflows of both signs, and `nothing`
# otherwise.
struct __CashflowResiduals{S, I, T, E, M, V, G}
    spline::S
    interpolation::I
    tenors::T
    extrapolation::E
    A::M
    times::V
    prices::V
    gross::G
end
function __discounts(r::__CashflowResiduals, z)
    curve = __trial_curve(r.spline, r.interpolation, z, r.tenors, r.extrapolation)
    return map(t -> discount(curve, t), r.times)
end
# The discount factors at the cashflow times as a callable, whose type tags the Jacobian's dual numbers.
struct __DiscountsAt{R}
    r::R
end
(f::__DiscountsAt)(z) = __discounts(f.r, z)

# The residuals at `z`, with what the quotes' sizes need (the discount factors, or nothing).
function __evaluate(r::__CashflowResiduals, z)
    d = __discounts(r, z)
    return r.A * d .- r.prices, d
end

# The Jacobian's configuration, built once per fit: its chunk follows the number of knots, so its type is
# known only at run time, and the solve asserts the Jacobian's type `T` to stay inferred (as the exact
# solve does).
__jacobian_config(r::__CashflowResiduals, z0) = ForwardDiff.JacobianConfig(__DiscountsAt(r), z0)
__residual_jacobian(r::__CashflowResiduals, z, config, ::Type{T}) where {T} = r.A * (ForwardDiff.jacobian(__DiscountsAt(r), z, config)::Matrix{T})

# The quotes must determine every knot rate: each knot moves some quote's price, and the Jacobian with its
# columns scaled to unit length is well conditioned. `J` is the weighted residuals' Jacobian.
function __check_determined(J, tenors, who)
    all(isfinite, J) || throw(
        FitConvergenceError(
            Optimization.SciMLBase.ReturnCode.InternalLinearSolveFailed,
            "$who: the Jacobian of the quotes' pricing residuals is not finite at the starting knot rates"
        )
    )
    norms = [norm(view(J, :, k)) for k in axes(J, 2)]
    k = findfirst(iszero, norms)
    k === nothing || throw(
        ArgumentError(
            "$who: the knot at $(tenors[k]) does not move any quote's price (with a nonzero weight), so the " *
                "quotes do not determine its rate. Remove the knot, or add quotes whose cashflows depend on it."
        )
    )
    κ = cond(J ./ permutedims(norms))
    κ <= 1.0e10 || throw(
        ArgumentError(
            "$who: the quotes barely determine the knot rates (the Jacobian of their pricing residuals has " *
                "condition number $κ). Use fewer knots, or knots where the quotes' cashflows are."
        )
    )
    return nothing
end

# Levenberg-Marquardt can stall when the reductions of the sum of squares ϕ reach floating-point
# precision. The fit then polishes with contracting Gauss-Newton steps, then Newton steps with the full
# Hessian of ϕ/2 (JᵀJ plus the residuals' curvature S). It accepts on the full Newton correction: a small
# Gauss-Newton correction can hide a substantial error in the knot rates when S nearly cancels JᵀJ. A fit
# that reprices its quotes (weighted residuals within `__EXACT_FIT_FTOL` of their sizes) has converged
# without this check: its residuals' curvature is zero to rounding.
const __LSQ_GTOL = 1.0e-12      # LM stops at this cosine between the residuals and the Jacobian's columns
const __LSQ_STEPTOL = 1.0e-15   # LM stops at this step, in knot rates
const __LSQ_FLAT = 1.0e-12      # relative change of ϕ treated as rounding (LM's stop; polishing's allowance)
const __LSQ_POLISH = 10         # Gauss-Newton steps, while they contract and exceed __LSQ_ACCEPT
const __LSQ_NEWTON = 10         # Newton steps, each kept if ϕ stays within __LSQ_FLAT and the gradient falls
const __LSQ_HALVINGS = 30       # halvings of a Newton step before it is given up
const __LSQ_ACCEPT = 1.0e-10    # largest accepted Newton correction, in knot rates (1e-6 bp)

# The largest cosine between the residuals `f` and a column of `J`: zero at a stationary point, whatever
# the scale of the residuals or of the knot rates.
function __gradient_cosine(J, f)
    nf = norm(f)
    iszero(nf) && return zero(nf)
    return maximum(k -> (c = view(J, :, k); abs(dot(c, f)) / (norm(c) * nf)), axes(J, 2))
end

# The least-squares solve's vector arithmetic, which does not depend on the quotes' types, kept out of
# line as the exact solve's is. `__lm_step` solves the damped normal equations (JᵀJ + μD²)δ = -Jᵀf as a
# least-squares problem (QR); `__gauss_newton_step` is its undamped step, from the pivoted QR factorization
# `F` of J.
@noinline __weighted(sw, x) = sw .* x
@noinline __lm_step(J, f, μ, D) = -(vcat(J, sqrt(μ) .* Diagonal(D)) \ vcat(f, zeros(length(D))))
@noinline __qr(J) = qr(J, ColumnNorm())
@noinline __gauss_newton_step(F, f) = -(F \ f)
@noinline __predicted_decrease(ϕ, f, J, δ) = ϕ - sum(abs2, f .+ J * δ)
@noinline __column_norms(J) = [norm(view(J, :, k)) for k in axes(J, 2)]
@noinline __widen_scaling!(D, J) = (D .= max.(D, __column_norms(J)))
@noinline __gradient(J, f) = J' * f
# The weighted residuals, their Jacobian and the gradient Jᵀf at `z`.
function __lsq_point(r, sw, z, config, ::Type{T}) where {T}
    f = __weighted(sw, first(__evaluate(r, z)))
    J = __weighted(sw, __residual_jacobian(r, z, config, T))
    return f, J, __gradient(J, f)
end

# Residual curvature, shared by the solver's Newton polishing and the implicit sensitivities.
# Only the discount factors carry knot duals; the residual weights stay fixed in this Hessian.
function __curvature(r, sw, f, z)
    c = r.A' * __weighted(sw, f)
    return ForwardDiff.hessian(x -> dot(c, __discounts(r, x)), z)
end
__lsq_accepted(δ) = maximum(abs, δ) <= __LSQ_ACCEPT

# The Newton step with JᵀJ + S, the Hessian of half the sum of squares, through the pivoted QR
# factorization `F` of J (J P = QR):
# with C = R⁻ᵀ PᵀSP R⁻¹, the Hessian is P Rᵀ(I + C)R Pᵀ, so the step solves (I + C) R Pᵀδ = -(Qᵀf)₁…ₙ, and
# the Hessian is positive definite where I + C is. This holds where JᵀJ itself is too ill-conditioned to
# factor. `nothing` when I + C is not positive definite, which the fit's acceptance requires.
@noinline function __lsq_newton_step(F, f, S)
    n = size(S, 1)
    R = UpperTriangular(F.R)
    y = -(F.Q' * f)[1:n]
    C = R' \ (R' \ S[F.p, F.p])'   # R⁻ᵀ (R⁻ᵀ Sₚ)ᵀ = R⁻ᵀ Sₚ R⁻¹, Sₚ symmetric
    M = (C .+ C') ./ 2 + I
    all(isfinite, M) || return nothing
    K = cholesky(Symmetric(M); check = false)
    issuccess(K) || return nothing
    δ = similar(y)
    δ[F.p] = R \ (K \ y)
    return δ
end
# Whether every weighted residual `f` is within `__EXACT_FIT_FTOL` of its quote's size: for quotes of
# nonnegative cashflows (no `gross` matrix), their value, checked without allocating, which stops at the
# first quote a fit with more quotes than knots does not reprice.
function __reprices(r::__CashflowResiduals{S, I, T, E, M, V, Nothing}, f, sw, ru, d) where {S, I, T, E, M, V}
    for i in eachindex(f, sw, ru)
        abs(f[i]) <= __EXACT_FIT_FTOL * sw[i] * abs(ru[i] + r.prices[i]) || return false
    end
    return true
end
function __reprices(r, f, sw, ru, d)
    sizes = r.gross * d
    for i in eachindex(f, sw, sizes)
        abs(f[i]) <= __EXACT_FIT_FTOL * sw[i] * sizes[i] || return false
    end
    return true
end

# Minimize the sum of the squared residuals `r`, each multiplied by `sw` (the square roots of the weights),
# from the seed `z0`, in at most `maxiters` iterations. Steps solve the damped normal equations with D the
# largest column norms of J seen so far (Marquardt's scaling); a step is kept if it decreases the sum, and
# μ follows the ratio of the actual to the predicted decrease (Nielsen's rule).
function __lm_solve(r, sw, z0, tenors, maxiters, who, config)
    ReturnCode = Optimization.SciMLBase.ReturnCode
    z = copy(z0)
    ru, d = __evaluate(r, z)
    T = eltype(ru)
    f = __weighted(sw, ru)
    J = __weighted(sw, __residual_jacobian(r, z, config, T))
    __check_determined(J, tenors, who)
    D = __column_norms(J)
    ϕ = sum(abs2, f)
    μ, ν = 1.0e-3, 2.0
    retcode = ReturnCode.MaxIters
    for _ in 1:maxiters
        __reprices(r, f, sw, ru, d) && (retcode = ReturnCode.Success; break)
        __gradient_cosine(J, f) <= __LSQ_GTOL && (retcode = ReturnCode.Stalled; break)
        δ = __lm_step(J, f, μ, D)
        maximum(abs, δ) <= __LSQ_STEPTOL && (retcode = ReturnCode.Stalled; break)
        predicted = __predicted_decrease(ϕ, f, J, δ)
        predicted <= __LSQ_FLAT * ϕ && (retcode = ReturnCode.Stalled; break)
        zt = __fit_trial(z, 1.0, δ)
        rut, dt = __evaluate(r, zt)
        ft = __weighted(sw, rut)
        ϕt = sum(abs2, ft)
        ρ = (ϕ - ϕt) / predicted
        if ρ > 1.0e-4
            z, ru, d, f, ϕ = zt, rut, dt, ft, ϕt
            J = __weighted(sw, __residual_jacobian(r, z, config, T))
            __widen_scaling!(D, J)
            μ *= max(1 / 3, 1 - (2ρ - 1)^3)
            ν = 2.0
        else
            μ *= ν
            ν *= 2
            μ > 1.0e16 && (retcode = ReturnCode.Stalled; break)
        end
    end
    return (; z, f, J, retcode)
end

# Contracting Gauss-Newton steps recover progress after the objective reaches rounding.
function __gauss_newton_polish(r, sw, z, f, J, config)
    T = eltype(f)
    F = __qr(J)
    δ = __gauss_newton_step(F, f)
    ϕ = sum(abs2, f)
    for _ in 1:__LSQ_POLISH
        step = maximum(abs, δ)
        __lsq_accepted(δ) && break
        zt = __fit_trial(z, 1.0, δ)
        ft = __weighted(sw, first(__evaluate(r, zt)))
        ϕt = sum(abs2, ft)
        ϕt <= ϕ * (1 + __LSQ_FLAT) || break
        Jt = __weighted(sw, __residual_jacobian(r, zt, config, T))
        Ft = __qr(Jt)
        δt = __gauss_newton_step(Ft, ft)
        maximum(abs, δt) <= step / 2 || break
        z, f, J, F, ϕ, δ = zt, ft, Jt, Ft, ϕt, δt
    end
    return (; z, f, J)
end

# Certify with the full AD curvature and polish while the gradient still decreases.
function __newton_polish(r, sw, z, f, J, config, who, maxiters)
    T = eltype(f)
    ϕ, g = sum(abs2, f), __gradient(J, f)
    δ = __gauss_newton_step(__qr(J), f)
    definite, kink = true, false
    for _ in 1:__LSQ_NEWTON
        curve = __trial_curve(r.spline, r.interpolation, z, r.tenors, r.extrapolation)
        kink = Yield.__near_kink(curve, z)
        kink && break
        δn = __lsq_newton_step(__qr(J), f, __curvature(r, sw, f, z))
        definite = !isnothing(δn)
        definite || break
        δ = δn
        # The correction must resolve the same smooth interpolation branch as its Hessian.
        kink = Yield.__near_kink(curve, z, -δ)
        kink && break
        __lsq_accepted(δ) && return __fit_trial(z, 1.0, δ)
        kept = false
        α = 1.0
        for _ in 1:__LSQ_HALVINGS
            zt = __fit_trial(z, α, δ)
            ft, Jt, gt = __lsq_point(r, sw, zt, config, T)
            if sum(abs2, ft) <= ϕ * (1 + __LSQ_FLAT) && norm(gt) < norm(g)
                z, f, J, ϕ, g, kept = zt, ft, Jt, sum(abs2, ft), gt, true
                break
            end
            α /= 2
        end
        kept || break
    end
    curve = __trial_curve(r.spline, r.interpolation, z, r.tenors, r.extrapolation)
    kink = kink || Yield.__near_kink(curve, z, -δ)
    throw(__least_squares_error(who, Optimization.SciMLBase.ReturnCode.Stalled, maxiters, maximum(abs, δ), __gradient_cosine(J, f), kink, definite))
end

function __solve_least_squares(r, sw, z0, tenors, maxiters, who)
    config = __jacobian_config(r, z0)
    point = __lm_solve(r, sw, z0, tenors, maxiters, who, config)
    point.retcode == Optimization.SciMLBase.ReturnCode.Success && return point.z
    if point.retcode == Optimization.SciMLBase.ReturnCode.Stalled
        point = __gauss_newton_polish(r, sw, point.z, point.f, point.J, config)
        return __newton_polish(r, sw, point.z, point.f, point.J, config, who, maxiters)
    end
    δ = __gauss_newton_step(__qr(point.J), point.f)
    throw(__least_squares_error(who, point.retcode, maxiters, maximum(abs, δ), __gradient_cosine(point.J, point.f), false, true))
end
function __least_squares_error(who, retcode, maxiters, step, cosine, kink, definite)
    stopped = retcode == Optimization.SciMLBase.ReturnCode.MaxIters ?
        "did not converge within $maxiters iterations (Fit.LeastSquares(maxiters = …) raises the limit)" :
        "cannot reduce its quotes' weighted squared pricing residuals further"
    where_ = definite ?
        "its last step from the fitted knot rates is $step (at most $(__LSQ_ACCEPT) accepted)" :
        "the Hessian of half the sum of squares is not positive definite where it stopped, so the fit fails its " *
        "required positive-definiteness check"
    cause = if kink
        "The fit stopped where the interpolation switches formulas, so the objective has a kink there: a " *
            "least-squares minimum of a monotone convex, PCHIP or Akima curve can lie on such a switch, where " *
            "this fit cannot converge. Use an interpolation that is smooth in its knot rates (such as " *
            "Spline.Linear() or Spline.Cubic()), knots at the quote maturities (Fit.Exact()), or a minimizer " *
            "that does not use derivatives: fit(curve, quotes, Fit.Loss(); optimizer = NelderMead())."
    else
        "Check the quotes for prices that no curve with this interpolation comes near, or use fewer knots."
    end
    return FitConvergenceError(
        retcode,
        "$who $stopped: $where_, and the largest cosine between the residuals and a column of their " *
            "Jacobian is $cosine (zero at a stationary point). " * cause
    )
end

# Fit the knot rates at `tenors` to `quotes` by weighted least squares (`Fit.LeastSquares`), with the knot
# grid validated and seeded as a spline fit's, and return the curve, differentiated by the implicit
# function theorem when the quotes or the weights carry dual numbers.
function __fit_least_squares(spline, tenors, quotes, primal_quotes, method::Fit.LeastSquares, extrapolation)
    who = "fit($(spline))"
    weights = __least_squares_weights(method.weights, length(primal_quotes))
    w = float.(__primal.(weights))
    weighted = count(>(0), w)
    weighted >= length(tenors) || throw(
        ArgumentError(
            "$who: a least-squares fit of $(length(tenors)) knots needs at least as many quotes with a " *
                "positive weight; got $weighted of $(length(primal_quotes)) quotes. Use fewer knots."
        )
    )
    any(x -> __ad_depth(x) > 0, __extrapolation_dual(extrapolation)) && throw(
        ArgumentError(
            "$who: a least-squares fit cannot differentiate with respect to its extrapolation forward. Fit " *
                "with a primal Yield.FlatForwardAt, or use knots at the quote maturities."
        )
    )
    grid0 = Yield.KnotGrid(__knot_fit_seed(spline, tenors), tenors, spline; who)
    interpolation = Yield.__interpolation(spline, length(grid0.tenors))
    all(q -> __fixed_cashflows(q.instrument), primal_quotes) || throw(
        ArgumentError("$who: Fit.LeastSquares supports fixed cashflows only (Cashflow and Bond.Fixed); use Fit.Loss() for other contracts")
    )
    A, times, gross = __cashflow_matrix(primal_quotes)
    prices = float.([q.price for q in primal_quotes])
    (__ad_depth(eltype(A)) == 0 && __ad_depth(eltype(prices)) == 0) || throw(__foreign_dual_error(who, __SPLINE_FIT_HINT))
    r = __CashflowResiduals(spline, interpolation, grid0.tenors, extrapolation, A, times, convert(typeof(times), prices), gross)
    rates = __solve_least_squares(r, sqrt.(w), grid0.rates, grid0.tenors, method.maxiters, who)
    curve = Yield.ZeroRateCurve(rates, grid0.tenors, spline; extrapolation)   # public result: validated (finite rates)
    # the weights' dual numbers are found by element: a vector of abstract element type can hold them
    (all(splat(===), zip(quotes, primal_quotes)) && !any(x -> __ad_depth(x) > 0, weights)) && return curve
    return __implicit_least_squares(curve, r, quotes, weights)
end

# The fitted `curve` with knot rates carrying the derivatives of the least-squares fit (see the top of
# this file) with respect to the dual numbers in `quotes` (prices and cashflow amounts) and `weights`.
# `r` holds the fit's primal cashflow matrix and prices.
function __implicit_least_squares(curve, r::__CashflowResiduals, quotes, weights)
    A, _, _ = __cashflow_matrix(quotes)
    p = [q.price for q in quotes]
    D = __common_dual_type(Iterators.flatten((SparseArrays.nonzeros(A), p, weights)))
    D === nothing && return curve
    z = FinanceCore.rate.(Yield.knot_rates(curve))
    Yield.__near_kink(curve, z) && throw(__fit_kink_error(curve))
    d = __discounts(r, z)
    ∂d = ForwardDiff.jacobian(z -> __discounts(r, z), z)
    w = float.(__primal.(weights))
    c = r.A' * (w .* (r.A * d .- r.prices))
    J = r.A * ∂d
    H = J' * (w .* J) .+ __curvature(r, sqrt.(w), sqrt.(w) .* (r.A * d .- r.prices), z)
    F = cholesky(Symmetric((H .+ H') ./ 2); check = false)
    (issuccess(F) && cond(H) <= 1.0e12) || throw(
        ArgumentError(
            "fit cannot differentiate this least-squares fit: the Hessian of its objective is not positive " *
                "definite, or is ill-conditioned (condition number $(cond(H))), at the fitted knot rates."
        )
    )
    # The derivatives are those of the minimum only as near it as the fit's acceptance: the Newton step
    # with the exact Hessian is at most `__LSQ_ACCEPT`, as the fit's own was.
    correction = F \ (∂d' * c)
    step = maximum(abs, correction)
    __lsq_accepted(correction) || throw(
        FitConvergenceError(
            Optimization.SciMLBase.ReturnCode.Stalled,
            "fit cannot differentiate this least-squares fit: its knot rates are a Newton step of $step from " *
                "the minimum, more than the $(__LSQ_ACCEPT) its derivatives need."
        )
    )
    # The optimality condition at the fitted, primal rates: its value is zero to the fit's precision, and
    # its partials are those of the caller's dual numbers.
    # (the weights in the dual type, so that a vector of abstract element type multiplies the sparse matrix)
    g = ∂d' * (A' * (map(x -> convert(D, x), weights) .* (A * d .- p)))
    K = ForwardDiff.npartials(D)
    P = [ForwardDiff.partials(convert(D, g[i]), k) for i in eachindex(g), k in 1:K]
    dZ = -(F \ P)
    V = promote_type(ForwardDiff.valtype(D), eltype(z), eltype(dZ))
    Dz = ForwardDiff.Dual{ForwardDiff.tagtype(D), V, K}
    zd = [Dz(V(z[i]), ForwardDiff.Partials{K, V}(ntuple(k -> V(dZ[i, k]), K))) for i in eachindex(z)]
    interpolation = Yield.__interpolation(curve)
    return Yield.__build_public(
        curve.spline, Yield.KnotGrid(zd, Yield.knot_tenors(curve), curve.spline; who = "ZeroRateCurve"), interpolation;
        extrapolation = curve.extrapolation
    )
end
