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
# with c = Aᵀ W r(z), where ∂g (the caller's partials of g) is plain arithmetic in their dual numbers
# at the fitted, primal z. The second term of H, the residuals' curvature, does not vanish when the
# residuals do not: without it the derivatives are those of the Gauss-Newton approximation.

# The weights of a least-squares fit of `n` quotes: one finite, non-negative number per quote.
__least_squares_weights(::Nothing, n) = ones(n)
function __least_squares_weights(w::AbstractVector, n)
    length(w) == n || throw(DimensionMismatch("Fit.LeastSquares: $(length(w)) weights for $n quotes"))
    all(x -> isfinite(x) && x >= 0, __primal.(w)) || throw(
        ArgumentError("Fit.LeastSquares: weights must be finite and non-negative; got $(w)")
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
    A = SparseArrays.sparse(rows, cols, convert(Vector{__float_eltype(amounts)}, amounts), length(quotes), length(times))
    return A, times
end

# The residuals of a least-squares fit of fixed-cashflow quotes, on primal (or the solve's own dual) knot
# rates `z`: their cashflow matrix `A` (primal), its distinct `times`, and the quotes' `prices`.
struct __CashflowResiduals{S, I, T, E, M, V}
    spline::S
    interpolation::I
    tenors::T
    extrapolation::E
    A::M
    times::V
    prices::V
end
function __discounts(r::__CashflowResiduals, z)
    curve = __trial_curve(r.spline, r.interpolation, z, r.tenors, r.extrapolation)
    return map(t -> discount(curve, t), r.times)
end
__residuals(r::__CashflowResiduals, z) = r.A * __discounts(r, z) .- r.prices
__residual_jacobian(r::__CashflowResiduals, z) = r.A * ForwardDiff.jacobian(z -> __discounts(r, z), z)

__quote_prices(r::__CashflowResiduals) = r.prices

# Quotes whose cashflows depend on the curve are priced one by one, as an exact fit prices them.
__residuals(r::__RepricingResiduals, z) = r(z)
__residual_jacobian(r::__RepricingResiduals, z) = ForwardDiff.jacobian(r, z)
__quote_prices(r::__RepricingResiduals) = [q.price for q in r.quotes]

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

# Levenberg-Marquardt: iteration limit; the cosine between the residuals and every Jacobian column at
# which it has converged; the largest cosine it accepts when it stops otherwise (a step or a damping at
# rounding); and the step (in knot rates) below which it stops. A fit whose weighted residuals are all
# within `__EXACT_FIT_FTOL` of their quotes' sizes has converged too: it reprices its quotes, and the
# cosine of residuals at rounding says nothing.
const __LSQ_MAXITERS = 200
const __LSQ_GTOL = 1.0e-12
const __LSQ_ACCEPT = 1.0e-8
const __LSQ_STEPTOL = 1.0e-15

# The largest cosine between the residuals `f` and a column of `J`: zero at a minimum, whatever the
# scale of the residuals or of the knot rates.
function __gradient_cosine(J, f)
    nf = norm(f)
    iszero(nf) && return zero(nf)
    return maximum(k -> (c = view(J, :, k); abs(dot(c, f)) / (norm(c) * nf)), axes(J, 2))
end

# Minimize the sum of the squared residuals `r`, each multiplied by `sw` (the square roots of the weights),
# from the seed `z0`. Steps solve the damped normal equations (JᵀJ + μD²)δ = -Jᵀf as a least-squares
# problem (QR), with D the largest column norms of J seen so far (Marquardt's scaling); a step is kept if
# it decreases the sum, and μ follows the ratio of the actual to the predicted decrease (Nielsen's rule).
function __solve_least_squares(r, sw, z0, tenors, who)
    ReturnCode = Optimization.SciMLBase.ReturnCode
    weighted(z) = sw .* __residuals(r, z)
    jacobian(z) = sw .* __residual_jacobian(r, z)
    z = copy(z0)
    r0 = __residuals(r, z)
    f = sw .* r0
    J = jacobian(z)
    # each quote's size, as an exact fit scales its residuals: its price, or its value at the seed
    scale = map(__quote_prices(r), r0) do p, rᵢ
        s = max(abs(p), abs(rᵢ + p))
        iszero(s) ? one(s) : s
    end
    reprices(f) = all(i -> abs(f[i]) <= __EXACT_FIT_FTOL * sw[i] * scale[i], eachindex(f))
    __check_determined(J, tenors, who)
    n = length(z)
    D = [norm(view(J, :, k)) for k in 1:n]
    ϕ = sum(abs2, f)
    μ, ν = 1.0e-3, 2.0
    retcode = ReturnCode.MaxIters
    cosine = __gradient_cosine(J, f)
    for _ in 1:__LSQ_MAXITERS
        (cosine <= __LSQ_GTOL || reprices(f)) && (retcode = ReturnCode.Success; break)
        δ = -(vcat(J, sqrt(μ) .* Diagonal(D)) \ vcat(f, zeros(n)))
        if maximum(abs, δ) <= __LSQ_STEPTOL
            retcode = cosine <= __LSQ_ACCEPT ? ReturnCode.Success : ReturnCode.Stalled
            break
        end
        zt = z .+ δ
        ft = weighted(zt)
        ϕt = sum(abs2, ft)
        predicted = ϕ - sum(abs2, f .+ J * δ)
        ρ = (ϕ - ϕt) / predicted   # NaN or negative for a non-finite or worse trial
        if predicted > 0 && ρ > 1.0e-4
            z, f, ϕ = zt, ft, ϕt
            J = jacobian(z)
            D .= max.(D, [norm(view(J, :, k)) for k in 1:n])
            cosine = __gradient_cosine(J, f)
            μ *= max(1 / 3, 1 - (2ρ - 1)^3)
            ν = 2.0
        else
            μ *= ν
            ν *= 2
            if μ > 1.0e16
                retcode = cosine <= __LSQ_ACCEPT ? ReturnCode.Success : ReturnCode.Stalled
                break
            end
        end
    end
    # Stopped where the sum of squares no longer registers a decrease (its rounding), short of the gradient
    # tolerance: Gauss-Newton steps still move the knot rates towards the minimum by the residuals'
    # curvature ratio each, which the cosine (not the sum) shows.
    if retcode == ReturnCode.Success && cosine > __LSQ_GTOL && !reprices(f)
        for _ in 1:5
            zt = z .- (J \ f)
            ft = weighted(zt)
            sum(abs2, ft) <= ϕ * (1 + 1.0e-12) || break
            Jt = jacobian(zt)
            ct = __gradient_cosine(Jt, ft)
            ct < cosine || break
            z, f, J, cosine, ϕ = zt, ft, Jt, ct, sum(abs2, ft)
            cosine <= __LSQ_GTOL && break
        end
    end
    retcode == ReturnCode.Success && return z
    stopped = retcode == ReturnCode.MaxIters ? "did not converge within $(__LSQ_MAXITERS) iterations" :
        "cannot reduce its quotes' weighted squared pricing residuals further"
    throw(
        FitConvergenceError(
            retcode,
            "$who $stopped: the residuals and the Jacobian of the knot rates are not orthogonal (largest " *
                "cosine $cosine, at most $(__LSQ_ACCEPT) accepted). Check the quotes for prices that no curve " *
                "with this interpolation comes near, or use fewer knots."
        )
    )
end

# Fit the knot rates at `tenors` to `quotes` by weighted least squares (`Fit.LeastSquares`), with the knot
# grid validated and seeded as a spline fit's, and return the curve, differentiated by the implicit
# function theorem when the quotes or the weights carry dual numbers.
function __fit_least_squares(spline, tenors, quotes, primal_quotes, weights, extrapolation)
    who = "fit($(spline))"
    length(primal_quotes) >= length(tenors) || throw(
        ArgumentError(
            "$who: a least-squares fit of $(length(tenors)) knots needs at least as many quotes; got " *
                "$(length(primal_quotes)). Use fewer knots."
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
    fixed = all(q -> __fixed_cashflows(q.instrument), primal_quotes)
    r = if fixed
        A, times = __cashflow_matrix(primal_quotes)
        prices = float.([q.price for q in primal_quotes])
        # a dual number in a contract that `__primal_contract` does not strip
        (__ad_depth(eltype(A)) == 0 && __ad_depth(eltype(prices)) == 0) || throw(__foreign_dual_error(who, __SPLINE_FIT_HINT))
        __CashflowResiduals(spline, interpolation, grid0.tenors, extrapolation, A, times, convert(typeof(times), prices))
    else
        __RepricingResiduals(spline, interpolation, grid0.tenors, primal_quotes, extrapolation, who)
    end
    rates = __solve_least_squares(r, sqrt.(float.(__primal.(weights))), grid0.rates, grid0.tenors, who)
    curve = Yield.ZeroRateCurve(rates, grid0.tenors, spline; extrapolation)   # public result: validated (finite rates)
    (all(splat(===), zip(quotes, primal_quotes)) && __ad_depth(eltype(weights)) == 0) && return curve
    fixed || throw(
        ArgumentError(
            "$who cannot differentiate a least-squares fit of quotes whose cashflows depend on the curve " *
                "(such as a Bond.Floating): it differentiates fits of fixed cashflows (Cashflow and Bond.Fixed " *
                "instruments, which every FinanceModels quote type is)."
        )
    )
    return __implicit_least_squares(curve, r, quotes, weights)
end

# The fitted `curve` with knot rates carrying the derivatives of the least-squares fit (see the top of
# this file) with respect to the dual numbers in `quotes` (prices and cashflow amounts) and `weights`.
# `r` holds the fit's primal cashflow matrix and prices.
function __implicit_least_squares(curve, r::__CashflowResiduals, quotes, weights)
    A, times = __cashflow_matrix(quotes)
    times == r.times || throw(ArgumentError("fit: the quotes' cashflow times changed between their primal and dual copies"))
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
    H = J' * (w .* J) .+ ForwardDiff.hessian(z -> dot(c, __discounts(r, z)), z)
    F = lu(H; check = false)
    (issuccess(F) && cond(H) <= 1.0e12) || throw(
        ArgumentError(
            "fit cannot differentiate this least-squares fit: the Hessian of its objective is singular or " *
                "ill-conditioned (condition number $(cond(H))) at the fitted knot rates."
        )
    )
    # The optimality condition at the fitted, primal rates: its value is zero to the fit's precision, and
    # its partials are those of the caller's dual numbers.
    g = ∂d' * (A' * (weights .* (A * d .- p)))
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
