# ── Least-squares fits ─────────────────────────────────────────────────────────────────
#
# `Fit.LeastSquares` fits any model to any quotes with one solve. A fit has a parameter vector θ and
# a rule that builds the model from it: `__KnotRates` (a knot curve from its knot rates) or
# `__UnitBox` (any other model from its variables). Everything else sees only the residuals
#
#     rᵢ(θ) = √wᵢ · (value(model(θ), quoteᵢ) − priceᵢ)
#
# and never the model or the contracts:
#
# 1. Levenberg–Marquardt on the residuals' ForwardDiff Jacobian (first derivatives only) takes the
#    parameters to where a Gauss–Newton step is at rounding, or no longer shrinks.
# 2. Newton's method with the exact Hessian of the sum of squares certifies a minimum there: the
#    Hessian is positive definite and a Newton correction moves no parameter by more than 1e-10.
# 3. When the quotes carry a caller's dual numbers, the implicit function theorem at that minimum
#    gives its derivatives, with the Hessian there.

# θ = the continuously compounded zero rates at the knots.
struct __KnotRates{S, I, T, E}
    spline::S
    interpolation::I   # `Yield.__interpolation(spline, length(tenors))`, so trial curves have one type
    tenors::T
    extrapolation::E
end
(p::__KnotRates)(z) = __trial_curve(p.spline, p.interpolation, z, p.tenors, p.extrapolation)
__map_parameters(f::F, p::__KnotRates) where {F} =
    __KnotRates(p.spline, p.interpolation, p.tenors, __map_extrapolation(f, p.extrapolation))
__parameter_duals(p::__KnotRates) = __extrapolation_dual(p.extrapolation)
__bounded(::__KnotRates) = false
__parameter_name(p::__KnotRates, k) = "the knot rate at $(p.tenors[k])"
# What to say when a fit of an interpolation with formula switches fails (nothing for the others).
__switch_hint(p::__KnotRates) = p.spline isa Union{Spline.MonotoneConvex, Spline.PCHIP, Spline.Akima} ?
    "$(p.spline) switches formulas as its knot rates change, and a minimum on a switch cannot be reached by " *
    "derivatives (on one, they may not even be finite). Use Spline.Linear() or Spline.Cubic(), or a loss fit " *
    "with an optimizer that does not use derivatives." : ""
__on_kink(p::__KnotRates, z, δ) = Yield.__near_kink(p(z), z, δ)

# θ = the model's variables in AccessibleModels' unit box, each interval mapped to (0, 1).
struct __UnitBox{P}
    variables::P
end
(p::__UnitBox)(x) = AccessibleModels.from_transformed(x, p.variables)
__map_parameters(f::F, p::__UnitBox) where {F} = p
__parameter_duals(::__UnitBox) = ()
__bounded(::__UnitBox) = true
__parameter_name(::__UnitBox, k) = "variable $k"
__switch_hint(::__UnitBox) = ""
__on_kink(::__UnitBox, x, δ) = false

const __LEAST_SQUARES_HINT = "Fit.LeastSquares differentiates quote prices, the coupon rates of Bond.Fixed " *
    "and Bond.Floating, the amounts of Cashflow and FX.BasisSwapLeg (also inside a Composite or a vector), " *
    "its weights, and a Yield.FlatForwardAt forward; a dual number in the model's own fields or inside " *
    "another contract type is not supported."

# The residuals of `quotes` under `parameters(θ)`, each times its `scale` (√weight) unless that is
# `nothing`. `Own` holds the tags of the solve's own dual numbers: any other dual number that reaches
# a residual throws (`__own_duals`).
struct __Residuals{P, Q, S, Own}
    parameters::P
    quotes::Q
    scale::S
    who::String
    __Residuals(parameters::P, quotes::Q, scale::S, ::Type{Own}, who) where {P, Q, S, Own} =
        new{P, Q, S, Own}(parameters, quotes, scale, who)
end
function (r::__Residuals{P, Q, S, Own})(θ::AbstractVector{X}) where {P, Q, S, Own, X}
    model = r.parameters(θ)
    out = Vector{X}(undef, length(r.quotes))
    for i in eachindex(out)
        q = r.quotes[i]
        v = __own_duals(m -> __quote_value(m, q) - q.price, model, X, Own, r.who, __LEAST_SQUARES_HINT)
        out[i] = r.scale === nothing ? v : r.scale[i] * v
    end
    return out
end

__column_norms(J) = [norm(view(J, :, k)) for k in axes(J, 2)]

# The size of a step `δ` from `θ`, with each parameter measured by `D`, its residual column's norm.
__step_size(D, δ, θ) = norm(D .* δ) / (norm(D .* θ) + norm(D))

const __LSQ_FLAT = 1.0e-12       # a relative change of the sum of squares that is rounding
const __LSQ_STEP = 1.0e-10       # the search ends at a Gauss–Newton step of this size
const __LSQ_PATIENCE = 10        # or after this many steps without a smaller one, once the sum of squares is flat
const __LSQ_NEWTON = 1.0e-10     # a certified minimum: a Newton correction of at most this in every parameter
const __LSQ_NEWTON_STEPS = 10

__undetermined(who, parameters, k) = ArgumentError(
    "$who: $(__parameter_name(parameters, k)) does not move any quote's price (with a nonzero weight), " *
        "so the quotes do not determine it. Remove it from the fit, or add quotes that depend on it."
)

# Levenberg–Marquardt from `θ`, with Marquardt's scaling by the largest column norms seen and
# Nielsen's update of the damping `μ`. It ends where the undamped Gauss–Newton step is at most
# `__LSQ_STEP`, or has stopped shrinking, and returns that point without judging it:
# `__certify_minimum` does.
function __levenberg_marquardt(r::R, θ::Vector{T}, cfg::C, bounded::Bool, maxiters::Int) where {R, T, C}
    n = length(θ)
    f = r(θ)
    J = ForwardDiff.jacobian(r, θ, cfg, Val(false))
    nonfinite(at) = FitConvergenceError(
        Optimization.SciMLBase.ReturnCode.Failure,
        String(rstrip("$(r.who): the pricing residuals or their Jacobian are not finite $at. $(__switch_hint(r.parameters))"))
    )
    (all(isfinite, f) && all(isfinite, J)) || throw(nonfinite("at the starting parameters"))
    m = length(f)
    ϕ = sum(abs2, f)
    D = __column_norms(J)
    iszero(maximum(D)) && throw(__undetermined(r.who, r.parameters, 1))
    # A column that is zero here (a Nelson–Siegel decay time with no slope or hump yet) is damped at
    # rounding of the largest until it moves.
    widen!(D, J) = (c = __column_norms(J); D .= max.(D, c, eps(T) * maximum(c)); D)
    widen!(D, J)
    gauss_newton(θ, f, J) = __step_size(D, qr(J, ColumnNorm()) \ f, θ)
    step = least = gauss_newton(θ, f, J)
    waited = 0
    A = Matrix{T}(undef, m + n, n)
    b = zeros(T, m + n)
    μ, ν = T(1.0e-3), T(2)
    iters = 0
    while step > __LSQ_STEP && iters < maxiters
        iters += 1
        # The damped step minimizes ‖f + Jδ‖² + μ‖Dδ‖², as one least-squares problem.
        fill!(A, 0)
        copyto!(view(A, 1:m, :), J)
        for k in 1:n
            A[m + k, k] = sqrt(μ) * D[k]
        end
        copyto!(b, 1, f, 1, m)
        δ = -(qr!(A) \ b)
        # A bounded parameter covers at most 90% of the distance to a bound, so it stays inside.
        bounded && (δ .= clamp.(δ, -0.9 .* θ, 0.9 .* (1 .- θ)))
        predicted = ϕ - sum(abs2, f .+ J * δ)
        # Below this the sum of squares cannot tell a better point from a worse one: steps are then
        # taken while it does not rise, until Gauss–Newton steps stop shrinking.
        flat = predicted <= __LSQ_FLAT * ϕ
        θt = θ .+ δ
        ft = r(θt)
        ϕt = sum(abs2, ft)
        ρ = (ϕ - ϕt) / predicted
        if (predicted > 0 && ρ > 1.0e-4) || (flat && ϕt <= ϕ * (1 + __LSQ_FLAT))
            θ, f, ϕ = θt, ft, ϕt
            J = ForwardDiff.jacobian(r, θ, cfg, Val(false))
            all(isfinite, J) || throw(nonfinite("after $iters iterations"))
            widen!(D, J)
            μ = flat ? μ / 3 : μ * max(one(T) / 3, 1 - (2ρ - 1)^3)
            ν = T(2)
            step = gauss_newton(θ, f, J)
            waited = step < least ? 0 : waited + 1
            least = min(least, step)
            (flat && waited >= __LSQ_PATIENCE) && break
        else
            μ *= ν
            ν *= 2
            μ > 1.0e16 && break
        end
    end
    return (; θ, J, step, iters)
end

# ── The minimum, and its derivatives ───────────────────────────────────────────────────
#
# The gradient of half the weighted sum of squares is
#
#     g(θ, p) = Jₑ(θ, p)ᵀ (w ∘ e(θ, p)),      e = the pricing errors, Jₑ = ∂e/∂θ,
#
# where p stands for a caller's numbers (prices, contract numbers, weights, an extrapolation
# forward), and its Jacobian H = ∂g/∂θ is the full Hessian: JₑᵀWJₑ and the residuals' curvature
# Σ wᵢeᵢ∇²eᵢ, which does not vanish where the residuals do not. A minimum is certified by Newton's
# method on g = 0 (`__certify_minimum`). There, by the implicit function theorem,
# dθ = −H⁻¹ (∂g/∂p) dp (`__implicit_least_squares`).
#
# Both derivatives of g are taken with the solve's own two tags, `inner` for Jₑ and `outer` around
# it. For ∂g/∂p the caller's dual numbers are first rewritten under `outer` (`__Retag`), so their
# tag never meets FinanceModels' own (see `implicit.jl`), and the partials are written back into the
# caller's dual type at the end.

function __least_squares_gradient(e::E, w, θ, inner) where {E}
    J = ForwardDiff.jacobian(e, θ, ForwardDiff.JacobianConfig(e, θ, ForwardDiff.Chunk(θ), inner), Val(false))
    f = e(θ)
    return J' * (w === nothing ? f : w .* f)
end

# The Hessian H at `θ`, as the Cholesky factor `C` of H equilibrated by `d = √diag(H)` (so that no
# check depends on the parameters' units); or, where H is not positive definite and well conditioned,
# the reason in words.
function __hessian_factors(e, w, θ, outer, inner)
    gradient(x) = __least_squares_gradient(e, w, x, inner)
    H = ForwardDiff.jacobian(gradient, θ, ForwardDiff.JacobianConfig(gradient, θ, ForwardDiff.Chunk(θ), outer), Val(false))
    H = (H .+ H') ./ 2
    h = diag(H)
    indefinite = "the Hessian of the sum of squares is not positive definite"
    (all(isfinite, H) && all(>(0), h)) || return indefinite
    d = sqrt.(h)
    Hs = H ./ (d .* d')
    C = cholesky(Symmetric(Hs); check = false)
    issuccess(C) || return indefinite
    κ = cond(Hs)
    κ <= 1.0e10 || return "the Hessian of the sum of squares is ill-conditioned (condition number $κ)"
    return C, d
end

# From where the search ended, Newton's method on g = 0 with the exact Hessian. Returns a certified
# minimum and the last Newton correction `δ`, which reached it. Otherwise it throws: the Hessian must
# be positive definite and well conditioned, and within `__LSQ_NEWTON_STEPS` steps a correction must
# be at most `__LSQ_NEWTON` in every parameter.
function __certify_minimum(e, w, searched, outer, inner, maxiters)
    parameters, who = e.parameters, e.who
    θ = searched.θ
    k = findfirst(iszero, __column_norms(searched.J))
    k === nothing || throw(__undetermined(who, parameters, k))
    stopped = searched.step > __LSQ_STEP && searched.iters >= maxiters
    failed(why) = FitConvergenceError(
        stopped ? Optimization.SciMLBase.ReturnCode.MaxIters : Optimization.SciMLBase.ReturnCode.Stalled,
        "$who " * (stopped ? "did not converge within $maxiters iterations" : "did not reach a minimum") *
            ": $why (after $(searched.iters) Levenberg–Marquardt iterations, where a Gauss–Newton step had relative " *
            "size $(searched.step)). " * __failure_hint(parameters, θ)
    )
    stopped && throw(failed("the search was still moving"))
    sum_of_squares(x) = (f = e(x); w === nothing ? sum(abs2, f) : sum(w .* abs2.(f)))
    ϕ = sum_of_squares(θ)
    local δ
    for steps in 0:__LSQ_NEWTON_STEPS
        factors = __hessian_factors(e, w, θ, outer, inner)
        factors isa String && throw(failed(factors))
        C, d = factors
        δ = -(C \ (__least_squares_gradient(e, w, θ, inner) ./ d)) ./ d
        (__bounded(parameters) && !all(x -> 0 < x < 1, θ .+ δ)) && throw(failed("a Newton step leaves a variable's interval"))
        maximum(abs, δ) <= __LSQ_NEWTON && return θ .+ δ, δ
        steps == __LSQ_NEWTON_STEPS && break
        # A step may not raise the sum of squares beyond rounding. It is halved until it does not:
        # where the curvature changes quickly (near a fold), a full step overshoots.
        t = one(eltype(θ))
        θt = θ .+ δ
        ϕt = sum_of_squares(θt)
        while !(ϕt <= ϕ * (1 + __LSQ_FLAT))
            t /= 2
            t < 1.0e-3 && throw(failed("a Newton step with the full Hessian raises the sum of squares"))
            θt = θ .+ t .* δ
            ϕt = sum_of_squares(θt)
        end
        θ, ϕ = θt, ϕt
    end
    throw(failed("after $(__LSQ_NEWTON_STEPS) Newton steps with the full Hessian, a correction is still $(maximum(abs, δ))"))
end

__failure_hint(p::__KnotRates, z) = isempty(__switch_hint(p)) ?
    "With fewer quotes than knots, or knots that move the quotes' prices only together, the quotes do not " *
    "determine the knot rates: use fewer knots, or knots where the quotes' cashflows are. Otherwise check " *
    "the quotes for prices that no curve with this interpolation comes near." : __switch_hint(p)
function __failure_hint(p::__UnitBox, x)
    k = findfirst(v -> min(v, 1 - v) <= 1.0e-6, x)
    k === nothing || return "Variable $k is at the $(x[k] < 0.5 ? "lower" : "upper") end of its interval: the quotes " *
        "pull it to or beyond that bound. Widen the interval in `variables`."
    return "With fewer quotes than variables, or variables that move the quotes' prices only together, the quotes " *
        "do not determine them. Otherwise check the quotes for prices that the model cannot come near, or start " *
        "from a model nearer the fit."
end

struct __Retag{T} end
(::__Retag{T})(x) where {T} = x
(::__Retag{T})(x::ForwardDiff.Dual{<:Any, V, K}) where {T, V, K} =
    ForwardDiff.Dual{T, V, K}(ForwardDiff.value(x), ForwardDiff.partials(x))

# The minimum `θ` of the primal residuals `e` and weights `w`, reached by the Newton correction `δ`,
# with the derivatives of the calibration with respect to the caller's dual numbers in `quotes`,
# `weights` and `parameters`; `θ` itself when there are none.
function __implicit_least_squares(outer::TO, inner, e::__Residuals{<:Any, <:Any, <:Any, Own}, w, parameters, θ, δ, quotes, weights) where {TO, Own}
    primal_parameters, who = e.parameters, e.who
    # The caller's dual type, from the pricing errors at the fit with their numbers in plain arithmetic.
    model = parameters(θ)
    errors = [__quote_value(model, q) - q.price for q in quotes]
    D = __common_dual_type(Iterators.flatten((errors, something(weights, ()), __parameter_duals(parameters))))
    D === nothing && return θ
    # On a kink of the interpolant the refit has no derivative; nor may the last correction have
    # crossed one, since the Hessian that certified the minimum was taken before it (see
    # `__implicit_knot_curve`).
    __on_kink(primal_parameters, θ, δ) && throw(__fit_kink_error(primal_parameters(θ)))
    # The Hessian at the minimum itself: near a fold it changes within the last correction.
    factors = __hessian_factors(e, w, θ, outer, inner)
    factors isa String && throw(ArgumentError("$who cannot differentiate this fit: $factors at the fitted parameters."))
    C, d = factors
    # g at the fit with the caller's numbers under `outer`: its partials are ∂g/∂p.
    retag = __Retag{TO}()
    K = ForwardDiff.npartials(D)
    V = promote_type(ForwardDiff.valtype(D), eltype(θ))
    et = __Residuals(__map_parameters(retag, parameters), map(q -> __map_quote(retag, q), quotes), nothing, Own, who)
    wt = weights === nothing ? nothing : map(retag, weights)
    g = __least_squares_gradient(et, wt, ForwardDiff.Dual{TO, V, K}.(θ), inner)
    P = [ForwardDiff.partials(g[i], k) / d[i] for i in eachindex(g), k in 1:K]
    dθ = -(C \ P) ./ d
    W = promote_type(V, eltype(dθ))
    Dθ = ForwardDiff.Dual{ForwardDiff.tagtype(D), W, K}
    return [Dθ(W(θ[i]), ForwardDiff.Partials{K, W}(ntuple(k -> W(dθ[i, k]), K))) for i in eachindex(θ)]
end

# Fit `parameters` from `θ0` to `quotes` with the weights and iteration limit of `method`, and return
# the fitted parameter vector. The solve runs on `primal_quotes`; when `quotes`, the weights, or the
# parameters' fixed numbers carry a caller's dual numbers, the vector carries the derivatives of the
# minimum. The solve holds its owner number until it returns or throws.
function __least_squares(parameters, θ0, quotes, primal_quotes, method, who)
    owner = __lease_fit_owner()
    try
        return __least_squares(Val(owner), parameters, θ0, quotes, primal_quotes, method.weights, method.maxiters, who)
    finally
        __return_fit_owner(owner)
    end
end
function __least_squares(::Val{O}, parameters, θ0, quotes, primal_quotes, weights, maxiters, who) where {O}
    θ0 = collect(float.(θ0))
    outer = ForwardDiff.Tag(__FitTagKey{O, :outer}(), eltype(θ0))
    inner = ForwardDiff.Tag(__FitTagKey{O, :inner}(), eltype(θ0))  # created second, so it nests inside `outer`
    Own = Tuple{typeof(outer), typeof(inner)}
    # (`map` over quotes and weights of different lengths would drop the extra ones.)
    weights === nothing || length(weights) == length(quotes) || throw(
        DimensionMismatch("$who: Fit.LeastSquares has $(length(weights)) weights for $(length(quotes)) quotes")
    )
    primal_parameters = __map_parameters(__primal, parameters)
    primal_weights = weights === nothing ? nothing : map(__primal, weights)
    scale = primal_weights === nothing ? nothing : map(sqrt, primal_weights)
    r = __Residuals(primal_parameters, primal_quotes, scale, Own, who)
    cfg = ForwardDiff.JacobianConfig(r, θ0, ForwardDiff.Chunk(θ0), inner)
    searched = __levenberg_marquardt(r, θ0, cfg, __bounded(parameters), maxiters)
    e = __Residuals(primal_parameters, primal_quotes, nothing, Own, who)
    θ, δ = __certify_minimum(e, primal_weights, searched, outer, inner, maxiters)
    differentiated = parameters !== primal_parameters ||
        any(((q, p),) -> q !== p, zip(quotes, primal_quotes)) ||
        (weights !== nothing && any(w -> __ad_depth(w) > 0, weights))
    differentiated || return θ
    return __implicit_least_squares(outer, inner, e, primal_weights, parameters, θ, δ, quotes, weights)
end

# ── `fit` methods ──────────────────────────────────────────────────────────────────────

# Any model: its `variables`, in AccessibleModels' unit box, from the model's own values (a variable
# on a bound starts from the middle of its interval, as in a loss fit).
function fit(mod0, quotes, method::Fit.LeastSquares; variables = default_variables(mod0))
    quotes = collect(quotes)
    parameters = __UnitBox(AccessibleModel(mod0, variables))
    x0 = [0 < x < 1 ? x : oftype(x, 0.5) for x in AccessibleModels.transformed_vec(parameters.variables)]
    who = "fit($(nameof(typeof(mod0))), Fit.LeastSquares)"
    return parameters(__least_squares(parameters, x0, quotes, map(__primal_quote, quotes), method, who))
end

# A spline: one knot at each quote maturity.
function fit(mod0::T, quotes, method::Fit.LeastSquares; extrapolation = :flat_forward) where {T <: Spline.SplineCurve}
    quotes, primal_quotes = __calibration_quotes(quotes)
    tenors = sort!([maturity(q) for q in primal_quotes])
    return __fit_knot_rates(mod0, tenors, quotes, primal_quotes, method, extrapolation)
end

# A knot curve: its knots, wherever they are, its interpolation method and its extrapolation policy.
function fit(c::Yield.AbstractInterpolatedZeroCurve, quotes, method::Fit.LeastSquares)
    quotes, primal_quotes = __calibration_quotes(quotes)
    return __fit_knot_rates(c.spline, Yield.knot_tenors(c), quotes, primal_quotes, method, c.extrapolation)
end

# The knot grid is validated once, up front, as in a loss fit, from the same starting rates.
function __fit_knot_rates(spline, tenors, quotes, primal_quotes, method::Fit.LeastSquares, extrapolation)
    who = "fit($(spline), Fit.LeastSquares)"
    grid0 = Yield.KnotGrid(__knot_fit_seed(spline, tenors), tenors, spline; who)
    parameters = __KnotRates(spline, Yield.__interpolation(spline, length(grid0.tenors)), grid0.tenors, extrapolation)
    z = __least_squares(parameters, grid0.rates, quotes, primal_quotes, method, who)
    # The public result is validated; knot rates that carry derivatives keep the caller's policy.
    eltype(z) <: ForwardDiff.Dual || return Yield.ZeroRateCurve(z, grid0.tenors, spline; extrapolation = __primal_extrapolation(extrapolation))
    return Yield.__build_public(
        spline, Yield.KnotGrid(z, grid0.tenors, spline; who = "ZeroRateCurve"), parameters.interpolation; extrapolation
    )
end
