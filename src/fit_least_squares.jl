# ── Least-squares fits ─────────────────────────────────────────────────────────────────
#
# `Fit.LeastSquares` fits a model to quotes with one solve for every model. A fit has a parameter
# vector θ and a rule that builds the model from it: `__KnotRates` (a knot curve from its knot rates),
# `__ForeignKnotRates` (an FX model from its foreign curve's) or `__Variables` (any other model from
# its variables). Everything else sees only the residuals
#
#     rᵢ(θ) = √wᵢ · (value(model(θ), quoteᵢ) − priceᵢ)
#
# and never the model or the contracts:
#
# 1. Levenberg–Marquardt on the residuals' ForwardDiff Jacobian (first derivatives only) takes the
#    parameters to where a Gauss–Newton step is at rounding, or no longer shrinks.
# 2. Newton's method with the exact Hessian of the sum of squares goes on to a point it can certify
#    as a minimum, which is the point returned: the Hessian there is positive definite, and the
#    Newton correction there is negligible against the size of the prices.
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
__bounds(::__KnotRates) = nothing
__parameter_name(p::__KnotRates, k) = "the knot rate at $(p.tenors[k])"
# What to say when a fit of an interpolation with formula switches fails (nothing for the others).
__switch_hint(p::__KnotRates) = p.spline isa Union{Spline.MonotoneConvex, Spline.PCHIP, Spline.Akima} ?
    "$(p.spline) switches formulas as its knot rates change, and a minimum on a switch cannot be reached by " *
    "derivatives (on one, they may not even be finite). Use Spline.Linear() or Spline.Cubic(), or a loss fit " *
    "with an optimizer that does not use derivatives." : ""
__on_kink(p::__KnotRates, z, δ) = Yield.__near_kink(p(z), z, δ)

# θ = the knot rates of the foreign curve of an `FX.Forwards`, which is the model priced: for quotes
# that have no implied quote on the foreign curve (see the `FX.Forwards` method of `fit`).
struct __ForeignKnotRates{P, S, D, K}
    pair::P
    spot::S
    domestic::D
    knots::K
end
(p::__ForeignKnotRates)(z) = FX.Forwards(p.pair, p.spot, p.domestic, p.knots(z))
__map_parameters(f::F, p::__ForeignKnotRates) where {F} =
    __ForeignKnotRates(p.pair, f(p.spot), p.domestic, __map_parameters(f, p.knots))
__parameter_duals(p::__ForeignKnotRates) = (p.spot, __parameter_duals(p.knots)...)
__bounds(::__ForeignKnotRates) = nothing
__parameter_name(p::__ForeignKnotRates, k) = __parameter_name(p.knots, k)
__switch_hint(p::__ForeignKnotRates) = __switch_hint(p.knots)
__on_kink(p::__ForeignKnotRates, z, δ) = __on_kink(p.knots, z, δ)

# θ = the model's variables, in their own units, strictly inside their intervals `lower .. upper`.
# (Not AccessibleModels' unit box, which the loss fit uses: a tolerance there is a fraction of the
# interval, and a wide interval leaves a coordinate too few digits for the variable.)
struct __Variables{P, B}
    model::P   # AccessibleModels' description of the variables: their optics and intervals
    lower::B
    upper::B
end
(p::__Variables)(x) = AccessibleModels.from_raw(x, p.model)
__map_parameters(f::F, p::__Variables) where {F} = p
__parameter_duals(::__Variables) = ()
__bounds(p::__Variables) = (p.lower, p.upper)
__parameter_name(::__Variables, k) = "variable $k"
__switch_hint(::__Variables) = ""
__on_kink(::__Variables, x, δ) = false

const __LEAST_SQUARES_HINT = "Fit.LeastSquares differentiates quote prices, the coupon rates of Bond.Fixed " *
    "and Bond.Floating, the amounts of Cashflow and FX.BasisSwapLeg, the strike of FX.Forward (also inside a " *
    "Composite or a vector), its weights, a Yield.FlatForwardAt forward, and the spot of FX.Forwards; a dual " *
    "number in the model's other fields or inside another contract type is not supported."

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

# Every tolerance is relative to the size of the prices, which the parameters' units, the quotes'
# notionals and the weights cannot change: the norm of the weighted prices and of the weighted
# residuals `f`, or (for quotes that are all priced at zero and fitted exactly) of the price changes
# `D .* θ` that the parameters themselves stand for. `D` measures each parameter in prices: its
# column norm in the residuals' Jacobian, or the square root of its diagonal in the Hessian.
__price_size(prices, f, D, θ) = max(prices + norm(f), norm(D .* θ))

const __LSQ_FLAT = 1.0e-12       # a relative change of the sum of squares that is rounding
const __LSQ_STEP = 1.0e-12       # the search ends at a Gauss–Newton step that moves the prices by this much of their size
const __LSQ_PATIENCE = 10        # or after this many steps without a clearly smaller one, once the sum of squares is flat
const __LSQ_NEWTON = 1.0e-10     # a certified minimum: no parameter's Newton correction moves the prices by more
const __LSQ_NEWTON_STEPS = 10

__undetermined(who, parameters, k) = ArgumentError(
    "$who: $(__parameter_name(parameters, k)) does not move any quote's price (with a nonzero weight), " *
        "so the quotes do not determine it. Remove it from the fit, or add quotes that depend on it."
)

# Levenberg–Marquardt from `θ`, with Marquardt's scaling by the largest column norms seen and
# Nielsen's update of the damping `μ`. It ends where the undamped Gauss–Newton step moves the prices
# by at most `__LSQ_STEP` of their size (`prices` is the norm of the weighted prices), or has stopped
# shrinking, and returns that point without judging it: `__certify_minimum` does.
function __levenberg_marquardt(r::R, θ::Vector{T}, cfg::C, bounds, prices, maxiters::Int) where {R, T, C}
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
    # No parameter moves a residual here, so there is no step to take. That is a minimum when the
    # residuals' curvature alone makes the Hessian positive definite, which is the certificate's call.
    iszero(maximum(D)) && return (; θ, step = zero(T), iters = 0)
    # A column that is zero here (a Nelson–Siegel decay time with no slope or hump yet) is damped at
    # rounding of the largest until it moves.
    widen!(D, J) = (c = __column_norms(J); D .= max.(D, c, eps(T) * maximum(c)); D)
    widen!(D, J)
    gauss_newton(θ, f, J) = norm(D .* (qr(J, ColumnNorm()) \ f)) / __price_size(prices, f, D, θ)
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
        bounds === nothing || (δ .= clamp.(δ, 0.9 .* (bounds[1] .- θ), 0.9 .* (bounds[2] .- θ)))
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
            # (A step that shrinks by a millionth each time is not converging: a smaller one counts
            # when it is smaller by a tenth.)
            if step < 0.9 * least
                least, waited = step, 0
            else
                waited += 1
            end
            (flat && waited >= __LSQ_PATIENCE) && break
        else
            μ *= ν
            ν *= 2
            μ > 1.0e16 && break
        end
    end
    return (; θ, step, iters)
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

function __hessian(e, w, θ, outer, inner)
    gradient(x) = __least_squares_gradient(e, w, x, inner)
    H = ForwardDiff.jacobian(gradient, θ, ForwardDiff.JacobianConfig(gradient, θ, ForwardDiff.Chunk(θ), outer), Val(false))
    return (H .+ H') ./ 2
end

# The Cholesky factor `C` of H equilibrated by `d = √diag(H)` (so that no check depends on the
# parameters' units); or, where H is not positive definite and well conditioned, the reason in words.
# The condition number is the 1-norm's, which needs no singular values and so works in any float type.
function __hessian_factors(H)
    h = diag(H)
    indefinite = "the Hessian of the sum of squares is not positive definite"
    (all(isfinite, H) && all(>(0), h)) || return indefinite
    d = sqrt.(h)
    Hs = H ./ (d .* d')
    C = cholesky(Symmetric(Hs); check = false)
    issuccess(C) || return indefinite
    κ = opnorm(Hs, 1) * opnorm(inv(C), 1)
    κ <= 1.0e10 || return "the Hessian of the sum of squares is ill-conditioned (condition number $κ)"
    return C, d
end

# From where the search ended, Newton's method on g = 0 with the exact Hessian, until it stands on a
# point it can certify. That point is returned, with the Newton correction `δ` left there and the
# factors `C`, `d` of its Hessian; the correction is not applied, so that everything certified is
# true of the parameters returned. A certified point has:
#
# - a Hessian that is positive definite and well conditioned;
# - a Newton correction that is negligible: for each parameter, `d` times its correction (the price
#   change it stands for) is at most `__LSQ_NEWTON` of the size of the prices. A rate that matters in
#   its eleventh decimal, because its cashflows are 1e11 years away, is measured as finely;
# - every bounded parameter farther from the ends of its interval than a correction of that size.
#
# Steps on the way may not raise the sum of squares.
function __certify_minimum(e, w, prices, searched, outer, inner, maxiters)
    parameters, who = e.parameters, e.who
    θ = searched.θ
    stopped = searched.step > __LSQ_STEP && searched.iters >= maxiters
    failed(why) = FitConvergenceError(
        stopped ? Optimization.SciMLBase.ReturnCode.MaxIters : Optimization.SciMLBase.ReturnCode.Stalled,
        "$who " * (stopped ? "did not converge within $maxiters iterations" : "did not reach a minimum") *
            ": $why (after $(searched.iters) Levenberg–Marquardt iterations, where a Gauss–Newton step moved the " *
            "prices by $(searched.step) of their size). " * __failure_hint(parameters, θ)
    )
    stopped && throw(failed("the search was still moving"))
    sum_of_squares(x) = (f = e(x); w === nothing ? sum(abs2, f) : sum(w .* abs2.(f)))
    bounds = __bounds(parameters)
    margin(x, k) = bounds === nothing ? oftype(x[k], Inf) : min(x[k] - bounds[1][k], bounds[2][k] - x[k])
    ϕ = sum_of_squares(θ)
    moved = zero(ϕ)
    for steps in 0:__LSQ_NEWTON_STEPS
        H = __hessian(e, w, θ, outer, inner)
        # Neither the prices nor their curvature depend on such a parameter.
        k = findfirst(iszero, diag(H))
        k === nothing || throw(__undetermined(who, parameters, k))
        factors = __hessian_factors(H)
        factors isa String && throw(failed(factors))
        C, d = factors
        δ = -(C \ (__least_squares_gradient(e, w, θ, inner) ./ d)) ./ d
        size = max(prices + sqrt(ϕ), norm(d .* θ))
        moved = maximum(d .* abs.(δ)) / size
        if all(d .* abs.(δ) .<= __LSQ_NEWTON * size)
            # A minimum this near a bound cannot be told from one on it.
            k = findfirst(k -> margin(θ, k) <= __LSQ_NEWTON * size / d[k], eachindex(θ))
            k === nothing || throw(failed("$(__parameter_name(parameters, k)) is within the fit's resolution of an end of its interval"))
            return θ, δ, C, d
        end
        steps == __LSQ_NEWTON_STEPS && break
        all(k -> margin(θ .+ δ, k) > 0, eachindex(θ)) || throw(failed("a Newton step leaves a variable's interval"))
        # A step may not raise the sum of squares beyond rounding. It is halved until it does not:
        # where the curvature changes quickly (near a fold), a full step overshoots.
        θt = θ .+ δ
        ϕt = sum_of_squares(θt)
        halvings = 0
        while !(ϕt <= ϕ * (1 + __LSQ_FLAT))
            halvings == 30 && throw(failed("a Newton step with the full Hessian raises the sum of squares"))
            halvings += 1
            δ = δ ./ 2
            θt = θ .+ δ
            ϕt = sum_of_squares(θt)
        end
        θ, ϕ = θt, ϕt
    end
    throw(failed("after $(__LSQ_NEWTON_STEPS) Newton steps with the full Hessian, a correction still moves the prices by $moved of their size"))
end

__failure_hint(p::__KnotRates, z) = isempty(__switch_hint(p)) ?
    "With fewer quotes than knots, or knots that move the quotes' prices only together, the quotes do not " *
    "determine the knot rates: use fewer knots, or knots where the quotes' cashflows are. Otherwise check " *
    "the quotes for prices that no curve with this interpolation comes near." : __switch_hint(p)
__failure_hint(p::__ForeignKnotRates, z) = __failure_hint(p.knots, z)
function __failure_hint(p::__Variables, x)
    width = p.upper .- p.lower
    k = argmin(min.(x .- p.lower, p.upper .- x) ./ width)
    below, above = x[k] - p.lower[k], p.upper[k] - x[k]
    min(below, above) <= 1.0e-6 * width[k] && return "Variable $k is at the $(below < above ? "lower" : "upper") end of " *
        "its interval: the quotes pull it to or beyond that bound. Widen the interval in `variables`."
    return "With fewer quotes than variables, or variables that move the quotes' prices only together, the quotes " *
        "do not determine them. Otherwise check the quotes for prices that the model cannot come near, or start " *
        "from a model nearer the fit."
end

struct __Retag{T} end
(::__Retag{T})(x) where {T} = x
(::__Retag{T})(x::ForwardDiff.Dual{<:Any, V, K}) where {T, V, K} =
    ForwardDiff.Dual{T, V, K}(ForwardDiff.value(x), ForwardDiff.partials(x))

# The certified minimum `θ` (its Newton correction `δ`, its Hessian's factors `C` and `d`) with the
# derivatives of the calibration with respect to the caller's dual numbers in `quotes`, `weights` and
# `parameters`; `θ` itself when there are none.
function __implicit_least_squares(outer::TO, inner, ::Type{Own}, parameters, primal_parameters, θ, δ, C, d, quotes, weights, who) where {TO, Own}
    # The caller's dual type, from the pricing errors at the fit with their numbers in plain arithmetic.
    model = parameters(θ)
    errors = [__quote_value(model, q) - q.price for q in quotes]
    D = __common_dual_type(Iterators.flatten((errors, something(weights, ()), __parameter_duals(parameters))))
    D === nothing && return θ
    # On a kink of the interpolant the refit has no derivative; nor may the minimum, a correction
    # away, lie across one from the point whose Hessian this is (see `__implicit_knot_curve`).
    __on_kink(primal_parameters, θ, .-δ) && throw(__fit_kink_error(primal_parameters(θ)))
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

# Whether a weight is part of the fit: a quote with a weight of zero is left out (it is not priced,
# converted, or given a knot), unless that zero carries derivatives, to which the minimum responds.
__has_partials(x) = false
__has_partials(x::ForwardDiff.Dual) = !iszero(ForwardDiff.partials(x)) || __has_partials(ForwardDiff.value(x))

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
    outer = ForwardDiff.Tag(__FitTagKey{O, :outer}(), eltype(θ0))
    inner = ForwardDiff.Tag(__FitTagKey{O, :inner}(), eltype(θ0))  # created second, so it nests inside `outer`
    Own = Tuple{typeof(outer), typeof(inner)}
    primal_parameters = __map_parameters(__primal, parameters)
    primal_weights = __primal_weights(weights, length(quotes), who)
    # A quote with a weight of zero is left out: it is not priced. Its row returns for the derivatives
    # only if its weight carries them, since the minimum does respond to that weight.
    if primal_weights !== nothing && !all(>(0), primal_weights)
        active = findall(>(0), primal_weights)
        kept = findall(i -> primal_weights[i] > 0 || __has_partials(weights[i]), eachindex(primal_weights))
        quotes, weights = quotes[kept], [weights[i] for i in kept]
        kept_primal_quotes = primal_quotes[kept]
        primal_quotes, primal_weights = primal_quotes[active], primal_weights[active]
    else
        kept_primal_quotes = primal_quotes
    end
    scale = primal_weights === nothing ? nothing : map(sqrt, primal_weights)
    prices = norm(scale === nothing ? [float(q.price) for q in primal_quotes] : [s * q.price for (s, q) in zip(scale, primal_quotes)])
    r = __Residuals(primal_parameters, primal_quotes, scale, Own, who)
    cfg = ForwardDiff.JacobianConfig(r, θ0, ForwardDiff.Chunk(θ0), inner)
    searched = __levenberg_marquardt(r, θ0, cfg, __bounds(parameters), prices, maxiters)
    e = __Residuals(primal_parameters, primal_quotes, nothing, Own, who)
    θ, δ, C, d = __certify_minimum(e, primal_weights, prices, searched, outer, inner, maxiters)
    differentiated = parameters !== primal_parameters ||
        any(((q, p),) -> q !== p, zip(quotes, kept_primal_quotes)) ||
        (weights !== nothing && any(w -> __ad_depth(w) > 0, weights))
    differentiated || return θ
    return __implicit_least_squares(outer, inner, Own, parameters, primal_parameters, θ, δ, C, d, quotes, weights, who)
end

# The weights without their dual numbers, validated: one per quote, finite and non-negative. (Rows
# with a weight of zero are skipped, which would drop a negative or NaN weight without a word; and
# `map` over quotes and weights of different lengths would drop the extra ones.)
__primal_weights(::Nothing, n, who) = nothing
function __primal_weights(weights, n, who)
    length(weights) == n || throw(DimensionMismatch("$who: Fit.LeastSquares has $(length(weights)) weights for $n quotes"))
    w = map(__primal, weights)
    all(x -> isfinite(x) && x >= 0, w) || throw(ArgumentError("$who: Fit.LeastSquares weights must be finite and non-negative"))
    return w
end

# ── `fit` methods ──────────────────────────────────────────────────────────────────────

# Any model with `variables`: their own values, strictly inside their intervals, from the model's (a
# variable on a bound starts from the middle of its interval, as in a loss fit). The solve works in
# Float64 at least, whatever the model's own number type: its tolerances assume that precision.
function fit(mod0, quotes, method::Fit.LeastSquares; variables = default_variables(mod0))
    quotes = collect(quotes)
    who = "fit($(nameof(typeof(mod0))), Fit.LeastSquares)"
    model = AccessibleModel(mod0, variables)
    # A variable given a set of values rather than an interval would be fitted between them.
    k = findfirst(d -> !(d isa AccessibleModels.Distributions.ContinuousUnivariateDistribution), model.distributions)
    k === nothing || throw(ArgumentError("$who: variable $k is not given an interval (lower .. upper) to vary in, which Fit.LeastSquares requires"))
    at(u) = AccessibleModels.itransform(map(_ -> u, model.distributions), model)
    lower, middle, upper = at(0.0), at(0.5), at(1.0)
    raw = map(__primal, AccessibleModels.raw_vec(model))
    T = promote_type(Float64, eltype(raw), eltype(lower), eltype(upper))
    x0 = T[lower[k] < raw[k] < upper[k] ? raw[k] : middle[k] for k in eachindex(raw)]
    parameters = __Variables(model, collect(T, lower), collect(T, upper))
    return parameters(__least_squares(parameters, x0, quotes, map(__primal_quote, quotes), method, who))
end

# A spline: one knot at the maturity of each quote that is part of the fit.
fit(mod0::T, quotes, method::Fit.LeastSquares; extrapolation = :flat_forward) where {T <: Spline.SplineCurve} =
    __fit_knots(mod0, quotes, method, identity; extrapolation)

# A knot curve: its knots, wherever they are, its interpolation method and its extrapolation policy.
fit(c::Yield.AbstractInterpolatedZeroCurve, quotes, method::Fit.LeastSquares) = __fit_knots(c, quotes, method, identity)

# `model` makes the fit's parameters from the curve's `__KnotRates`: themselves, or an FX model's.
function __fit_knots(spline::Spline.SplineCurve, quotes, method, model::M; extrapolation = :flat_forward) where {M}
    quotes, primal_quotes = __calibration_quotes(quotes)
    tenors = sort!([maturity(q) for q in __knot_quotes(primal_quotes, method.weights, "fit($(spline), Fit.LeastSquares)")])
    return __fit_knot_rates(spline, tenors, quotes, primal_quotes, method, extrapolation, model)
end
function __fit_knots(c::Yield.AbstractInterpolatedZeroCurve, quotes, method, model::M) where {M}
    quotes, primal_quotes = __calibration_quotes(quotes)
    return __fit_knot_rates(c.spline, Yield.knot_tenors(c), quotes, primal_quotes, method, c.extrapolation, model)
end

# The quotes that get a knot: those that are part of the fit. A weight of zero that carries
# derivatives has none here: any positive weight would add a knot at its quote's maturity, so the
# fitted curve does not change smoothly with it.
__knot_quotes(quotes, ::Nothing, who) = quotes
function __knot_quotes(quotes, weights, who)
    w = __primal_weights(weights, length(quotes), who)
    any(i -> iszero(w[i]) && __has_partials(weights[i]), eachindex(w)) && throw(
        ArgumentError(
            "$who cannot differentiate with respect to a weight of zero: a positive weight would add a knot at " *
                "its quote's maturity. Fit a curve with knots of its own: " *
                "fit(ZeroRateCurve(rates, knots, spline), quotes, Fit.LeastSquares(weights = w))."
        )
    )
    return [q for (q, x) in zip(quotes, w) if x > 0]
end

# The knot grid is validated once, up front, as in a loss fit, from the same starting rates.
function __fit_knot_rates(spline, tenors, quotes, primal_quotes, method::Fit.LeastSquares, extrapolation, model::M) where {M}
    who = "fit($(spline), Fit.LeastSquares)"
    grid0 = Yield.KnotGrid(__knot_fit_seed(spline, tenors), tenors, spline; who)
    knots = __KnotRates(spline, Yield.__interpolation(spline, length(grid0.tenors)), grid0.tenors, extrapolation)
    z = __least_squares(model(knots), collect(float.(grid0.rates)), quotes, primal_quotes, method, who)
    # The public result is validated; knot rates that carry derivatives keep the caller's policy.
    eltype(z) <: ForwardDiff.Dual || return Yield.ZeroRateCurve(z, grid0.tenors, spline; extrapolation = __primal_extrapolation(extrapolation))
    return Yield.__build_public(
        spline, Yield.KnotGrid(z, grid0.tenors, spline; who = "ZeroRateCurve"), knots.interpolation; extrapolation
    )
end

# `FX.Forwards` with a spline or a knot curve as the foreign curve. Quotes left out by their weights
# are dropped first, so that they are not converted either.
#
# Outrights and basis-swap legs have implied quotes on the foreign curve, and the curve is fitted to
# those like the other methods' (`__fit_fx_via_implied`), which keeps the derivatives with respect to
# the domestic curve. An implied quote's residual is its own quote's pricing residual divided by
# `FX.__implied_residual_scale`, so each weight is multiplied by that scale squared: the fit then
# minimizes the weighted pricing residuals of the quotes as given, whatever mix they are.
#
# Any other contract (a `Composite` of forwards, for one) is priced under the FX model itself, whose
# foreign knot rates are the parameters (`__ForeignKnotRates`).
function fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.LeastSquares; kwargs...) where {P, S, D, F <: Union{Spline.SplineCurve, Yield.AbstractInterpolatedZeroCurve}}
    quotes = collect(quotes)
    weights = method.weights
    if weights !== nothing
        w = __primal_weights(weights, length(quotes), "fit(FX.Forwards, Fit.LeastSquares)")
        rows = findall(i -> w[i] > 0 || __has_partials(weights[i]), eachindex(w))
        quotes, weights = quotes[rows], [weights[i] for i in rows]
    end
    if all(q -> q isa FinanceCore.Quote{<:Any, <:Union{FX.Forward, FX.BasisSwapLeg}}, quotes)
        scales = [FX.__implied_residual_scale(mod0, q) for q in quotes]
        weights = weights === nothing ? scales .^ 2 : [x * s^2 for (x, s) in zip(weights, scales)]
        return __fit_fx_via_implied(mod0, quotes, Fit.LeastSquares(; weights, maxiters = method.maxiters); kwargs...)
    end
    model(knots) = __ForeignKnotRates(mod0.pair, mod0.spot, mod0.domestic, knots)
    foreign = __fit_knots(mod0.foreign, quotes, Fit.LeastSquares(; weights, maxiters = method.maxiters), model; kwargs...)
    return @set mod0.foreign = foreign
end
