# ── Least-squares fits ─────────────────────────────────────────────────────────────────
#
# `Fit.LeastSquares` fits a model to quotes with one solve for every model. A fit has a parameter
# vector θ and a rule that builds the model from it: `__KnotRates` (a knot curve from its knot rates),
# `__ForeignKnotRates` (an FX model from its foreign curve's) or `__Variables` (any other model from
# its variables). The search and the certificate use the same weighted pricing objective with residuals
#
#     rᵢ(θ) = √wᵢ · (value(model(θ), quoteᵢ) − priceᵢ)
#
# The model and contract also determine the price resolution used in the rounding checks.
#
# 1. Levenberg–Marquardt on the residuals' ForwardDiff Jacobian (first derivatives only) takes the
#    parameters to where a Gauss–Newton step is at rounding, or no longer shrinks.
# 2. Newton's method with the exact Hessian of the sum of squares goes on to a point it can certify
#    as a minimum, which is the point returned: the Hessian there is positive definite, and the
#    Newton correction meets the price and slope tolerances for each quote.
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
function __check_fit_kink(p::__KnotRates, z, δ)
    curve = p(z)
    Yield.__near_kink(curve, z, δ) && throw(__fit_kink_error(curve))
    return nothing
end

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
__check_fit_kink(p::__ForeignKnotRates, z, δ) = __check_fit_kink(p.knots, z, δ)

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
__check_fit_kink(::__Variables, x, δ) = nothing

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

__quote_size(price, residual) = max(abs(price), abs(price + residual))
__cancellation_allowance(model, q, price) = zero(price)
__price_spacing(model, q, price, size) = max(4 * eps(size), __cancellation_allowance(model, q, price))

# `__price_move` measures the largest first-order change in any one quote's price from any one
# parameter's step (`J[i, k] * δ[k]`), as a fraction of that quote's size.
# A quote that a parameter does not move says nothing about that
# parameter, however large it is, and neither the parameters' units nor the notionals matter.
#
# A quote's size is the larger of its price and its value now (`f` holds the pricing errors), which
# is also the scale of its error's rounding. A quote priced at zero has its value under the starting
# parameters as a `floor`, as a bootstrap step does; one that has no size at all cannot be measured,
# and a parameter that moves it is never found to have converged.
function __price_move(J, δ, sizes)
    worst = zero(float(real(eltype(J))))
    for i in axes(J, 1)
        for k in axes(J, 2)
            move = abs(J[i, k] * δ[k])
            iszero(move) || (worst = max(worst, move / sizes[i]))
        end
    end
    return worst
end

const __LSQ_FLAT = 1.0e-12       # a relative change of the sum of squares that is rounding
const __LSQ_STEP = 1.0e-12       # the search ends at a Gauss–Newton step that moves no price by more than this of its size
const __LSQ_PATIENCE = 10        # or after this many steps without a clearly smaller one, once the sum of squares is flat
const __LSQ_NEWTON = 1.0e-10     # a certified minimum: its Newton correction moves no price by more than this of its size
const __LSQ_SLOPE = 1.0e-8       # maximum weighted change of a pricing gradient in the inverse Hessian metric
const __LSQ_NEWTON_STEPS = 10

__undetermined(who, parameters, k) = ArgumentError(
    "$who: $(__parameter_name(parameters, k)) does not move any quote's price (with a nonzero weight), " *
        "so the quotes do not determine it. Remove it from the fit, or add quotes that depend on it."
)

# Compare the predicted and observed reductions before summing. An immovable quote can
# contribute an arbitrarily large constant to the loss and must cancel from both reductions
# and from the scale used to decide whether a trial step improved it.
function __lm_reductions(f, ft, linear_step)
    predicted = actual = level = zero(eltype(f))
    for i in eachindex(f)
        change = f[i] - ft[i]
        linear = linear_step[i]
        iszero(change) && iszero(linear) && continue
        predicted -= linear * (2 * f[i] + linear)
        actual += change * (f[i] + ft[i])
        level += abs2(f[i])
    end
    return predicted, actual, level
end

# Levenberg–Marquardt from `θ`, with Marquardt's scaling by the largest column norms seen and
# Nielsen's update of the damping `μ`. It ends where the undamped Gauss–Newton step moves no price by
# more than `__LSQ_STEP` of its size (`__price_move`; `prices` and `floor` are weighted as the
# residuals are), or has stopped shrinking, and returns that point without judging it:
# `__certify_minimum` does.
function __levenberg_marquardt(r::R, θ::Vector{T}, cfg::C, bounds, prices, floor, maxiters::Int) where {R, T, C}
    n = length(θ)
    f = r(θ)
    J = ForwardDiff.jacobian(r, θ, cfg, Val(false))
    nonfinite(at) = FitConvergenceError(
        Optimization.SciMLBase.ReturnCode.Failure,
        String(rstrip("$(r.who): the pricing residuals or their Jacobian are not finite $at. $(__switch_hint(r.parameters))"))
    )
    (all(isfinite, f) && all(isfinite, J)) || throw(nonfinite("at the starting parameters"))
    m = length(f)
    D = __column_norms(J)
    # No parameter moves a residual here, so there is no step to take. That is a minimum when the
    # residuals' curvature alone makes the Hessian positive definite, which is the certificate's call.
    iszero(maximum(D)) && return (; θ, step = zero(T), iters = 0)
    # A column that is zero here (a Nelson–Siegel decay time with no slope or hump yet) is damped at
    # rounding of the largest until it moves.
    widen!(D, J) = (c = __column_norms(J); D .= max.(D, c, eps(T) * maximum(c)); D)
    widen!(D, J)
    gauss_newton(f, J) = __price_move(J, qr(J, ColumnNorm()) \ f, max.(__quote_size.(prices, f), floor))
    step = least = gauss_newton(f, J)
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
        θt = θ .+ δ
        ft = r(θt)
        predicted, actual, level = __lm_reductions(f, ft, J * δ)
        flat = predicted <= __LSQ_FLAT * level
        ρ = actual / predicted
        if (predicted > 0 && ρ > 1.0e-4) || (flat && actual >= -__LSQ_FLAT * level)
            θ, f = θt, ft
            J = ForwardDiff.jacobian(r, θ, cfg, Val(false))
            all(isfinite, J) || throw(nonfinite("after $iters iterations"))
            widen!(D, J)
            μ = flat ? μ / 3 : μ * max(one(T) / 3, 1 - (2ρ - 1)^3)
            ν = T(2)
            step = gauss_newton(f, J)
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

function __pricing_jacobian(e, x, inner)
    cfg = ForwardDiff.JacobianConfig(e, x, ForwardDiff.Chunk(x), inner)
    return ForwardDiff.jacobian(e, x, cfg, Val(false))
end
__least_squares_gradient(J, w, f) = J' * (w === nothing ? f : w .* f)

# H = J'WJ + Σ wᵢ eᵢ ∇²PVᵢ. At a turning point, rounding eᵢ to zero can make a
# singular minimum appear differentiable. R bounds the resulting uncertainty in each
# Hessian entry by Σ noiseᵢ |∂²PVᵢ/∂θⱼ∂θₖ|, with noiseᵢ = weight × price resolution.
# Each outer lane seeds one parameter. Its partials give a column of H from the
# gradient and the matching column of R from J, without storing a Hessian per quote.
function __hessian(e, w, θ, outer, inner, noise, ::ForwardDiff.Chunk{N} = ForwardDiff.Chunk(θ)) where {N}
    n = length(θ)
    H = similar(θ, n, n)
    R = similar(H)
    for firstcol in 1:N:n
        x = map(eachindex(θ)) do i
            v = θ[i]
            seed = ntuple(lane -> i == firstcol + lane - 1 ? one(v) : zero(v), Val(N))
            ForwardDiff.Dual{typeof(outer)}(v, seed)
        end
        J = __pricing_jacobian(e, x, inner)
        g = __least_squares_gradient(J, w, e(x))
        for k in firstcol:min(firstcol + N - 1, n)
            lane = k - firstcol + 1
            for j in axes(H, 1)
                H[j, k] = ForwardDiff.partials(g[j], lane)
                R[j, k] = sum(i -> noise[i] * abs(ForwardDiff.partials(J[i, j], lane)), axes(J, 1); init = zero(eltype(R)))
            end
        end
    end
    return (H .+ H') ./ 2, R
end

# If H = D U'U D, a symmetric error E with |Eⱼₖ| ≤ Rⱼₖ preserves positive
# definiteness whenever ‖U⁻ᵀ D⁻¹ E D⁻¹ U⁻¹‖₂ < 1. The 1-norm of the symmetric
# entrywise bound is a conservative upper bound, shared by all sensitivity directions.
function __resolved_curvature(R, C, d)
    A = abs.(inv(C.U))
    bound = A' * (max.(R, R') ./ (d .* d')) * A
    return opnorm(bound, 1) < 1
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

# Compare changes in individual squared residuals before summing. Constant cashflows then
# cancel exactly, even when their value dwarfs every instrument the model can move.
function __acceptable_step(e, model, w, prices, f, ft)
    change = slack = zero(eltype(ft))
    for i in eachindex(ft)
        delta = ft[i] - f[i]
        iszero(delta) && continue
        weight = w === nothing ? 1 : w[i]
        p = prices[i]
        size = max(__quote_size(p, f[i]), abs(p + ft[i]))
        u = __price_spacing(model, e.quotes[i], p, size)
        change += weight * delta * (ft[i] + f[i])
        slack += weight * (2u * (abs(ft[i]) + abs(f[i])) + 2u^2)
    end
    return isfinite(change) && change <= slack
end

# A bound is indistinguishable from the fitted point when both the price and its
# local slope variation over the remaining distance are within resolution.
function __unresolved_bound(bounds, θ, J, sizes, jacobian, outer, slope_norm)
    bounds === nothing && return nothing
    lower, upper = bounds
    for k in eachindex(θ)
        m = min(θ[k] - lower[k], upper[k] - θ[k])
        m <= 0 && return k
        isfinite(m) || continue
        __price_move(view(J, :, k:k), m, sizes) <= __LSQ_NEWTON || continue
        # A distant bound can itself be a stationary point. Measure the variation at θ.
        x = [ForwardDiff.Dual{typeof(outer)}(v, i == k ? m : zero(m)) for (i, v) in pairs(θ)]
        slopes = map(v -> ForwardDiff.partials(v, 1), jacobian(x))
        slope_norm(slopes) <= __LSQ_SLOPE && return k
    end
    return nothing
end

function __backtrack_newton(e, model, w, prices, θ, f, δ, θt)
    for halvings in 0:30
        ft = e(θt)
        __acceptable_step(e, model, w, prices, f, ft) && return (θt, ft, iszero(halvings))
        halvings == 30 && break
        δ = δ ./ 2
        θt = θ .+ δ
    end
    return nothing
end

# From where the search ended, Newton's method on g = 0 with the exact Hessian, until it stands on a
# point it can certify. That point is returned, with the Newton correction `δ` left there and the
# factors `C`, `d` of its Hessian; the correction is not applied, so that everything certified is
# true of the parameters returned. A certified point has:
#
# - a Hessian that is positive definite and well conditioned;
# - a Newton correction that moves no quote's price by more than `__LSQ_NEWTON` of that quote's size
#   (`__price_move`), parameter by parameter. A rate that matters in its eleventh decimal, because
#   its cashflows are 1e11 years away, is measured as finely; a parameter that the prices are flat in
#   where it stands is measured as loosely as they are flat, which is why there is the next test;
# - stable pricing slopes across that correction, measured in the Hessian metric. When this
#   check calls for refinement, a full Newton step must also have stopped improving the correction;
# - every bounded parameter strictly inside its interval, and farther from its ends than a move the
#   price test would not notice.
#
# The correction is a local estimate of the distance to the minimum, not a bound on it. Steps on the
# way may not raise the sum of squares.
function __certify_minimum(e, w, prices, floor, searched, outer, inner, maxiters)
    parameters, who = e.parameters, e.who
    θ = searched.θ
    stopped = searched.step > __LSQ_STEP && searched.iters >= maxiters
    failed(why, at) = FitConvergenceError(
        stopped ? Optimization.SciMLBase.ReturnCode.MaxIters : Optimization.SciMLBase.ReturnCode.Stalled,
        "$who " * (stopped ? "did not converge within $maxiters iterations" : "did not reach a minimum") *
            ": $why (after $(searched.iters) Levenberg–Marquardt iterations, where a Gauss–Newton step moved a " *
            "price by $(searched.step) of its size). " * __failure_hint(parameters, at)
    )
    stopped && throw(failed("the search was still moving", θ))
    bounds = __bounds(parameters)
    jacobian(x) = __pricing_jacobian(e, x, inner)
    f = e(θ)
    refining, previous = false, nothing
    local reason
    for steps in 0:__LSQ_NEWTON_STEPS
        model = parameters(θ)
        sizes = __quote_size.(prices, f)
        spacing = __price_spacing.(Ref(model), e.quotes, prices, sizes)
        H, uncertainty = __hessian(e, w, θ, outer, inner, w === nothing ? spacing : w .* spacing)
        k = findfirst(iszero, diag(H))
        k === nothing || throw(__undetermined(who, parameters, k))
        factors = __hessian_factors(H)
        factors isa String && throw(failed(factors, θ))
        C, d = factors

        J = jacobian(θ)
        δ = -(C \ (__least_squares_gradient(J, w, f) ./ d)) ./ d
        sizes .= max.(sizes, floor)
        moved = __price_move(J, δ, sizes)
        priced = moved <= __LSQ_NEWTON
        θt = θ .+ δ
        inside = bounds === nothing || all(k -> bounds[1][k] < θt[k] < bounds[2][k], eachindex(θt))
        if priced && inside
            slope_norm(change) = __slope_norm(change, w, C, d)
            slope = slope_norm(jacobian(θt) .- J)
            decrement = dot(δ, H, δ)
            settled = !refining || iszero(decrement) || (previous !== nothing && decrement >= previous / 4)
            if slope <= __LSQ_SLOPE && settled
                k = __unresolved_bound(bounds, θ, J, sizes, jacobian, outer, slope_norm)
                k === nothing || throw(failed("$(__parameter_name(parameters, k)) is on an end of its interval, or within the fit's resolution of one", θ))
                return (; θ, δ, C, d, uncertainty)
            end
            refining = true
        end
        if steps == __LSQ_NEWTON_STEPS
            reason = priced ? "the pricing slopes or their Newton correction have not settled" : "a correction still moves a price by $moved of its size"
            break
        end
        inside || throw(failed("a Newton step leaves a variable's interval", θ))
        trial = __backtrack_newton(e, model, w, prices, θ, f, δ, θt)
        trial === nothing && throw(failed("a Newton step with the full Hessian raises the sum of squares", θ))
        θ, f, full_step = trial
        previous = full_step ? dot(δ, H, δ) : nothing
    end
    throw(failed("after $__LSQ_NEWTON_STEPS Newton steps with the full Hessian, $reason", θ))
end

# Largest weighted change of an individual pricing gradient in the inverse Hessian metric.
# Whitening by the existing Cholesky factor makes the test independent of parameter coordinates.
function __slope_norm(change, w, C, d)
    V = C.L \ (change' ./ d)
    return maximum(i -> (w === nothing ? 1 : sqrt(w[i])) * norm(view(V, :, i)), axes(change, 1))
end

__failure_hint(p::__KnotRates, z) = isempty(__switch_hint(p)) ?
    "Fewer quotes than knots, or knots that move prices only together, can leave the quotes unable to " *
    "determine the knot rates: use fewer knots, or knots where the quotes' cashflows are. Otherwise check " *
    "the quotes for prices that no curve with this interpolation comes near." : __switch_hint(p)
__failure_hint(p::__ForeignKnotRates, z) = __failure_hint(p.knots, z)
function __failure_hint(p::__Variables, x)
    width = p.upper .- p.lower
    k = argmin(min.(x .- p.lower, p.upper .- x) ./ width)
    below, above = x[k] - p.lower[k], p.upper[k] - x[k]
    min(below, above) <= 1.0e-6 * width[k] && return "Variable $k is at the $(below < above ? "lower" : "upper") end of " *
        "its interval: the quotes pull it to or beyond that bound. Widen the interval in `variables`."
    return "Fewer quotes than variables, or variables that move prices only together, can leave the quotes " *
        "unable to determine them. Otherwise check the quotes for prices that the model cannot come near, or start " *
        "from a model nearer the fit."
end

struct __Retag{T} end
(::__Retag{T})(x) where {T} = x
(::__Retag{T})(x::ForwardDiff.Dual{<:Any, V, K}) where {T, V, K} =
    ForwardDiff.Dual{T, V, K}(ForwardDiff.value(x), ForwardDiff.partials(x))

# The certified minimum `θ` (its Newton correction `δ`, its Hessian's factors `C` and `d`) with the
# derivatives of the calibration with respect to the caller's dual numbers in `quotes`, `weights` and
# `parameters`; `θ` itself when there are none.
function __implicit_least_squares(outer::TO, inner, ::Type{Own}, parameters, primal_parameters, certificate, quotes, weights, who) where {TO, Own}
    # The caller's dual type, from the pricing errors at the fit with their numbers in plain arithmetic.
    (; θ, δ, C, d, uncertainty) = certificate
    model = parameters(θ)
    errors = (__quote_value(model, q) - q.price for q in quotes)
    D = __common_dual_type(Iterators.flatten((errors, something(weights, ()), __parameter_duals(parameters))))
    D === nothing && return θ
    __resolved_curvature(uncertainty, C, d) || throw(
        FitConvergenceError(Optimization.SciMLBase.ReturnCode.Stalled, "$who: a two-sided sensitivity cannot be resolved at the precision of the quoted prices")
    )
    # On a kink of the interpolant the refit has no derivative; nor may the minimum, a correction
    # away, lie across one from the point whose Hessian this is (see `__implicit_knot_curve`).
    __check_fit_kink(primal_parameters, θ, .-δ)
    # g at the fit with the caller's numbers under `outer`: its partials are ∂g/∂p.
    retag = __Retag{TO}()
    K = ForwardDiff.npartials(D)
    V = promote_type(ForwardDiff.valtype(D), eltype(θ))
    et = __Residuals(__map_parameters(retag, parameters), map(q -> __map_quote(retag, q), quotes), nothing, Own, who)
    wt = weights === nothing ? nothing : map(retag, weights)
    x = ForwardDiff.Dual{TO, V, K}.(θ)
    g = __least_squares_gradient(__pricing_jacobian(et, x, inner), wt, et(x))
    P = [ForwardDiff.partials(g[i], k) / d[i] for i in eachindex(g), k in 1:K]
    dθ = -(C \ P) ./ d
    W = promote_type(V, eltype(dθ))
    Dθ = ForwardDiff.Dual{ForwardDiff.tagtype(D), W, K}
    return [Dθ(W(θ[i]), ForwardDiff.Partials{K, W}(ntuple(k -> W(dθ[i, k]), K))) for i in eachindex(θ)]
end

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
    primal_weights = weights === nothing ? nothing : map(__primal, weights)
    if primal_weights !== nothing && any(iszero, primal_weights)
        active = findall(>(0), primal_weights)
        primal_quotes, primal_weights = primal_quotes[active], primal_weights[active]
    end
    scale = primal_weights === nothing ? nothing : map(sqrt, primal_weights)
    r = __Residuals(primal_parameters, primal_quotes, scale, Own, who)
    e = __Residuals(primal_parameters, primal_quotes, nothing, Own, who)
    # A zero-price quote also uses its absolute starting PV. Dividing its cancellation
    # allowance by the Newton tolerance puts both floors in units of quote size.
    # The search uses a tighter tolerance on that scale so it can resolve a one-spacing move.
    prices = [eltype(θ0)(q.price) for q in primal_quotes]
    initial = primal_parameters(θ0)
    floor = [
        max(iszero(p) ? abs(v) : zero(v), __cancellation_allowance(initial, q, p) / __LSQ_NEWTON)
            for (q, p, v) in zip(primal_quotes, prices, e(θ0))
    ]
    cfg = ForwardDiff.JacobianConfig(r, θ0, ForwardDiff.Chunk(θ0), inner)
    scaled(x) = scale === nothing ? x : scale .* x
    searched = __levenberg_marquardt(r, θ0, cfg, __bounds(parameters), scaled(prices), scaled(floor), maxiters)
    certificate = __certify_minimum(e, primal_weights, prices, floor, searched, outer, inner, maxiters)
    return __implicit_least_squares(outer, inner, Own, parameters, primal_parameters, certificate, quotes, weights, who)
end

# Drop an unseeded zero before conversion, pricing or maturity lookup. A zero with a
# nonzero partial stays in the implicit weight derivative only; the primal solve omits it.
# An inferred knot grid rejects that derivative because a positive weight would add a knot.
function __included_quotes(quotes, method, who)
    qs = collect(quotes)
    weights = method.weights
    weights === nothing && return qs, method
    length(weights) == length(qs) || throw(DimensionMismatch("$who: Fit.LeastSquares has $(length(weights)) weights for $(length(qs)) quotes"))
    w = map(__primal, weights)
    all(x -> isfinite(x) && x >= 0, w) || throw(ArgumentError("$who: Fit.LeastSquares weights must be finite and non-negative"))
    all(>(0), w) && return qs, method
    rows = findall(i -> w[i] > 0 || __has_partials(weights[i]), eachindex(w))
    return qs[rows], Fit.LeastSquares(; weights = [weights[i] for i in rows], maxiters = method.maxiters)
end

# ── `fit` methods ──────────────────────────────────────────────────────────────────────

# Any model with `variables`: their own values, strictly inside their intervals, from the model's (a
# variable on a bound starts from the middle of its interval, as in a loss fit). The solve works in
# Float64 at least, whatever the model's own number type: its tolerances assume that precision.
function fit(mod0, quotes, method::Fit.LeastSquares; variables = default_variables(mod0))
    who = "fit($(nameof(typeof(mod0))), Fit.LeastSquares)"
    quotes, method = __included_quotes(quotes, method, who)
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
function fit(spline::Spline.SplineCurve, quotes, method::Fit.LeastSquares; extrapolation = :flat_forward)
    qs, method = __included_quotes(quotes, method, "fit($(spline), Fit.LeastSquares)")
    return __fit_knots(spline, qs, method, identity; extrapolation)
end

function fit(c::Yield.AbstractInterpolatedZeroCurve, quotes, method::Fit.LeastSquares)
    qs, method = __included_quotes(quotes, method, "fit($(c.spline), Fit.LeastSquares)")
    return __fit_knots(c, qs, method, identity)
end

function __fit_knots(spline::Spline.SplineCurve, quotes, method, model::M; extrapolation = :flat_forward) where {M}
    if method.weights !== nothing && any(w -> iszero(__primal(w)), method.weights)
        throw(ArgumentError("fit($(spline), Fit.LeastSquares) cannot differentiate with respect to a weight of zero: a positive weight would add a knot at its quote's maturity. Fit a curve with knots of its own."))
    end
    quotes, primal_quotes = __calibration_quotes(quotes)
    tenors = sort!(maturity.(primal_quotes))
    return __fit_knot_rates(spline, tenors, quotes, primal_quotes, method, extrapolation, model, Float64)
end
function __fit_knots(c::Yield.AbstractInterpolatedZeroCurve, quotes, method, model::M) where {M}
    quotes, primal_quotes = __calibration_quotes(quotes)
    T = typeof(__primal(rate(first(Yield.knot_rates(c)))))
    return __fit_knot_rates(c.spline, Yield.knot_tenors(c), quotes, primal_quotes, method, c.extrapolation, model, T)
end

# The knot grid is validated once, up front, as in a loss fit, from the same starting rates.
function __fit_knot_rates(spline, tenors, quotes, primal_quotes, method::Fit.LeastSquares, extrapolation, model::M, ::Type{T}) where {M, T}
    who = "fit($(spline), Fit.LeastSquares)"
    grid0 = Yield.KnotGrid(__knot_fit_seed(spline, tenors), tenors, spline; who)
    knots = __KnotRates(spline, Yield.__interpolation(spline, length(grid0.tenors)), grid0.tenors, extrapolation)
    z = __least_squares(model(knots), collect(promote_type(Float64, T, eltype(grid0.rates)), grid0.rates), quotes, primal_quotes, method, who)
    # The public result is validated; knot rates that carry derivatives keep the caller's policy.
    eltype(z) <: ForwardDiff.Dual || return Yield.ZeroRateCurve(z, grid0.tenors, spline; extrapolation = __primal_extrapolation(extrapolation))
    return Yield.__build_public(
        spline, Yield.KnotGrid(z, grid0.tenors, spline; who = "ZeroRateCurve"), knots.interpolation; extrapolation
    )
end

# On the foreign curve an outright receives spot at maturity and pays the discounted strike
# today. Keeping the two legs identifiable gives this quote its own cancellation allowance.
struct __ForwardLegs{C, K} <: FinanceCore.AbstractContract
    receipt::C
    discounted_strike::K
end
FinanceCore.maturity(c::__ForwardLegs) = c.receipt.time
__quote_value(curve, q::Quote{<:Any, <:__ForwardLegs}) = present_value(curve, q.instrument.receipt) - q.instrument.discounted_strike
function __map_contract(f::F, c::__ForwardLegs) where {F}
    receipt, strike = __map_contract(f, c.receipt), f(c.discounted_strike)
    return receipt === c.receipt && strike === c.discounted_strike ? c : __ForwardLegs(receipt, strike)
end

function __foreign_lsq_quote(m, q::Quote{<:Any, <:FX.Forward})
    FX.__check_quote_pair(m, q)
    c = q.instrument
    strike = c.strike * discount(m.domestic, c.time)
    FX.__implied_df(q, m.spot, strike)  # validate the implied discount factor
    return Quote(q.price, __ForwardLegs(Cashflow(m.spot, c.time), strike))
end
__foreign_lsq_quote(m, q::Quote{<:Any, <:FX.BasisSwapLeg}) = FX.__implied_foreign_quote(m, q)

__leg_roundoff(strike, p) = 4 * eps(max(abs(oftype(p, strike)), abs(oftype(p, strike) + p)))
__cancellation_allowance(curve, q::Quote{<:Any, <:__ForwardLegs}, p) = __leg_roundoff(q.instrument.discounted_strike, p)
function __cancellation_allowance(m::FX.Forwards, q::Quote{<:Any, <:FX.Forward}, p)
    c = q.instrument
    return __leg_roundoff(c.strike * discount(m.domestic, c.time), p)
end

function fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.LeastSquares; kwargs...) where {P, S, D, F <: Union{Spline.SplineCurve, Yield.AbstractInterpolatedZeroCurve}}
    qs, method = __included_quotes(quotes, method, "fit(FX.Forwards, Fit.LeastSquares)")
    if all(q -> q isa Quote{<:Any, <:Union{FX.Forward, FX.BasisSwapLeg}}, qs)
        # Conversion validates each outright's implied discount factor and carries the
        # domestic curve derivative through its discounted strike. Other contracts need
        # the full FX model for pricing.
        foreign = __fit_knots(mod0.foreign, map(q -> __foreign_lsq_quote(mod0, q), qs), method, identity; kwargs...)
    else
        model(knots) = __ForeignKnotRates(mod0.pair, mod0.spot, mod0.domestic, knots)
        foreign = __fit_knots(mod0.foreign, qs, method, model; kwargs...)
    end
    return @set mod0.foreign = foreign
end
