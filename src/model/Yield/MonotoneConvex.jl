"""
    MonotoneConvex(rates, tenors; extrapolation=:flat_forward)

A Monotone Convex yield curve model implementing the Hagan-West interpolation method.
`ZeroRateCurve(rates, tenors)` (whose default method is `Spline.MonotoneConvex()`) returns the
same curve, and loss-fitting `Spline.MonotoneConvex()` returns one.

With the default extrapolation, this interpolation method guarantees:
- Continuous forward rates (including at and beyond the last knot, where the
  forward is extrapolated flat at the boundary instantaneous forward `f(t_n)`; unlike
  the DataInterpolations-backed curves, which anchor `:flat_forward` on the last discrete
  forward, MonotoneConvex keeps its own boundary forward)
- Positive forward rates (when input rates imply positive discrete forwards)
- Monotone convex forward curves that match discrete forward rates at knot points

The first interval runs from `t = 0`, with the instantaneous forward at the origin set by the
Hagan-West boundary condition, so the short end is part of the construction.

The `extrapolation` keyword also accepts `:flat_zero`, `:linear`, and
[`FlatForwardAt`](@ref). These change only the tail strictly beyond the last knot;
`:flat_zero` and a supplied forward usually introduce a forward jump there.
`:extension` is unavailable because this model has no DataInterpolations polynomial
piece to continue.

Negative input rates are supported: the Hagan-West positivity collar is applied
symmetrically (bounding node forwards between 0 and twice the adjacent discrete
forwards), though the positivity guarantee is only meaningful when the discrete
forwards themselves are positive.

The implementation follows the Hagan-West method as described in WILMOTT magazine.

# Examples
```julia
prices = [0.98, 0.955, 0.92, 0.88, 0.830]
tenors = [1, 2, 3, 4, 5]
rates = @. -log(prices) / tenors

c = Yield.MonotoneConvex(rates, tenors)
zero(c, 2.5)  # Get the zero rate at t=2.5
discount(c, 2.5)  # Get the discount factor at t=2.5
```

`Yield.MonotoneConvex` is a [`Yield.AbstractInterpolatedZeroCurve`](@ref): it copies and
validates its inputs like every knot curve (a single knot is allowed and gives a flat curve),
its knots are read with [`knot_rates`](@ref)/[`knot_tenors`](@ref), and it changes only through
[`reconstruct`](@ref), which recomputes the node forwards. [`Yield.instantaneous_forward`](@ref) gives
the interpolated instantaneous forward.

# Derivatives with respect to the knot rates

The interpolant switches formula where two adjacent discrete forwards are equal (every flat
stretch of the curve) and where a node forward meets its positivity bound, so it is only
piecewise smooth in its knot rates. With ForwardDiff dual knot rates, derivatives are exact away
from those points; at them, each partial is the limit of a centered bump in its direction. See
"Sensitivities Through Calibration" in the documentation.

# References
- Hagan & West, "Interpolation Methods for Curve Construction", WILMOTT magazine
- Dehlbom, "Interpolation of the yield curve" (http://uu.diva-portal.org/smash/get/diva2:1477828/FULLTEXT01.pdf)
"""
struct MonotoneConvex{T, U, P, E} <: AbstractInterpolatedZeroCurve
    spline::Sp.MonotoneConvex    # stored so every knot curve has the same public properties
    rates::ReadOnlyVector{T}     # continuously-compounded zero rates at the knots
    tenors::ReadOnlyVector{U}    # finite, ≥ 0, strictly increasing
    extrapolation::P             # validated long-end policy
    _f::Vector{T}                # node instantaneous forwards, derived from (rates, tenors)
    _fᵈ::Vector{T}               # discrete forwards, derived from (rates, tenors)
    _tail::E                     # boundary quantities derived from the policy and knots
    # The only inner constructor: takes an owned grid (validated, or an optimizer trial grid) and
    # computes the forwards from it, so they are consistent with the knots by construction.
    function MonotoneConvex(g::KnotGrid{T, U}; extrapolation = :flat_forward) where {T, U}
        extrapolation = __monotone_extrapolation_method(extrapolation)
        f, fᵈ = __monotone_convex_fs(g.rates, g.tenors)
        t, z = last(g.tenors), last(g.rates)
        # `:flat_forward` keeps the native Hagan-West boundary forward f(tₙ); `:linear` uses
        # the zero-rate slope it implies, (f(tₙ) - zₙ)/tₙ.
        forward = () -> last(f)
        slope = () -> iszero(t) ? zero(z) : (last(f) - z) / t
        tail = __build_tail(extrapolation, t, z, forward, slope)
        return new{T, U, typeof(extrapolation), typeof(tail)}(
            Sp.MonotoneConvex(), ReadOnlyVector(g.rates), ReadOnlyVector(g.tenors), extrapolation, f, fᵈ, tail
        )
    end
end

# Public form: copies, promotes and validates through the shared knot-grid path.
MonotoneConvex(rates, tenors; extrapolation = :flat_forward) =
    MonotoneConvex(KnotGrid(rates, tenors, Sp.MonotoneConvex(); who = "Yield.MonotoneConvex"); extrapolation)


function __issector1(g0, g1)
    a = (g0 > 0) && (g1 >= -2 * g0) && (-0.5 * g0 >= g1)
    b = (g0 < 0) && (g1 <= -2 * g0) && (-0.5 * g0 <= g1)
    return a || b
end

# `>=`/`<=` (not `>`/`<`) so the boundary g0 == 0 (left node forward equals the
# discrete forward) is classified here, giving η == 1. Otherwise it falls through
# to sector (iii) with η == 3g1/(g1-g0) == 3 > 1, where the [0, η] sub-region the
# formula assumes within [0, 1] does not exist, so ∫₀¹ g ≠ 0 and the interval
# fails to reprice its right knot (a small zero-rate discontinuity). For g0 != 0
# the bounds are unchanged. (The η == 1 singularity at x == 1 is handled by the
# boundary guard in `g`/`g_rate`.)
function __issector2(g0, g1)
    a = (g0 >= 0) && (g1 < -2 * g0)
    b = (g0 <= 0) && (g1 > -2 * g0)
    return a || b
end


# Deviations within `16 eps` of the discrete forward are treated as zero: a flat interval (both
# negligible) takes the linear fallback, avoiding the 0/0 of the sector formulas.
__mc_tol(fᵈ) = 16 * eps(float(fᵈ))
_mc_negligible(g0, g1, fᵈ) = abs(g0) <= __mc_tol(fᵈ) && abs(g1) <= __mc_tol(fᵈ)

# ── Hagan-West sector formulas ─────────────────────────────────────────────────────────
# On an interval with boundary deviations g0 = g(0) and g1 = g(1) from the discrete forward,
# `__g_*` is the deviation g(x) = f(x) - fᵈ of the instantaneous forward and `__G_*` its
# integral G(x) = ∫₀ˣ g(u) du, at an interior position 0 < x < 1.
__g_i(x, g0, g1) = g0 * (1 - 4 * x + 3 * x^2) + g1 * (-2 * x + 3 * x^2)
__G_i(x, g0, g1) = g0 * (x - 2 * x^2 + x^3) + g1 * (-x^2 + x^3)

function __g_ii(x, g0, g1)
    η = (g1 + 2 * g0) / (g1 - g0)
    return x < η ? g0 : g0 + (g1 - g0) * ((x - η) / (1 - η))^2
end
function __G_ii(x, g0, g1)
    η = (g1 + 2 * g0) / (g1 - g0)
    return x < η ? g0 * x : g0 * x + ((g1 - g0) * (x - η)^3 / (1 - η)^2) / 3
end

function __g_iii(x, g0, g1)
    η = 3 * g1 / (g1 - g0)
    return x > η ? g1 : g1 + (g0 - g1) * ((η - x) / η)^2
end
function __G_iii(x, g0, g1)
    η = 3 * g1 / (g1 - g0)
    return x > η ? (2 * g1 + g0) / 3 * η + g1 * (x - η) : g1 * x - (g0 - g1) * ((η - x)^3 / η^2 - η) / 3
end

function __g_iv(x, g0, g1)
    η = g1 / (g1 + g0)
    α = -g0 * g1 / (g1 + g0)
    return x < η ? α + (g0 - α) * ((η - x) / η)^2 : α + (g1 - α) * ((x - η) / (1 - η))^2
end
function __G_iv(x, g0, g1)
    η = g1 / (g1 + g0)
    α = -g0 * g1 / (g1 + g0)
    return if x < η
        α * x - (g0 - α) * ((η - x)^3 / η^2 - η) / 3
    else
        (2 * α + g0) / 3 * η + α * (x - η) + (g1 - α) / 3 * (x - η)^3 / (1 - η)^2
    end
end

# Linear fallback for a flat interval: g(x) = g0(1-x) + g1·x and its integral.
__g_flat(x, g0, g1) = g0 * (1 - x) + g1 * x
__G_flat(x, g0, g1) = g0 * x + (g1 - g0) * x^2 / 2  # ∫₀ˣ [g0(1-u) + g1·u] du

# The deviation (`__MC_DEVIATION`) and integrated-deviation (`__MC_INTEGRAL`) formulas.
const __MC_DEVIATION = (i = __g_i, ii = __g_ii, iii = __g_iii, iv = __g_iv, flat = __g_flat)
const __MC_INTEGRAL = (i = __G_i, ii = __G_ii, iii = __G_iii, iv = __G_iv, flat = __G_flat)

function __mc_sector(k, x, g0, g1)
    return if sign(g0) == sign(g1)
        k.iv(x, g0, g1)
    elseif __issector1(g0, g1)
        k.i(x, g0, g1)
    elseif __issector2(g0, g1)
        k.ii(x, g0, g1)
    else
        k.iii(x, g0, g1)
    end
end

__mc_kernel(k, x, g0, g1, fᵈ) = _mc_negligible(g0, g1, fᵈ) ? k.flat(x, g0, g1) : __mc_sector(k, x, g0, g1)

# ── Derivatives at kinks ───────────────────────────────────────────────────────────────
# The interpolant is piecewise smooth in the knot rates: it has kinks where a boundary
# deviation is zero, i.e. where two adjacent discrete forwards are equal (a flat interval when
# both are). There, knot-rate derivatives report the limit of a centered bump in the direction
# of each partial. The value is always the primal value.
# - One deviation zero: two sectors meet along the axis (g0 = 0: sectors (ii) and (iv);
#   g1 = 0: (iii) and (iv)). The centered limit is the average of the two sectors' partials,
#   which is linear in the direction, so key-rate sensitivities still add up to the parallel one.
# - Both zero (a flat interval): g and G are odd and positively homogeneous of degree one in
#   (g0, g1), so a bump in direction d moves them by the formula evaluated at d, equally up and
#   down. This is not linear in d: sensitivities to single knots need not add up to the parallel
#   sensitivity.
__mc_kernel(k, x, g0, g1, fᵈ, ktol) = __mc_kernel(k, x, g0, g1, fᵈ)
function __mc_kernel(k, x, g0::D, g1::D, fᵈ, ktol = nothing) where {D <: ForwardDiff.Dual}
    tol = max(__mc_tol(ForwardDiff.value(fᵈ)), something(ktol, zero(ForwardDiff.value(fᵈ))))
    v0, v1 = ForwardDiff.value(g0), ForwardDiff.value(g1)
    n0, n1 = abs(v0) <= tol, abs(v1) <= tol
    (n0 || n1) || return __mc_sector(k, x, g0, g1)
    __ad_depth(D) == 1 || throw(
        __mc_kink_error("second derivatives do not exist where adjacent discrete forwards are equal")
    )
    xv = x isa ForwardDiff.Dual ? ForwardDiff.value(x) : x
    value = __mc_kernel(k, xv, v0, v1, ForwardDiff.value(fᵈ))
    partials = if n0 && n1
        p0, p1 = ForwardDiff.partials(g0), ForwardDiff.partials(g1)
        ForwardDiff.Partials(ntuple(j -> __mc_tangent(k, xv, p0[j], p1[j]), Val(ForwardDiff.npartials(D))))
    else
        a, b = n0 ? (k.ii, k.iv) : (k.iii, k.iv)
        (ForwardDiff.partials(a(x, g0, g1)) + ForwardDiff.partials(b(x, g0, g1))) / 2
    end
    return D(value, partials)
end

# The flat-interval formula evaluated at a tangent (d0, d1): zero for a parallel direction.
__mc_tangent(k, x, d0, d1) = iszero(d0) && iszero(d1) ? zero(d0 * x) : __mc_sector(k, x, d0, d1)

# Rounding in the discrete forwards (differences of t·z over short intervals) can leave the
# deviations of a flat interval slightly nonzero, and the sector formulas then differentiate
# rounding noise. With dual knot rates, deviations of interval `i` within this bound count as
# zero. It covers the intervals that the node forwards at either end of `i` average over (the
# first, from 0 to a first knot at t = 0, has no width and no rounding).
function __mc_kink_tol(rates, times, i)
    s = zero(float(__primal(first(rates))))
    for j in max(1, i - 1):min(lastindex(times), i + 1)
        tp, zp = j == 1 ? (zero(first(times)), zero(s)) : (times[j - 1], __primal(rates[j - 1]))
        times[j] > tp || continue
        s = max(s, (abs(times[j] * __primal(rates[j])) + abs(tp * zp)) / (times[j] - tp))
    end
    return 16 * eps(s)
end
__mc_kink_tol(mc, i) = eltype(mc.rates) <: ForwardDiff.Dual ? __mc_kink_tol(mc.rates, mc.tenors, i) : nothing

__mc_kink_error(why) = ArgumentError(
    "Yield.MonotoneConvex cannot differentiate with respect to its knot rates here: $why. Use a " *
        "smooth interpolation (e.g. Spline.Linear() or Spline.Cubic()), key-rate shifts added to the " *
        "curve (ActuaryUtilities `KeyRates`), or finite bumps."
)

"""
    g(x, f⁻, f, fᵈ)

Compute the deviation of the instantaneous forward rate from the discrete forward
rate at normalized position x ∈ [0, 1] within an interval. Following Hagan-West,
g(x) = f(x) - fᵈ with boundary conditions g₀ = f⁻ - fᵈ and g₁ = f - fᵈ that determine
the sector-specific polynomial used for interpolation.
"""
function g(x, f⁻, f, fᵈ, ktol = nothing)
    g0 = f⁻ - fᵈ
    g1 = f - fᵈ
    # Interval-boundary deviations are exact: g(0) = g0 and g(1) = g1. Returning
    # them directly avoids the 0/0 the sector formulas produce at a degenerate η:
    # sector (iii) has `/η` and is singular when g1 == 0 (η == 0, hit at every
    # interior knot, x == 0); sector (ii) has `/(1-η)` and is singular when g0 == 0
    # (η == 1, hit at the last knot, x == 1). An `iszero(η)` guard is *not* safe
    # under ForwardDiff — η carries nonzero partials despite a zero value, so the
    # singular branch still runs and yields NaN gradients. `x` is the normalized
    # interval position — a plain Float at the boundaries, never the AD variable.
    iszero(x) && return g0
    isone(x) && return g1
    return __mc_kernel(__MC_DEVIATION, x, g0, g1, fᵈ, ktol)
end

"""
    g_rate(x, f⁻, f, fᵈ)

Compute the integrated deviation G(x) = ∫₀ˣ g(u) du, which captures how the
instantaneous forward curve deviates from the discrete forward across an interval.
This quantity feeds into the zero-rate relation r(t) = fᵈ + (Δt / t) ⋅ G(x) used by
the Hagan-West construction.
"""
function g_rate(x, f⁻, f, fᵈ, ktol = nothing)
    g0 = f⁻ - fᵈ
    g1 = f - fᵈ
    # G(0) = ∫₀⁰ g = 0 and G(1) = ∫₀¹ g = 0 are both identities (the interpolated
    # forward integrates to the discrete forward over each interval). Returning
    # zero directly avoids the 0/0 the sector formulas produce at a degenerate η:
    # sector (iii) `/η²` is singular at g1 == 0 (η == 0, interior knots, x == 0),
    # sector (ii) `/(1-η)²` at g0 == 0 (η == 1, last knot, x == 1). See `g` for why
    # an `iszero(η)` guard is unsafe under ForwardDiff and `x` is the discriminator.
    (iszero(x) || isone(x)) && return zero(g0 * x)
    return __mc_kernel(__MC_INTEGRAL, x, g0, g1, fᵈ, ktol)
end

"""
    instantaneous_forward(mc::MonotoneConvex, t)

The instantaneous (continuously-compounded) forward rate of the Hagan-West
interpolant at time `t`. Beyond the last knot it follows `mc.extrapolation`.
The default `:flat_forward` holds the boundary forward constant. `:flat_zero`
holds the forward at the last zero rate, `:linear` derives it from the linearly
extended zero rate, and `FlatForwardAt(f)` holds it at `f`.
At the last knot itself this method returns the interior (left) forward.

Note this is distinct from `forward(curve, from, to)`, which is the *discrete*
forward `Rate` between two times and is defined for every yield model.
"""
function instantaneous_forward(mc::MonotoneConvex, t)
    __check_time(t, "instantaneous_forward")
    f, fᵈ, times = mc._f, mc._fᵈ, mc.tenors
    lt = last(times)

    # At the boundary report the interior (left) forward. Policies that impose a
    # different forward begin strictly after the knot, consistently with `zero`.
    if t > lt
        return __tail_forward(mc._tail, t)
    elseif t == lt
        return f[end]
    end

    i_time = __i_time(t, times)

    if i_time == 1
        # First interval: from 0 to times[1]
        x = t / times[1]
        return fᵈ[1] + g(x, f[1], f[2], fᵈ[1], __mc_kink_tol(mc, 1))
    else
        # Interval from times[i_time-1] to times[i_time]
        t_prev = times[i_time - 1]
        t_curr = times[i_time]
        x = (t - t_prev) / (t_curr - t_prev)
        return fᵈ[i_time] + g(x, f[i_time], f[i_time + 1], fᵈ[i_time], __mc_kink_tol(mc, i_time))
    end
end

"""
    returns the index associated with the time t, an initial rate vector, and a time vector
"""
function __i_time(t, times)
    # first interval whose right endpoint is ≥ t, i.e. the index after the last
    # knot that is ≤ t (times sorted ascending); O(log n) vs the prior findfirst
    i_time = searchsortedlast(times, t) + 1
    return min(i_time, lastindex(times))
end

# Hagan-West positivity collar, generalized for negative discrete forwards:
# bound the node forward between 0 and twice the adjacent discrete forward(s).
# For positive forwards this is the paper's clamp(f, 0, 2m); for negative ones
# the interval flips to [2m, 0] rather than producing an inverted (lo > hi) clamp.
__collar(f, m) = clamp(f, min(zero(m), 2 * m), max(zero(m), 2 * m))
__collar(f, a, b) = __collar(f, min(a, b))

# With dual knot rates the collar is a kink where the raw forward meets the bound exactly: there
# the partials are the average of the raw forward's and the bound's (the limit of a centered
# bump). A tie that coincides with another kink (a zero discrete forward, tied neighbors, or an
# adjacent interval whose deviation is zero) has no such limit and throws.
__collar(f::D, m::D) where {D <: ForwardDiff.Dual} = __dual_collar(f, (m,))
__collar(f::D, a::D, b::D) where {D <: ForwardDiff.Dual} = __dual_collar(f, (a, b))
function __dual_collar(f::D, ms) where {D <: ForwardDiff.Dual}
    fv, mvs = ForwardDiff.value(f), map(ForwardDiff.value, ms)
    mv = minimum(mvs)
    lo, hi = min(zero(mv), 2 * mv), max(zero(mv), 2 * mv)
    lo < fv < hi && return f   # the collar does not bind
    tied = [m for m in ms if ForwardDiff.value(m) == mv]
    (iszero(mv) || !allequal(map(ForwardDiff.partials, tied))) &&
        throw(__mc_kink_error("the positivity collar binds at a node whose discrete forwards are zero or tied"))
    m = first(tied)
    bound = ((fv <= lo) == (mv < 0)) ? 2 * m : zero(m)
    (fv < lo || fv > hi) && return bound
    __ad_depth(D) == 1 || throw(__mc_kink_error("second derivatives do not exist where the positivity collar binds exactly"))
    any(m -> abs(fv - ForwardDiff.value(m)) <= __mc_tol(ForwardDiff.value(m)), ms) &&
        throw(__mc_kink_error("the positivity collar binds exactly at a node whose adjacent interval is flat"))
    return D(fv, (ForwardDiff.partials(f) + ForwardDiff.partials(bound)) / 2)
end

"""
    returns a pair of vectors (f and fᵈ) used in Monotone Convex Yield Curve fitting
"""
function __monotone_convex_fs(rates, times)
    f, fᵈ = __mc_raw_forwards(rates, times)
    length(rates) == 1 && return f, fᵈ
    # Collar each node against the discrete forwards of its adjacent intervals.
    # With f[j] = f(t_{j-1}), node t_{j-1} adjoins intervals j-1 and j, so the
    # interior bound is min(fᵈ[j-1], fᵈ[j]); the endpoints have one neighbor each.
    f[1] = __collar(f[1], fᵈ[1])
    f[end] = __collar(f[end], fᵈ[end])
    for j in 2:length(times)
        f[j] = __collar(f[j], fᵈ[j - 1], fᵈ[j])
    end
    return f, fᵈ
end

# The discrete forwards `fᵈ` and the node forwards `f` before the positivity collar.
function __mc_raw_forwards(rates, times)
    # step 1
    fᵈ = copy(rates)
    for i in 2:length(times)
        fᵈ[i] = (times[i] * rates[i] - times[i - 1] * rates[i - 1]) / (times[i] - times[i - 1])
    end
    # step 2
    # Convention: f[j] is the instantaneous forward at node t_{j-1} (with t_0 = 0),
    # so f has length n+1, f[1] = f(0) and f[end] = f(t_n).
    f = similar(rates, length(rates) + 1)
    if length(rates) == 1
        # single interval: flat forward (the boundary adjustments below would
        # otherwise read f[2] before it is initialized)
        f[1] = fᵈ[1]
        f[2] = fᵈ[1]
        return f, fᵈ
    end
    # fill in middle elements first, then do 1st and last
    for i in 1:(length(rates) - 1)
        t_prior = if i == 1
            0
        else
            times[i - 1]
        end

        weight1 = (times[i] - t_prior) / (times[i + 1] - t_prior)
        weight2 = (times[i + 1] - times[i]) / (times[i + 1] - t_prior)
        f[i + 1] = weight1 * fᵈ[i + 1] + weight2 * fᵈ[i]
    end
    # step 3: boundary node forwards (the positivity collar follows in `__monotone_convex_fs`).
    f[1] = fᵈ[1] - 0.5 * (f[2] - fᵈ[1])
    f[end] = fᵈ[end] - 0.5 * (f[end - 1] - fᵈ[end])
    return f, fᵈ
end

# Kinks in the knot rates `z`: the boundary deviations of each interval, and each node's raw
# forward against its collar bounds. A single knot gives a flat curve for every `z`: no kink.
function __kink_quantities(::Sp.MonotoneConvex, z, tenors)
    length(z) == 1 && return eltype(z)[]
    f, fᵈ = __monotone_convex_fs(z, tenors)
    raw, _ = __mc_raw_forwards(z, tenors)
    q = eltype(f)[]
    for i in eachindex(fᵈ)
        push!(q, f[i] - fᵈ[i], f[i + 1] - fᵈ[i])
    end
    for j in eachindex(raw)
        m = j == 1 ? fᵈ[1] : j == lastindex(raw) ? fᵈ[end] : min(fᵈ[j - 1], fᵈ[j])
        push!(q, raw[j] - min(zero(m), 2 * m), raw[j] - max(zero(m), 2 * m))
    end
    return q
end


function Base.zero(mc::MonotoneConvex, t)
    __check_time(t, "zero")
    f, fᵈ, rates, times = mc._f, mc._fᵈ, mc.rates, mc.tenors
    lt = last(times)

    # Handle t=0 case (limit is instantaneous forward at t=0)
    if t == 0
        return Continuous(f[1])
    end

    if t > lt
        return Continuous(mc._tail(t))
    end

    i_time = __i_time(t, times)

    # Calculate normalized position x in interval and interval bounds
    if i_time == 1
        # First interval: from 0 to times[1]
        x = t / times[1]
        G = g_rate(x, f[1], f[2], fᵈ[1], __mc_kink_tol(mc, 1))
        # r(t) = (1/t) * [t * fᵈ[1] + times[1] * G(x)]
        return Continuous(fᵈ[1] + times[1] * G / t)
    else
        # Interval from times[i_time-1] to times[i_time]
        t_prev = times[i_time - 1]
        t_curr = times[i_time]
        x = (t - t_prev) / (t_curr - t_prev)
        G = g_rate(x, f[i_time], f[i_time + 1], fᵈ[i_time], __mc_kink_tol(mc, i_time))
        # r(t) = (1/t) * [t_prev * r_prev + (t - t_prev) * fᵈ[i] + (t_curr - t_prev) * G(x)]
        return Continuous((t_prev * rates[i_time - 1] + (t - t_prev) * fᵈ[i_time] + (t_curr - t_prev) * G) / t)
    end
end
function FinanceCore.discount(mc::MonotoneConvex, t)
    __check_time(t, "discount")
    return isinf(t) ? __discount_at_infinity(mc._tail) : _discount_from_zero(mc, t)
end
