## Blends of curves

"""
    Yield.DiscountFactors()

Blend the curves' discount factors. See [`Blend`](@ref FinanceModels.Yield.Blend).
"""
struct DiscountFactors end

"""
    Yield.ZeroRates()

Blend the curves' continuous zero rates. See [`Blend`](@ref FinanceModels.Yield.Blend).
"""
struct ZeroRates end

"""
    Yield.ForwardRates(period)
    Yield.ForwardRates()

Blend the curves' instantaneous forward rates. A weight that varies with tenor is read once per
`period`, at the period's start. `ForwardRates()` is for constant weights. See
[`Blend`](@ref FinanceModels.Yield.Blend).
"""
struct ForwardRates{P}
    period::P
    function ForwardRates(period::Real)
        p = __primal(period)
        (isfinite(p) && p > 0) || throw(ArgumentError("a `ForwardRates` period must be finite and positive, got $period"))
        return new{typeof(period)}(period)
    end
    ForwardRates() = new{Nothing}(nothing)
end

const __BlendSpace = Union{DiscountFactors, ZeroRates, ForwardRates}
const __NumberWeights = Union{Tuple{Vararg{Real}}, ReadOnlyVector{<:Real}}

"""
    Yield.Blend(curves, weights, space)
    Yield.Blend(a, b, w, space)

A curve that combines `curves` with `weights` in `space`:
[`DiscountFactors()`](@ref FinanceModels.Yield.DiscountFactors),
[`ZeroRates()`](@ref FinanceModels.Yield.ZeroRates) or
[`ForwardRates(period)`](@ref FinanceModels.Yield.ForwardRates). The two-curve form weights `a`
by `w` and `b` by `1 - w`.

`curves` and `weights` are tuples or vectors of the same length; vectors are copied. A weight is a
number, or a function of tenor `t -> w`.

| Space | The blend | Weights |
|:------|:----------|:--------|
| `DiscountFactors()` | ``D(t) = \\sum_i \\tilde w_i D_i(t)``, with ``\\tilde w_i = w_i / \\sum_j w_j`` | in [0, 1], summing to 1 |
| `ZeroRates()` | ``z(t) = \\sum_i w_i z_i(t)``, so ``-\\log D(t) = \\sum_i w_i \\cdot (-\\log D_i(t))`` | summing to 1 |
| `ForwardRates(p)` | ``f(t) = \\sum_i w_i(kp) f_i(t)`` for ``kp \\le t < (k+1)p`` | summing to 1 |

Number weights are checked when the blend is built, and weights that vary with tenor each time they
are read: they must be finite and sum to 1, to rounding, and discount-factor weights must lie in
[0, 1]. Otherwise an `ArgumentError` is thrown.

- **Discount factors.** The weights are shares of their sum. With constant weights, the present
  value of cashflows valued from time 0 is linear in the curves: `pv(blend, cfs) ≈ Σ wᵢ·pv(cᵢ, cfs)`.
  The long end follows the curve with the lowest rate. Mathematically, a derivative with respect to
  one of N weights is `(Dᵢ - D)/Σw`; in the two-curve form, `∂D/∂w = D_a - D_b`. Numerical limitations
  of automatic differentiation are described below.
  - A blend rebased to time `s` (`ForwardStarting`) mixes its curves with the weights
    `w̃ᵢ·Dᵢ(s)/D(s)`, not `w̃ᵢ`.
  - An expected discount factor gives the expected present value only of cashflows that don't vary
    by scenario.
- **Zero rates.** The weights are used as given, so negative weights extrapolate. With constant
  weights, `Blend(a, b, w, ZeroRates())` is the same curve as `w * a + (1 - w) * b`.
- **Forward rates.** The blend is built from the curves' interval factors, so it starts at
  `discount(blend, 0) == 1` whatever a curve's own value at 0. With constant weights it equals the
  zero-rate blend of curves with `discount(c, 0) == 1`. A weight that varies with tenor costs about
  `t / period` interval factors of each curve per evaluation at `t`.

With weights that vary with tenor, the limits at `t = Inf` are `NaN`, since the weights there are not
known.

!!! warning "Automatic differentiation at extreme scales"
    In discount-factor blends, extremely small weights or widely separated component log-discounts
    can produce inaccurate derivatives, including spurious zeros, `Inf` or `NaN`, even when the
    price and mathematical derivative are finite. This affects differentiation with respect to
    weights, time and curve parameters. A derivative that works for a direct discount calculation
    may fail through zero rates, interval discounts or curve composition. Verify sensitivities in
    these regimes against analytic or suitably higher-precision references.

# Examples

```julia
market = ZeroRateCurve([0.03, 0.035, 0.04], [1.0, 10.0, 20.0])
ultimate = Yield.Constant(Continuous(0.045))

# grade the market curve's annual forwards into the ultimate rate between 20 and 60 years
grade(t) = clamp((60 - t) / 40, 0, 1)
Yield.Blend(market, ultimate, grade, Yield.ForwardRates(1.0))

# the expected discount factor over three scenarios
Yield.Blend((s1, s2, s3), (0.5, 0.3, 0.2), Yield.DiscountFactors())
```
"""
struct Blend{C, W, S <: __BlendSpace} <: AbstractYieldModel
    curves::C
    weights::W
    space::S
    function Blend(
            curves::Union{Tuple{Vararg{AbstractYieldModel}}, AbstractVector{<:AbstractYieldModel}}, weights,
            space::__BlendSpace
        )
        cs, ws = __owned(curves), __owned(weights)
        length(cs) == length(ws) || throw(DimensionMismatch("a blend of $(length(cs)) curves has $(length(ws)) weights"))
        if ws isa __NumberWeights
            __check_weights(map(__primal, ws), space)
        elseif space isa ForwardRates{Nothing}
            throw(ArgumentError("weights that vary with tenor need a period: `ForwardRates(period)`"))
        end
        return new{typeof(cs), typeof(ws), typeof(space)}(cs, ws, space)
    end
end
Blend(a::AbstractYieldModel, b::AbstractYieldModel, w, space::__BlendSpace) = Blend((a, b), (w, __complement(w)), space)

# A blend's curves and weights: a tuple as it is, a vector copied into a `ReadOnlyVector`, so that
# neither the caller's array nor the stored field can change what was checked.
__owned(x::Tuple) = x
__owned(x::ReadOnlyVector) = x
__owned(x::AbstractVector) = ReadOnlyVector(collect(x))

# The second weight of the two-curve form.
__complement(w::Real) = 1 - w
__complement(w) = __Complement(w)
struct __Complement{F}
    w::F
end
(c::__Complement)(t) = 1 - c.w(t)

# Weights, as primal values, must be finite and sum to 1 to rounding. The bound on the sum's error
# grows with the weights' size, but is never wider than √eps, so offsetting weights such as
# (1e16, -1e16), whose sum is 0, are rejected. A sum with a non-finite weight fails the bound too.
# Discount-factor weights must also lie in [0, 1].
function __check_weights(ws, space)
    F = float(mapreduce(typeof, promote_type, ws))
    s = sum(ws)
    tol = min(sqrt(eps(F)), 4 * length(ws) * eps(F) * max(one(F), sum(abs, ws)))
    abs(s - 1) <= tol || throw(ArgumentError("blend weights must be finite and sum to 1, got $(Tuple(ws)), summing to $s"))
    space isa DiscountFactors && !all(w -> 0 <= w <= 1, ws) &&
        throw(ArgumentError("discount-factor weights must lie in [0, 1], got $(Tuple(ws))"))
    return nothing
end

# The weights at tenor `t`: number weights as they are (checked at construction), and weights that
# vary with tenor read and checked at `t`.
__weight_at(w::Real, t) = w
__weight_at(w, t) = w(t)
__weights(b::Blend{<:Any, <:__NumberWeights}, t) = b.weights
function __weights(b::Blend, t)
    ws = map(w -> __weight_at(w, t), b.weights)
    __check_weights(map(__primal, ws), b.space)
    return ws
end

# A zero-rate blend with number weights combines like `CompositeYield` and `ScaledYield`: its
# log-discount, zero rate, forward, intervals and tail are each the weighted sum of its curves' (see
# `__CombinedYield` in Composition.jl, whose methods are more specific than `Blend`'s below).
const __ZeroBlend = Blend{<:Any, <:__NumberWeights, ZeroRates}

__curve_log_discounts(b, t) = map(@inline(c -> __log_discount(c, t)), b.curves)
# Σ f(xᵢ, yᵢ, …). Over tuples it is unrolled and inlined, left to right, so that each curve's own
# method inlines at the call as it does for `CompositeYield`; over vectors it is `sum`.
@inline __mapsum(f::F, xs::Vararg{Tuple, N}) where {F, N} = __mapsum_from(f, f(map(first, xs)...), map(Base.tail, xs)...)
@inline __mapsum_from(f::F, acc, ::Vararg{Tuple{}, N}) where {F, N} = acc
@inline __mapsum_from(f::F, acc, xs::Vararg{Tuple, N}) where {F, N} =
    __mapsum_from(f, acc + f(map(first, xs)...), map(Base.tail, xs)...)
__mapsum(f::F, xs::Vararg{Any, N}) where {F, N} = sum(map(f, xs...))

# ── Discount factors ─────────────────────────────────────────────────────────────────────

# The discount factor D and the log-discount L = -log D of the mixture Σ w̃ᵢ·exp(-Lᵢ), with shares
# w̃ᵢ = wᵢ/s and s = Σw, all in the promoted type of the weights and the log-discounts. Products with
# a weight are strong zeros (`__strong_zero_mul`), and weighted exponentials w·eˣ are `__wexp`, which
# evaluates eˣ·∂w and (w·eˣ)·∂x separately to reduce intermediate overflow and underflow. The curves
# with a positive share P are taken relative to the one with the smallest Lᵢ, m. The curves with a
# zero share Z add 0 in value but carry their weights' partials.
#
# D = S·e^{-m} + Σ_Z w̃ᵢ·e^{-Lᵢ}, with S = Σ_P wᵢ·e^{m-Lᵢ}/s, which lies between the reference's share
# and 1. Scaling can keep the value and some derivatives representable when a component discount
# factor overflows. It does not preserve every finite derivative at extreme scales: relative terms
# can underflow, and differentiating the reference can introduce cancellation. See the known AD
# limitations in test/Blend.jl. When every Lᵢ = 0, the sum in S is s itself, so D = 1 exactly.
#
# L = L_P - log1p(Σ_Z w̃ᵢ·e^{L_P-Lᵢ}). L_P = m - log1p(x), with x = Σ_P w̃ᵢ·expm1(m - Lᵢ) - Σ_Z w̃ᵢ
# (the last sum is 0 in value but carries the partials of Σ_P w̃ᵢ - 1, which log1p would otherwise
# drop), or, where x < -1/2 (the reference has a small share), a log-sum-exp of log w̃ᵢ + m - Lᵢ,
# with the logs of the shares as `__strong_zero_log`, so that a tiny constant share has zero partials.
# At t = 0, L = 0 exactly.
function __mixture(ws, Ls)
    T = promote_type(mapreduce(typeof, promote_type, ws), mapreduce(typeof, promote_type, Ls))
    w, L = map(T, ws), map(T, Ls)
    s = sum(w)
    w̃ = map(x -> x / s, w)
    pos = map(x -> __primal(x) > 0, w̃)
    j = argmin(i -> pos[i] ? __primal(L[i]) : oftype(__primal(L[i]), Inf), eachindex(L))
    m = L[j]
    S = sum(i -> pos[i] ? __wexp(w[i], m - L[i]) : zero(T), eachindex(L)) / s
    D = __wexp(S, -m) + sum(i -> pos[i] ? zero(T) : __wexp(w̃[i], -L[i]), eachindex(L))
    x = sum(i -> pos[i] ? __strong_zero_mul(w̃[i], expm1(m - L[i])) : -w̃[i], eachindex(L))
    LP = if __primal(x) >= -1 / 2
        m - log1p(x)
    else
        lt(i) = __strong_zero_log(w̃[i]) + (m - L[i])
        k = argmax(i -> pos[i] ? __primal(lt(i)) : oftype(__primal(L[i]), -Inf), eachindex(L))
        M = lt(k)
        m - (M + log(sum(i -> pos[i] ? exp(lt(i) - M) : zero(T), eachindex(L))))
    end
    Lb = LP - log1p(sum(i -> pos[i] ? zero(T) : __wexp(w̃[i], LP - L[i]), eachindex(L)))
    return (; D, L = Lb, shares = w̃, Ls = L)
end
__mixture(b::Blend, t) = __mixture(__weights(b, t), __curve_log_discounts(b, t))

# A discount factor that is a normal float, where a ratio of two keeps its precision.
__is_normal(x) = (p = __primal(x); floatmin(p) <= p < Inf)

__blend_discount(::DiscountFactors, b, t) = __mixture(b, t).D
__blend_log_discount(::DiscountFactors, b, t) = __mixture(b, t).L
__blend_log_interval(::DiscountFactors, b, from, to) = __mixture(b, to).L - __mixture(b, from).L
# The ratio of the discount factors where both are normal floats: D(0) is exactly 1, with zero
# partials, so an interval from 0 is `discount(b, t)` itself. Otherwise the difference of the
# log-discounts, which stays finite where a discount factor underflows or overflows.
# Reconstructing the interval discount from L can lose a finite derivative when a partial of L
# overflows, even if the direct discount calculation preserves that derivative.
function __blend_interval_discount(::DiscountFactors, b, from, to)
    (isinf(from) || isinf(to)) && return exp(-__log_interval(b, from, to))
    x, y = __mixture(b, from), __mixture(b, to)
    return __is_normal(x.D) && __is_normal(y.D) ? y.D / x.D : exp(x.L - y.L)
end
# The forward Σ pᵢ·fᵢ, where pᵢ = w̃ᵢ·Dᵢ/D = w̃ᵢ·e^{L-Lᵢ} (`__wexp`) is a curve's share of the
# discount factor: each curve's own forward, so it is exact at t = 0. At t = Inf it is the forward of
# the curves that decide the tail. With weights that vary with tenor, the derivative of L.
function __blend_force(::DiscountFactors, b, t)
    b.weights isa __NumberWeights || return __log_discount_derivative(b, t)
    __at_infinity(t) && return __mixture_force_at_infinity(b)
    x = __mixture(b, t)
    return __mapsum(@inline((w, Li, c) -> __strong_zero_mul(__wexp(w, x.L - Li), force_of_interest(c, t))), x.shares, x.Ls, b.curves)
end
# The limit of the forward: that of the curve that decides the tail, or of the curves tied with it,
# weighted by their limiting shares w̃ᵢ·e^{-a0ᵢ} where their forwards differ; NaN where the order of
# the tails is unknown.
function __mixture_force_at_infinity(b)
    tails = map(__log_tail, b.curves)
    w̃ = __tail_shares(b.weights, tails)
    fs = map(c -> force_of_interest(c, Inf), b.curves)
    T = promote_type(eltype(w̃), mapreduce(typeof, promote_type, fs))
    d = __dominant_tails(w̃, tails)
    d === nothing && return T(NaN)
    j, tied = d
    all(i -> !tied(i) || fs[i] == fs[j], eachindex(fs)) && return T(fs[j])
    r = __tied_reference(w̃, tails, tied)
    q = map(i -> tied(i) ? __wexp(T(w̃[i]), T(tails[r].a0) - T(tails[i].a0)) : zero(T), eachindex(fs))
    return sum(i -> __strong_zero_mul(q[i], T(fs[i])), eachindex(fs)) / sum(q)
end
function __blend_tail(::DiscountFactors, b)
    tails = map(__log_tail, b.curves)
    return __mixture_tail(__tail_shares(b.weights, tails), tails)
end
# The discount factor at t = Inf. Where the mixture's limit is finite (the deciding curves' L tends to
# a constant), it is the sum Σ w̃ᵢ·e^{-a0ᵢ} over the tied curves, as `__mixture` sums at a finite time,
# which can preserve weight derivatives where e^{-a0ᵢ} underflows. Otherwise it is
# the tail's limit: 0, Inf, or NaN where unknown.
__blend_discount_at_infinity(space, b) = __discount_at_infinity(b)
function __blend_discount_at_infinity(::DiscountFactors, b)
    b.weights isa __NumberWeights || return __discount_at_infinity(b)
    tails = map(__log_tail, b.curves)
    w̃ = __tail_shares(b.weights, tails)
    tail = __mixture_tail(w̃, tails)
    __tail_limit_sign(tail) == 0 || return __discount_at_infinity(tail)
    j, tied = __dominant_tails(w̃, tails)
    T = typeof(tail.a0)
    a0r = T(tails[__tied_reference(w̃, tails, tied)].a0)
    pos(i) = __primal(w̃[i]) > 0
    S = sum(i -> tied(i) && pos(i) ? __wexp(T(w̃[i]), a0r - T(tails[i].a0)) : zero(T), eachindex(tails))
    return __wexp(S, -a0r) + sum(i -> tied(i) && !pos(i) ? __wexp(T(w̃[i]), -T(tails[i].a0)) : zero(T), eachindex(tails))
end
# Number weights as shares, in the type of the tails' coefficients.
function __tail_shares(ws, tails)
    T = promote_type(mapreduce(typeof, promote_type, ws), __tail_type(tails))
    w = map(T, ws)
    s = sum(w)
    return map(x -> x / s, w)
end

# ── Zero rates, with weights that vary with tenor ────────────────────────────────────────

__blend_log_discount(::ZeroRates, b, t) = __mapsum(__strong_zero_mul, __weights(b, t), __curve_log_discounts(b, t))
__blend_discount(s::ZeroRates, b, t) = exp(-__blend_log_discount(s, b, t))
__blend_log_interval(::ZeroRates, b, from, to) = __log_discount(b, to) - __log_discount(b, from)
__blend_interval_discount(::ZeroRates, b, from, to) = exp(-__log_interval(b, from, to))
__blend_force(::ZeroRates, b, t) = __log_discount_derivative(b, t)
__blend_zero(::ZeroRates, b, t) =
    Continuous(__mapsum(@inline((w, c) -> __strong_zero_mul(w, __continuous(Base.zero(c, t)))), __weights(b, t), b.curves))

# ── Forward rates ────────────────────────────────────────────────────────────────────────

# The period K of `t` on the grid k·p: K·p ≤ t < (K + 1)·p for the grid points as computed.
function __period_index(t, p)
    x, q = __primal(t), __primal(p)
    K = floor(Int, x / q)
    while (K + 1) * q <= x
        K += 1
    end
    while K * q > x
        K -= 1
    end
    return K
end

# The log-discount of a forward blend over [from, to], for finite from ≤ to: each curve's own
# interval log-discount, weighted. With number weights that is one interval per curve. With weights
# that vary with tenor, each period k's piece [a, e] of the interval takes the weights at the period's
# start, k·p. The last piece is kept even when it is empty, so that a time derivative at a grid point
# is the right-hand one. A log-native curve's interval is L(e) - L(a) (`__log_native_interval`), so
# its L is evaluated once per grid point instead, with the same result.
function __forward_interval(b, from, to)
    b.weights isa __NumberWeights &&
        return __mapsum(@inline((w, c) -> __strong_zero_mul(w, __log_interval(c, from, to))), b.weights, b.curves)
    p = b.space.period
    T = promote_type(typeof(from), typeof(to), typeof(p))
    K0, K1 = __period_index(from, p), __period_index(to, p)
    a = T(from)
    La = __native_log_discounts(b.curves, a)
    e = K0 == K1 ? T(to) : T((K0 + 1) * p)
    Le = __native_log_discounts(b.curves, e)
    total = __forward_piece(b, K0 * p, a, e, La, Le)
    for k in (K0 + 1):K1
        a, La = e, Le
        e = k == K1 ? T(to) : T((k + 1) * p)
        Le = __native_log_discounts(b.curves, e)
        total += __forward_piece(b, k * p, a, e, La, Le)
    end
    return total
end
# The piece [a, e] of the period starting at `start`, with the log-native curves' L at a and e.
function __forward_piece(b, start, a, e, La, Le)
    return __mapsum(__weights(b, start), b.curves, La, Le) do w, c, la, le
        __strong_zero_mul(w, la === nothing ? __log_interval(c, a, e) : le - la)
    end
end
__native_log_discounts(curves, t) = map(c -> __native_log_discount(c, t), curves)
# L(t) of a log-native curve as its interval from `t` uses it (0 at an exact 0), or `nothing`.
function __native_log_discount(c, t)
    __log_native(c) || return nothing
    L = __log_discount(c, t)
    return iszero(t) ? zero(L) : L
end

function __blend_log_interval(::ForwardRates, b, from, to)
    return __primal(from) <= __primal(to) ? __forward_interval(b, from, to) : -__forward_interval(b, to, from)
end
# From 0 by the ordered interval, so a negative time reaches the curves, which may reject it.
__blend_log_discount(s::ForwardRates, b, t) = __blend_log_interval(s, b, zero(t), t)
__blend_discount(s::ForwardRates, b, t) = exp(-__blend_log_discount(s, b, t))
__blend_interval_discount(::ForwardRates, b, from, to) = exp(-__log_interval(b, from, to))
function __blend_force(s::ForwardRates, b, t)
    f(ws) = __mapsum(@inline((w, c) -> __strong_zero_mul(w, force_of_interest(c, t))), ws, b.curves)
    b.weights isa __NumberWeights && return f(b.weights)
    # the weights at infinity are not known
    __at_infinity(t) && return oftype(__mapsum(c -> force_of_interest(c, t), b.curves), NaN)
    return f(__weights(b, __period_index(t, s.period) * s.period))
end
# Each curve's tail from its interval factors, L(t) - L(0).
__blend_tail(::ForwardRates, b) =
    __mapsum((w, c) -> __strong_zero_mul(w, __rebase_tail(__log_tail(c), __log_discount(c, 0))), b.weights, b.curves)

# ── The curve interface ──────────────────────────────────────────────────────────────────

# Every space except the zero rates with number weights (a `__CombinedYield`) dispatches here.
function __log_discount(b::Blend, t)
    __at_infinity(t) && return __log_discount_at_infinity(b)
    return __blend_log_discount(b.space, b, t)
end
function FinanceCore.discount(b::Blend, t)
    __at_infinity(t) && return __blend_discount_at_infinity(b.space, b)
    return __blend_discount(b.space, b, t)
end
function __log_interval(b::Blend, from, to)
    (isinf(from) || isinf(to)) && return __log_discount(b, to) - __log_discount(b, from)
    return __blend_log_interval(b.space, b, from, to)
end
function FinanceCore.discount(b::Blend, from, to)
    d = __blend_interval_discount(b.space, b, from, to)
    return from == to ? one(d) : d
end
force_of_interest(b::Blend, t) = __blend_force(b.space, b, t)
Base.zero(b::Blend, t) = b.space isa ZeroRates ? __blend_zero(b.space, b, t) : __zero_from_log_discount(b, t)
# With weights that vary with tenor the tail is unknown.
__log_tail(b::Blend) = b.weights isa __NumberWeights ? __blend_tail(b.space, b) : __nan_tail(map(__log_tail, b.curves))
