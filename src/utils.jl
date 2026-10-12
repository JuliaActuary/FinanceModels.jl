# all we need is cumulative normal, so avoid Distributions.jl dependency
# https://www.johndcook.com/blog/cpp_phi/

function __ϕ(x)
    return 0.5 * SpecialFunctions.erfc(-x * √(0.5))
end

__N(x) = __ϕ(x)

# the standard normal density
__normal_density(x) = exp(-x^2 / 2) / sqrt(2 * oftype(float(x), π))

# Black's formula: the prices `(call, put)` of European options on an asset whose forward price at
# expiry is `F` and whose strike is `K`, both discounted to today, where `v ≥ 0` is the volatility of
# the log forward price to expiry (σ√τ under Black–Scholes). The option is at the money when F = K:
# the forward price equals the strike, with discounting and dividends included in both.
#
# Derivatives are Black's Greeks, not the formula differentiated: ∂/∂F is Φ(d₊) for the call and
# -Φ(-d₊) for the put, ∂/∂K is -Φ(d₋) and Φ(-d₋), and ∂/∂v is F·φ(d₊) for both. Differentiating
# the formula itself cancels terms of size m/v² (m the log-moneyness), which loses every digit as
# v → 0 and overflows to NaN below about v = 1e-155. Each derivative term is a strong-zero product
# (`__strong_zero_mul`), so a zero Greek contributes nothing where an input's derivative is infinite
# (√τ at τ = 0).
#
# Prices are accurate to a few roundings of F and K in absolute terms, not relative to the price: near
# the money F·Φ(d₊) and K·Φ(d₋) are each about F/2, so a time value of about F·v/√(2π) at a tiny v
# loses digits, half of them by v = 1e-8, and can round to zero below v ≈ 1e-16.
#
# Zero volatility is the intrinsic value, and the derivatives there are the limits as v → 0⁺ at
# fixed F and K: in the money Φ = 1 and no vega, out of the money Φ = 0 and no vega, and at the
# money Φ = 1/2 and vega F/√(2π), as the price grows like F·v/√(2π) there. The price is not jointly
# smooth in v and the moneyness at v = 0, so these are one-sided: a forward a rounding error away
# from the strike is in or out of the money, with no vega. Second derivatives are the Greeks
# differentiated; at v = 0 they are the limits (zero) away from the money and, at the money, in v
# alone, while gamma there is infinite (it grows like 1/v) and second derivatives that move the
# moneyness are not finite.
__black(F, K, v) = __black(promote(float(F), float(K), float(v))...)
function __black(F::T, K::T, v::T) where {T <: Real}
    # a negative volatility is not |v|: d₊ and d₋ change places and an at-the-money call is negative
    __primal(v) < 0 && throw(DomainError(v, "Black's formula needs a non-negative volatility"))
    d₊, d₋ = __black_d(F, K, v)
    return F * __black_Φ(d₊) - K * __black_Φ(d₋), K * __black_Φ(-d₋) - F * __black_Φ(-d₊)
end
function __black(F::D, K::D, v::D) where {Tag, D <: ForwardDiff.Dual{Tag}}
    f, k, s = ForwardDiff.value(F), ForwardDiff.value(K), ForwardDiff.value(v)
    call, put = __black(f, k, s)
    d₊, d₋ = __black_d(f, k, s)
    Φ₊, Φ₋, vega = __black_Φ(d₊), __black_Φ(d₋), f * __black_φ(d₊)
    Φ₊ᶜ, Φ₋ᶜ = __black_Φ(-d₊), __black_Φ(-d₋)
    ∂(w, p) = __strong_zero_mul(w, p)
    pF, pK, pv = ForwardDiff.partials(F).values, ForwardDiff.partials(K).values, ForwardDiff.partials(v).values
    ∂call = map((dF, dK, dv) -> ∂(Φ₊, dF) - ∂(Φ₋, dK) + ∂(vega, dv), pF, pK, pv)
    ∂put = map((dF, dK, dv) -> ∂(Φ₋ᶜ, dK) - ∂(Φ₊ᶜ, dF) + ∂(vega, dv), pF, pK, pv)
    return ForwardDiff.Dual{Tag}(call, ForwardDiff.Partials(∂call)), ForwardDiff.Dual{Tag}(put, ForwardDiff.Partials(∂put))
end

# d₊ = m/v + v/2 and d₋ = d₊ - v for the log-moneyness m = log(F/K). Within a factor of 2 of the money
# F - K is exact (Sterbenz), so log1p((F - K)/K) has the relative accuracy of m itself, and d keeps
# its digits where v is tiny; log(F/K) has an error of a rounding of 1. At v = 0, d₊ = d₋ is ±∞ by the
# sign of m, and 0 at the money: the ratio m/v of an exact zero m is zero (`__ratio`).
function __black_d(F, K, v)
    m = __primal(K) / 2 <= __primal(F) <= 2 * __primal(K) ? log1p((F - K) / K) : log(F / K)
    d₊ = __ratio(m, v) + v / 2
    return d₊, d₊ - v
end

# Φ(d) and φ(d) for Black's formula. A dual number's partials, φ(d)·d′ and -d·φ(d)·d′, are strong-zero
# products in φ(d): where φ(d) underflows to zero (|d| > 38.6, including the infinite d at v = 0 away
# from the money) a term is dropped even when d′ ≈ -m/v² has overflowed, where differentiating Φ or φ
# formed 0·∞ = NaN. A dropped term (before the vega's factor F) is below 4e-14 wherever d′ is finite,
# and away from the money, where d′ overflows, |d| exceeds 1e140 and the term is zero beyond any
# precision. The strong zero is in φ(d) only: -d·d′ is 0·∞ at the money at v = 0, where a second
# derivative that moves the moneyness is not finite, and at a subnormal v, where d′ = 1/(F·v) overflows.
__black_Φ(d::Real) = __N(d)
function __black_Φ(d::ForwardDiff.Dual{Tag}) where {Tag}
    x = ForwardDiff.value(d)
    w = __black_φ(x)
    return ForwardDiff.Dual{Tag}(__black_Φ(x), ForwardDiff.Partials(map(p -> __strong_zero_mul(w, p), ForwardDiff.partials(d).values)))
end
__black_φ(d::Real) = __normal_density(d)
function __black_φ(d::ForwardDiff.Dual{Tag}) where {Tag}
    x = ForwardDiff.value(d)
    w = __black_φ(x)
    return ForwardDiff.Dual{Tag}(w, ForwardDiff.Partials(map(p -> __strong_zero_mul(w, -x * p), ForwardDiff.partials(d).values)))
end

# Black–Scholes prices `(call, put)`: Black's formula on the discounted forward S·e^{-qτ} and strike
# K·e^{-rτ}, with the total volatility σ√τ. σ is checked itself, as σ√τ is zero at τ = 0 whatever its
# sign, where a negative σ gives a τ-derivative of the wrong sign. σ√τ is a strong-zero product: √τ
# has an infinite slope at τ = 0, which zero volatility doesn't multiply.
function __black_scholes(S, K, τ, r, σ, q)
    __primal(σ) < 0 && throw(DomainError(σ, "the volatility σ must be non-negative"))
    r, q = __continuous(r), __continuous(q)
    return __black(S * exp(-q * τ), K * exp(-r * τ), __strong_zero_mul(σ, __strong_zero_sqrt(τ)))
end

"""
    eurocall(;S=1.,K=1.,τ=1,r,σ,q=0.)

Calculate the Black-Scholes implied option price for a european call, where:

- `S` is the current asset price
- `K` is the strike or exercise price
- `τ` is the time remaining to maturity (can be typed with \\tau[tab])
- `r` is the continuously compounded risk free rate
- `σ` is the (implied) volatility (can be typed with \\sigma[tab])
- `q` is the continuously paid dividend rate

`r` and `q` are numbers, read as continuously compounded, or `Rate`s, converted to their
continuously compounded values.

Rates should be input as rates (not percentages), e.g.: `0.05` instead of `5` for a rate of five percent.

A zero volatility or time to maturity gives the intrinsic value of the forward,
`max(S e^{-qτ} - K e^{-rτ}, 0)`, and a negative `σ` throws a `DomainError`. ForwardDiff derivatives are
the Black–Scholes Greeks, accurate and finite down to `σ = 0`, where they are their limits as
`σ → 0⁺`: at the money, when the forward price `S e^{(r-q)τ}` equals `K`, the delta is `e^{-qτ}/2` and
the vega `S e^{-qτ}√τ/√(2π)`; otherwise the delta is `e^{-qτ}` or 0 and the vega 0. The price is not
smooth in `σ` and the moneyness together at `σ = 0`, so these derivatives are one-sided, and a forward
a rounding error away from `K` is in or out of the money.

Prices are accurate to a few roundings of `S e^{-qτ}` and `K e^{-rτ}`, not relative to the price
itself: near the money at a tiny volatility an option is worth about `S e^{-qτ}σ√τ/√(2π)`, which has
lost half its digits by `σ√τ = 1e-8` and can round to zero below `σ√τ ≈ 1e-16`.

!!! warning "Experimental"
    This function is well-tested, but the derivatives functionality (API) may change in a future version of FinanceModels.

# Extended Help

This is the same as the formulation presented in the [dividend extension of the BS model in Wikipedia](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model#Black%E2%80%93Scholes_equation).

## Other general comments:

- Swap/OIS curves are generally better sources for `r` than government debt (e.g. US Treasury) due to the collateralized nature of swap instruments.
- (Implied) volatility is characterized by a curve that is a function of the strike price (among other things), so take care when using
- A `Rate` in another convention is converted for you, e.g. `r = Periodic(0.05, 1)`.

"""
eurocall(; S = 1.0, K = 1.0, τ = 1, r, σ, q = 0.0) = first(__black_scholes(S, K, τ, r, σ, q))

"""
    europut(;S=1.,K=1.,τ=1,r,σ,q=0.)

Calculate the Black-Scholes implied option price for a european put, where:

- `S` is the current asset price
- `K` is the strike or exercise price
- `τ` is the time remaining to maturity (can be typed with \\tau[tab])
- `r` is the continuously compounded risk free rate
- `σ` is the (implied) volatility (can be typed with \\sigma[tab])
- `q` is the continuously paid dividend rate

`r` and `q` are numbers, read as continuously compounded, or `Rate`s, converted to their
continuously compounded values.

Rates should be input as rates (not percentages), e.g.: `0.05` instead of `5` for a rate of five percent.

A zero volatility or time to maturity gives the intrinsic value of the forward,
`max(K e^{-rτ} - S e^{-qτ}, 0)`, and derivatives at `σ = 0` are the one-sided limits described for
[`eurocall`](@ref); the delta at the money is `-e^{-qτ}/2`. Prices have the absolute accuracy described
there.

!!! warning "Experimental"
    This function is well-tested, but the derivatives functionality (API) may change in a future version of FinanceModels.

# Extended Help

This is the same as the formulation presented in the [dividend extension of the BS model in Wikipedia](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model#Black%E2%80%93Scholes_equation).

## Other general comments:

- Swap/OIS curves are generally better sources for `r` than government debt (e.g. US Treasury) due to the collateralized nature of swap instruments.
- (Implied) volatility is characterized by a curve that is a function of the strike price (among other things), so take care when using
- A `Rate` in another convention is converted for you, e.g. `r = Periodic(0.05, 1)`.

"""
europut(; S = 1.0, K = 1.0, τ = 1, r, σ, q = 0.0) = last(__black_scholes(S, K, τ, r, σ, q))

"""
    ReadOnlyVector(v::Vector)
    ReadOnlyVector{Rate{T, Continuous}, T}(v::Vector{T})

Internal read-only view over an owned `Vector`. It defines no `setindex!`, so indexed assignment
— and therefore `.=`, `fill!`, `sort!`, `reverse!`, and writes through `view` — throws Base's
`CanonicalIndexError`; reads behave as a normal `AbstractVector` (indexing, iteration,
`searchsortedlast`, broadcasting, `==`/`isequal`/`hash` identical to the equivalent `Vector`).
`copy`/`collect` return a mutable `Vector`. To change a curve's knots, use
`Accessors.@set curve.rates[i] = x` or `reconstruct`, which rebuild it.

The second form reads each stored number as a `Continuous` rate, without a copy: indexing builds
the rate, and `copy`/`collect` return a `Vector` of rates. [`knot_rates`](@ref) views a curve's
numeric knot rates this way.

Used by the knot curves (`Yield.Spline`, `Yield.MonotoneConvex`), which cache state derived
from their knot vectors, so that ordinary array operations on the public fields cannot
desynchronise the cache. This is Julia's conventional privacy, not literal immutability: the
backing `Vector` is the internal field `_data`, and code that mutates it is unsupported.
"""
struct ReadOnlyVector{T, S} <: AbstractVector{T}
    _data::Vector{S}
end
ReadOnlyVector(v::Vector{T}) where {T} = ReadOnlyVector{T, T}(v)

Base.size(v::ReadOnlyVector) = size(getfield(v, :_data))
Base.IndexStyle(::Type{<:ReadOnlyVector}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(v::ReadOnlyVector{T}, i::Int) where {T} = __element(T, getfield(v, :_data)[i])

# A stored value as an element of a `ReadOnlyVector` with element type `T`: the value itself, or a
# stored number read as a continuous rate.
__element(::Type, x) = x
__element(::Type{Rate{T, Continuous}}, x) where {T} = Rate{T, Continuous}(x, Continuous())

# The continuously compounded value of a rate input: a number is taken as continuous, and a `Rate`
# is converted from its own convention.
__continuous(x::Real) = x
__continuous(x::Rate) = FinanceCore.rate(convert(Continuous(), x))

# Every `frequency` input, a `Periodic` or an integer number of payments per period, as a `Periodic`.
__frequency(f::Periodic) = f
__frequency(n::Integer) = Periodic(n)

# The nominal rate of a coupon, margin, or interest-rate strike on a contract paid `frequency`
# times per period: a number as it is, or a `Periodic` rate of that frequency. Reading the nominal
# rate of another frequency and converting it give different coupons, so a mismatch throws.
__nominal(x::Real, frequency) = x
function __nominal(x::Rate{<:Any, Periodic}, frequency)
    n = __frequency(frequency).frequency
    x.compounding.frequency == n || throw(
        ArgumentError(
            "rate $x has frequency $(x.compounding.frequency), but the contract's frequency is $n; convert it explicitly."
        )
    )
    return FinanceCore.rate(x)
end

# The polynomial with exact (Rational) coefficients `cs` at `x`, each coefficient rounded once to
# the precision of `x` (for a dual number, of its primal value): Float32 stays Float32 and BigFloat
# keeps its precision. Inlined, so that for a concrete float type the rounding of a constant table
# folds at compile time.
@inline __evalpoly_exact(x, cs) = evalpoly(x, map(c -> convert(typeof(float(__primal(x))), c), cs))

# x·y where an exact zero factor gives zero even when the other is infinite or NaN ("strong zero"),
# applied to each product of the product rule, so a dual number's partials of every order follow the
# same rule: a term whose weight is zero contributes nothing when its other factor, or a derivative
# of it, has overflowed. Other products are x·y.
__strong_zero_mul(x, y) = __strong_zero_mul(promote(x, y)...)
__strong_zero_mul(x::T, y::T) where {T <: Real} = iszero(x) || iszero(y) ? zero(T) : x * y
function __strong_zero_mul(x::D, y::D) where {D <: ForwardDiff.Dual}
    vx, vy = ForwardDiff.value(x), ForwardDiff.value(y)
    p = map((dx, dy) -> __strong_zero_mul(vx, dy) + __strong_zero_mul(dx, vy), ForwardDiff.partials(x).values, ForwardDiff.partials(y).values)
    return D(__strong_zero_mul(vx, vy), ForwardDiff.Partials(p))
end

# p/v without forming 1/v, so that a ratio of two tiny numbers (1e-310/1e-310) is not ∞·tiny; 0 where
# p is 0. A dual number's partials are (∂p - r·∂v)/v with r = p/v, by the same rule, so no v² is formed.
__ratio(p, v) = __ratio(promote(p, v)...)
__ratio(p::T, v::T) where {T <: Real} = iszero(p) ? zero(T) : p / v
function __ratio(p::D, v::D) where {D <: ForwardDiff.Dual}
    p0, v0 = ForwardDiff.value(p), ForwardDiff.value(v)
    r = __ratio(p0, v0)
    q = map((dp, dv) -> __ratio(dp - __strong_zero_mul(r, dv), v0), ForwardDiff.partials(p).values, ForwardDiff.partials(v).values)
    return D(r, ForwardDiff.Partials(q))
end

# log(x), with partials ∂x/x as `__ratio`: a zero partial stays zero where 1/x overflows (a constant
# weight of 1e-310 under a time derivative), and a tiny one keeps its ratio to x.
__strong_zero_log(x::Real) = log(x)
function __strong_zero_log(x::D) where {D <: ForwardDiff.Dual}
    v = ForwardDiff.value(x)
    return D(__strong_zero_log(v), ForwardDiff.Partials(map(p -> __ratio(p, v), ForwardDiff.partials(x).values)))
end

# √x, with partials ∂x/(2√x) as `__ratio`: a zero partial stays zero at x = 0, where the slope 1/(2√x)
# is infinite and `sqrt` gives 0·∞ = NaN. A state that is exactly zero with zero partials (the CIR
# diffusion's clipped state) contributes nothing to a derivative of any order; a nonzero partial at
# x = 0 still gives ±Inf, the infinite slope. The result keeps the input's tag and takes the numeric
# type `sqrt` gives, so an integer or rational dual number promotes to floating point.
__strong_zero_sqrt(x::Real) = sqrt(x)
function __strong_zero_sqrt(x::ForwardDiff.Dual{Tag}) where {Tag}
    s = __strong_zero_sqrt(ForwardDiff.value(x))
    return ForwardDiff.Dual{Tag}(s, ForwardDiff.Partials(map(p -> __ratio(p, 2s), ForwardDiff.partials(x).values)))
end

# w·eˣ, where eˣ alone, or a derivative of either factor, can overflow while the product doesn't. The
# value is the product where that is a normal float, and sign(w)·exp(log|w| + x) otherwise; 0 for
# w = 0. A dual number's partials are eˣ·∂w + (w·eˣ)·∂x, the first again as `__wexp` and the second a
# strong-zero product with the value, so that each term is formed from factors that are finite
# wherever the term is representable, at every order.
__wexp(w, x) = __wexp(promote(w, x)...)
function __wexp(w::T, x::T) where {T <: Real}
    iszero(w) && return zero(T)
    p = w * exp(x)
    return floatmin(T) <= abs(p) < Inf ? p : sign(w) * exp(log(abs(w)) + x)
end
function __wexp(w::D, x::D) where {D <: ForwardDiff.Dual}
    w0, x0 = ForwardDiff.value(w), ForwardDiff.value(x)
    v = __wexp(w0, x0)
    q = map((dw, dx) -> __wexp(dw, x0) + __strong_zero_mul(v, dx), ForwardDiff.partials(w).values, ForwardDiff.partials(x).values)
    return D(v, ForwardDiff.Partials(q))
end

# The float element type of a collection of reals, for a copy or an accumulator that keeps the values'
# numeric type (BigFloat, dual numbers): the declared element type when concrete, otherwise the
# promotion of the values' types. An untyped empty collection gives `Float64`: the promotion starts
# from `Bool`, which every real type absorbs.
__float_eltype(v) = float(isconcretetype(eltype(v)) ? eltype(v) : mapreduce(typeof, promote_type, v; init = Bool))
