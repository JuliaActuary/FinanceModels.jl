"""
Default intensities and the default assumptions that defaultable contracts read by key.

A hazard curve gives the intensity of default `λ(t)`, the rate of default at time `t` of an issuer
that has survived to `t`. Its integral is the cumulative hazard `Λ`, survival from `from` to `to`
is `exp(-(Λ(to) - Λ(from)))`, and the default probability is the complement. Two curves are
provided:

- [`Credit.ConstantHazard`](@ref): one intensity at every time;
- [`Credit.PiecewiseHazard`](@ref): constant between tenors, flat beyond the last.

Intensities are nonnegative. A [`Credit.Assumptions`](@ref) pairs a hazard curve with a recovery
fraction, and is the entry a defaultable contract reads from a [`Models`](@ref FinanceModels.Models)
store.

```julia
h = Credit.PiecewiseHazard([0.01, 0.02, 0.03], [1.0, 5.0, 10.0])
Credit.survival(h, 7.0)                    # exp(-(0.01 + 0.08 + 0.06))
Credit.default_probability(h, 7.0, 8.0)    # 1 - exp(-0.03)
```

The functions share their names with MortalityTables' (`hazard`, `cumhazard`, `survival`) but are
FinanceModels' own; call them qualified, `Credit.survival`, when both packages are loaded.
"""
module Credit
import ..FinanceCore
import ..AbstractModel
import ..__primal
import ..ReadOnlyVector
import ..Yield: Unchecked

export ConstantHazard, PiecewiseHazard, Assumptions, Defaultable

"""
    Credit.AbstractHazard

The supertype of hazard curves. A curve defines [`Credit.cumhazard`](@ref)`(c, from, to)` and
[`Credit.hazard`](@ref)`(c, t)`; the other functions derive from them.
"""
abstract type AbstractHazard <: AbstractModel end

# A hazard rate must be finite and nonnegative: a negative one gives survival above one.
function __check_hazard(λ)
    p = __primal(λ)
    (isfinite(p) && p >= 0) || throw(DomainError(λ, "a hazard rate must be finite and nonnegative"))
    return λ
end

function __check_time(t)
    __primal(t) >= 0 || throw(DomainError(t, "hazard curves start at time 0; got a negative time"))
    return t
end

"""
    Credit.ConstantHazard(λ)

A constant default intensity `λ ≥ 0`: survival to `t` is `exp(-λt)`.

```julia
h = Credit.ConstantHazard(0.02)
Credit.survival(h, 5.0)              # exp(-0.1)
Credit.default_probability(h, 1.0)   # 1 - exp(-0.02)
```
"""
struct ConstantHazard{T <: Real} <: AbstractHazard
    rate::T
    ConstantHazard(λ::T) where {T <: Real} = new{T}(__check_hazard(λ))
end

"""
    Credit.PiecewiseHazard(rates, tenors)

A default intensity that is `rates[i]` from `tenors[i-1]` to `tenors[i]` (from 0 for the
first), and `rates[end]` beyond the last tenor. The intensity is right-continuous: at a tenor
it is the next interval's. Tenors are positive and strictly increasing, rates finite and
nonnegative; both vectors are copied.

`Credit.PiecewiseHazard()` has no tenors: it is the template a calibration fills in (see
[`fit`](@ref FinanceModels.fit)), and evaluating it throws.

```julia
h = Credit.PiecewiseHazard([0.01, 0.02], [2.0, 5.0])
Credit.hazard(h, 2.0)       # 0.02
Credit.cumhazard(h, 6.0)    # 0.01·2 + 0.02·4
```
"""
struct PiecewiseHazard{R <: Real, T <: Real, C <: Real} <: AbstractHazard
    rates::ReadOnlyVector{R, R}
    tenors::ReadOnlyVector{T, T}
    _cumhazard::Vector{C}   # Λ at each tenor
    function PiecewiseHazard(rates::AbstractVector{<:Real}, tenors::AbstractVector{<:Real})
        n = length(tenors)
        length(rates) == n ||
            throw(DimensionMismatch("a PiecewiseHazard needs one rate per tenor, got $(length(rates)) rates for $n tenors"))
        n >= 1 || throw(ArgumentError("a PiecewiseHazard needs at least one tenor"))
        ts, λs = collect(tenors), collect(rates)
        foreach(__check_hazard, λs)
        for (i, t) in enumerate(ts)
            p, prev = __primal(t), i == 1 ? zero(__primal(t)) : __primal(ts[i - 1])
            (isfinite(p) && p > prev) || throw(
                ArgumentError("PiecewiseHazard tenors must be finite, positive and strictly increasing; got $ts")
            )
        end
        C = __cumulative(λs, ts)
        return new{eltype(λs), eltype(ts), eltype(C)}(ReadOnlyVector(λs), ReadOnlyVector(ts), C)
    end
    # The template of a calibration (see `fit`): no tenors yet, so it cannot be evaluated.
    PiecewiseHazard() = new{Float64, Float64, Float64}(ReadOnlyVector(Float64[]), ReadOnlyVector(Float64[]), Float64[])
    # A calibration's trial curve: neither copied nor validated, so a solver may try a negative intensity.
    function PiecewiseHazard(::Unchecked, rates::Vector{R}, tenors::Vector{T}) where {R <: Real, T <: Real}
        C = __cumulative(rates, tenors)
        return new{R, T, eltype(C)}(ReadOnlyVector(rates), ReadOnlyVector(tenors), C)
    end
end

# Λ at each tenor, accumulated interval by interval.
function __cumulative(λs, ts)
    C = Vector{typeof(first(λs) * first(ts))}(undef, length(ts))
    acc = zero(eltype(C))
    for i in eachindex(ts)
        acc += λs[i] * (ts[i] - (i == 1 ? zero(ts[i]) : ts[i - 1]))
        C[i] = acc
    end
    return C
end

"""
    Credit.hazard(c, t)

The default intensity of `c` at time `t`, right-continuous at a [`Credit.PiecewiseHazard`](@ref)'s tenors.
"""
hazard(c::ConstantHazard, t) = (__check_time(t); c.rate)
function hazard(c::PiecewiseHazard, t)
    __check_time(t)
    return c.rates[min(searchsortedlast(c.tenors, t) + 1, length(c.rates))]
end

"""
    Credit.cumhazard(c, t)
    Credit.cumhazard(c, from, to)

The integrated default intensity of `c` from 0 (or `from`) to `t` (or `to`). The interval form
integrates directly, so it keeps its relative precision over a short interval at a late time.
At `t = Inf` it is `Inf`, unless the last intensity is zero.
"""
cumhazard(c::ConstantHazard, t) = (__check_time(t); __flat_integral(c.rate, t))
cumhazard(c::ConstantHazard, from, to) = (__check_time(from); __check_time(to); __flat_integral(c.rate, to - from))
function cumhazard(c::PiecewiseHazard, t)
    __check_time(t)
    ts, λs, n = c.tenors, c.rates, length(c.rates)
    k = searchsortedlast(ts, t)
    k == 0 && return λs[1] * t
    k >= n && return c._cumhazard[n] + __flat_integral(λs[n], t - ts[n])
    return c._cumhazard[k] + λs[k + 1] * (t - ts[k])
end
function cumhazard(c::PiecewiseHazard, from, to)
    __check_time(from)
    __check_time(to)
    to < from && return -cumhazard(c, to, from)
    ts, λs, n = c.tenors, c.rates, length(c.rates)
    # the interval of `from`, i.e. the number of tenors at or before it
    k = searchsortedlast(ts, from)
    k >= n && return __flat_integral(λs[n], to - from)
    total = zero(first(λs) * (to - from))
    lo = from
    for i in (k + 1):n
        hi = min(to, ts[i])
        total += λs[i] * (hi - lo)
        to <= ts[i] && return total
        lo = ts[i]
    end
    return total + __flat_integral(λs[n], to - ts[n])
end
# λ·Δ, with a zero intensity over an infinite interval integrating to zero
__flat_integral(λ, Δ) = iszero(λ) ? zero(λ * one(Δ)) : λ * Δ

"""
    Credit.survival(c, t)
    Credit.survival(c, from, to)

The probability of no default from 0 (or `from`) to `t` (or `to`): `exp(-cumhazard)`.
"""
survival(c::AbstractHazard, t) = exp(-cumhazard(c, t))
survival(c::AbstractHazard, from, to) = exp(-cumhazard(c, from, to))

"""
    Credit.default_probability(c, t)
    Credit.default_probability(c, from, to)

The probability of default from 0 (or `from`) to `t` (or `to`), given survival to its start:
`1 - survival`, computed as `-expm1(-cumhazard)` so that small probabilities keep their digits.
"""
default_probability(c::AbstractHazard, t) = -expm1(-cumhazard(c, t))
default_probability(c::AbstractHazard, from, to) = -expm1(-cumhazard(c, from, to))

# The default time of an issuer whose cumulative hazard reaches `E`: Λ⁻¹(E), the first time Λ = E,
# and `Inf` when it never does. With `E` an exponential draw, it samples the default time.
function __default_time(c::ConstantHazard, E)
    return iszero(c.rate) ? oftype(float(E / one(c.rate)), Inf) : E / c.rate
end
function __default_time(c::PiecewiseHazard, E)
    ts, λs, C = c.tenors, c.rates, c._cumhazard
    i = searchsortedfirst(C, E)   # the first interval whose end reaches E
    if i > length(ts)
        λ = λs[end]
        return iszero(λ) ? oftype(float(E), Inf) : ts[end] + (E - C[end]) / λ
    end
    start, base = i == 1 ? (zero(ts[1]), zero(C[1])) : (ts[i - 1], C[i - 1])
    return start + (E - base) / λs[i]
end

# The times at which the intensity can jump.
__breakpoints(::ConstantHazard) = ()
__breakpoints(c::PiecewiseHazard) = c.tenors

"""
    Credit.Assumptions(hazard; recovery)

The default assumptions a defaultable contract reads from a valuation context under its key:
the issuer's `hazard` curve, and the `recovery` fraction of the principal outstanding
immediately before default, paid at default. Recovery excludes accrued interest, and no
contractual payment follows a default. `recovery` lies in `[0, 1]` and has no default.

```julia
acme = Credit.Assumptions(Credit.ConstantHazard(0.02); recovery = 0.4)
ctx = Models(curve, Dict("ACME" => acme))
```
"""
struct Assumptions{H <: AbstractHazard, R <: Real} <: AbstractModel
    hazard::H
    recovery::R
    function Assumptions(hazard::H, recovery::R) where {H <: AbstractHazard, R <: Real}
        p = __primal(recovery)
        0 <= p <= 1 || throw(DomainError(recovery, "a recovery fraction lies in [0, 1]"))
        return new{H, R}(hazard, recovery)
    end
end
Assumptions(hazard; recovery) = Assumptions(hazard, recovery)

"""
    Credit.Defaultable(asset, key)

`asset` subject to the default of its issuer, whose [`Credit.Assumptions`](@ref) a valuation
context holds under `key`. `asset` has principal: a `Bond.Fixed`, `Bond.Amortizing` or
`Bond.Pool`. On default, at any time between payments, the holder recovers the `recovery`
fraction of the principal outstanding immediately before it, and no later payment is made.

Its cashflow projection is its expected cashflows: each payment weighted by survival to it, and
the expected recoveries between payments at quadrature nodes (4-point Gauss–Legendre on each
interval between payments and the hazard curve's tenors, in parts no longer than a month). Its
value is the value of those cashflows, so it discounts on the context's curve and composes with
portfolios. The nodes depend only on the times, so a value with dual numbers equals the value
without them.

`key` identifies the issuer: contracts with the same key default together in sampled projections.

```julia
bond = Credit.Defaultable(Bond.Fixed(0.05, Periodic(2), 10), "ACME")
ctx = Models(Yield.Constant(0.04), Dict("ACME" => Credit.Assumptions(Credit.ConstantHazard(0.02); recovery = 0.4)))
present_value(ctx, bond)
collect(Projection(bond, ctx))   # expected payments and recoveries
```
"""
struct Defaultable{A <: FinanceCore.AbstractContract, K} <: FinanceCore.AbstractContract
    asset::A
    key::K
end
FinanceCore.maturity(d::Defaultable) = FinanceCore.maturity(d.asset)

end
