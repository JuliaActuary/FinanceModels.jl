# Assets with an exercise right. Included in `module Bond`.

"""
    Bond.Callable(asset, times, prices)

`asset` with the issuer's right to redeem it at each of `times`, immediately after the payment
due then, for `prices[i]` times the principal remaining after that payment (a number applies to
every time). The issuer calls when that is worth less than continuing, so the holder's value is
never above `asset`'s.

`asset` is a `Bond.Fixed`, `Bond.Amortizing` or `Bond.Pool`, possibly a
[`Credit.Defaultable`](@ref FinanceModels.Credit.Defaultable) one; `Credit.Defaultable` may also
wrap the callable asset, with the same meaning. Each exercise time must be a payment date of
`asset` before its maturity, and there is no notice period. A pool is redeemed after its scheduled
payment and prepayment.

An exercisable asset is valued on a lattice: its valuation model must be a
[`ShortRate.HullWhiteTree`](@ref FinanceModels.ShortRate.HullWhiteTree), and a plain curve throws.
First derivatives follow the exercise decisions the lattice makes. Second derivatives through an
exercise decision throw: the second derivative of the lattice's value along its decisions leaves
out the effect of the exercise boundary moving, which is where a callable bond's negative
convexity comes from. Measure convexity by revaluing at bumped inputs.

```julia
bond = Bond.Callable(Bond.Fixed(0.05, Periodic(2), 10), 3.0:0.5:9.5, 1.0)   # callable at par after 3 years
present_value(ShortRate.HullWhiteTree(hw), bond)
```

See also [`Bond.Putable`](@ref), [`exercise_policy`](@ref FinanceModels.exercise_policy).
"""
struct Callable{A <: AbstractContract, T <: Real, P <: Real} <: AbstractBond
    asset::A
    times::Vector{T}
    prices::Vector{P}
    function Callable(asset::A, times, prices) where {A <: AbstractContract}
        ts, ps = __exercise_schedule(times, prices)
        return new{A, eltype(ts), eltype(ps)}(asset, ts, ps)
    end
end

"""
    Bond.Putable(asset, times, prices)

`asset` with the holder's right to sell it back at each of `times`, immediately after the payment
due then, for `prices[i]` times the principal remaining after that payment (a number applies to
every time). The holder puts when that is worth more than continuing, so the value is never below
`asset`'s. Otherwise as [`Bond.Callable`](@ref).
"""
struct Putable{A <: AbstractContract, T <: Real, P <: Real} <: AbstractBond
    asset::A
    times::Vector{T}
    prices::Vector{P}
    function Putable(asset::A, times, prices) where {A <: AbstractContract}
        ts, ps = __exercise_schedule(times, prices)
        return new{A, eltype(ts), eltype(ps)}(asset, ts, ps)
    end
end

# Validates and copies an exercise schedule: positive, strictly increasing times, one price each
# (a number for all).
function __exercise_schedule(times, prices)
    ts = collect(times)
    ps = prices isa Real ? fill(prices, length(ts)) : collect(prices)
    isempty(ts) && throw(ArgumentError("an exercise schedule needs at least one exercise time"))
    length(ps) == length(ts) ||
        throw(DimensionMismatch("an exercise schedule needs one price per time; got $(length(ts)) times and $(length(ps)) prices"))
    for (i, t) in enumerate(ts)
        p, prev = __primal(t), i == 1 ? zero(__primal(t)) : __primal(ts[i - 1])
        (isfinite(p) && p > prev) || throw(ArgumentError("exercise times must be finite, positive and strictly increasing; got $ts"))
    end
    return ts, ps
end

maturity(c::Union{Callable, Putable}) = maturity(c.asset)
