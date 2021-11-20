module Yields

import Interpolations
import ForwardDiff
using LinearAlgebra
using UnicodePlots

# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCurve
export rate, discount, accumulation, forward, Yield, Rate, rate, spot

abstract type CompoundingFrequency end
Base.Broadcast.broadcastable(x::T) where {T<:CompoundingFrequency} = Ref(x)

""" 
    Continuous()

A type representing continuous interest compounding frequency.

# Examples

```julia-repl
julia> Rate(0.01,Continuous())
Rate(0.01, Continuous())
```

See also: [`Periodic`](@ref)
"""
struct Continuous <: CompoundingFrequency end


""" 


```julia-repl
julia> Continuous(0.01)
Rate(0.01, Continuous())
```

See also: [`Periodic`](@ref)
"""
Continuous(rate) = Rate(rate, Continuous())

""" 
    Periodic(frequency)

A type representing periodic interest compounding with the given frequency

# Examples

Creating a semi-annual bond equivalent yield:

```julia-repl
julia> Rate(0.01,Periodic(2))
Rate(0.01, Periodic(2))
```

See also: [`Continuous`](@ref)
"""
struct Periodic <: CompoundingFrequency
    frequency::Int
end

""" 
    Periodic(rate,frequency)

A convenience constructor for Rate(rate,Periodic(frequency)).

# Examples

Creating a semi-annual bond equivalent yield:

```julia-repl
julia> Periodic(0.01,2)
Rate(0.01, Periodic(2))
```

See also: [`Continuous`](@ref)
"""
Periodic(x, frequency) = Rate(x, Periodic(frequency))

struct Rate{N<:Real,T<:CompoundingFrequency}
    value::N
    compounding::T
end

# Base.:==(r1::Rate,r2::Rate) = (r1.value == r2.value) && (r1.compounding == r2.compounding)

"""
    Rate(rate[,frequency=1])
    Rate(rate,frequency::CompoundingFrequency)

Rate is a type that encapsulates an interest `rate` along with its compounding `frequency`.

Periodic rates can be constructed via `Rate(rate,frequency)` or `Rate(rate,Periodic(frequency))`.

Continuous rates can be constructed via `Rate(rate, Inf)` or `Rate(rate,Continuous())`.

# Examples

```julia-repl
julia> Rate(0.01,Continuous())
Rate(0.01, Continuous())

julia> Rate(0.01,Periodic(2))
Rate(0.01, Periodic(2))

julia> Rate(0.01)
Rate(0.01, Periodic(1))

julia> Rate(0.01,2)
Rate(0.01, Periodic(2))

julia> Rate(0.01,Periodic(4))
Rate(0.01, Periodic(4))

julia> Rate(0.01,Inf)
Rate(0.01, Continuous())

julia> Rate(0.01,Continuous())
Rate(0.01, Continuous())
```
"""
Rate(rate) = Rate(rate, Periodic(1))
Rate(x, frequency::T) where {T<:Real} = isinf(frequency) ? Rate(x, Continuous()) : Rate(x, Periodic(frequency))

"""
    convert(T::CompoundingFrequency,r::Rate)

Returns a `Rate` with an equivalent discount but represented with a different compounding frequency.

# Examples

```julia-repl
julia> r = Rate(Periodic(12),0.01)
Rate(0.01, Periodic(12))

julia> convert(Periodic(1),r)
Rate(0.010045960887181016, Periodic(1))

julia> convert(Continuous(),r)
Rate(0.009995835646701251, Continuous())
```
"""
function Base.convert(T::CompoundingFrequency, r::Rate{<:Real,<:CompoundingFrequency})
    convert(T, r, r.compounding)
end
function Base.convert(to::Continuous, r, from::Continuous)
    return r
end

function Base.convert(to::Periodic, r, from::Continuous)
    return Rate(to.frequency * (exp(r.value / to.frequency) - 1), to)
end

function Base.convert(to::Continuous, r, from::Periodic)
    return Rate(from.frequency * log(1 + r.value / from.frequency), to)
end

function Base.convert(to::Periodic, r, from::Periodic)
    c = convert(Continuous(), r, from)
    return convert(to, c, Continuous())
end

function rate(r::Rate{<:Real,<:CompoundingFrequency})
    r.value
end

function Base.isapprox(a::Rate{N,T}, b::Rate{N,T}; atol::Real = 0, rtol::Real = atol > 0 ? 0 : √eps()) where {T<:Periodic,N<:Real}
    return (a.compounding.frequency == b.compounding.frequency) && isapprox(rate(a), rate(b); atol, rtol)
end

function Base.isapprox(a::Rate{N,T}, b::Rate{N,T}; atol::Real = 0, rtol::Real = atol > 0 ? 0 : √eps()) where {T<:Continuous,N<:Real}
    return isapprox(rate(a), rate(b); atol, rtol)
end

# the fallback for rates not of the same type
function Base.isapprox(a::T, b::N; atol::Real = 0, rtol::Real = atol > 0 ? 0 : √eps()) where {T<:Rate,N<:Rate}
    return isapprox(convert(b.compounding, a), b; atol, rtol)
end

"""
An AbstractYield is an object which can be used as an argument to:

- zero-coupon spot rates via [`zero`](@ref)
- discount factor via [`discount`](@ref)
- accumulation factor via [`accumulation`](@ref)

It can be be constructed via:

- zero rate curve with [`Zero`](@ref)
- forward rate curve with [`Forward`](@ref)
- par rate curve with [`Par`](@ref)
- typical OIS curve with [`OIS`](@ref)
- typical constant maturity treasury (CMT) curve with [`CMT`](@ref)
"""
abstract type AbstractYield end

# make interest curve broadcastable so that you can broadcast over multiple`time`s in `interest_rate`
Base.Broadcast.broadcastable(ic::T) where {T<:AbstractYield} = Ref(ic)

struct YieldCurve{T,U,V} <: AbstractYield
    rates::T
    maturities::U
    zero::V # function time -> continuous zero rate
end

# Forward curves

"""
    ForwardStarting(curve,forwardstart)

Rebase a `curve` so that `discount`/`accumulation`/etc. are re-based so that time zero from the new curves perspective is the given `forwardstart` time.

# Examples

```julia-repl
julia> zero = [5.0, 5.8, 6.4, 6.8] ./ 100
julia> maturity = [0.5, 1.0, 1.5, 2.0]
julia> curve = Yields.Zero(zero, maturity)
julia> fwd = Yields.ForwardStarting(curve, 1.0)

julia> discount(curve,1,2)
0.9275624570410582

julia> discount(fwd,1) # `curve` has effectively been reindexed to `1.0`
0.9275624570410582
```

# Extended Help

While `ForwardStarting` could be nested so that, e.g. the third period's curve is the one-period forward of the second period's curve, it will be more efficient to reuse the initial curve from a runtime and compiler perspective.

`ForwardStarting` is not used to construct a curve based on forward rates. See  [`Forward`](@ref) instead.
"""
struct ForwardStarting{T,U} <: AbstractYield
    curve::U
    forwardstart::T
end

function discount(c::ForwardStarting, to)
    discount(c.curve, c.forwardstart, to + c.forwardstart)
end

"""
    zero(curve,time)
    zero(curve,time,CompoundingFrequency)

Return the zero rate for the curve at the given time. If not specified, will use `Periodic(1)` compounding.
"""
Base.zero(c::YieldCurve, time) = convert(Periodic(1), Continuous(c.zero(time)))
function Base.zero(c::YieldCurve, time, cf::Periodic)
    z = Continuous(c.zero(time))
    i = convert(cf, z)
    return i
end

function Base.zero(c::YieldCurve, time, cf::Continuous)
    return Continuous(c.zero(time))
end


"""
    Constant(rate)

Construct a yield object where the spot rate is constant for all maturities. If `rate` is not a `Rate` type, will assume `Periodic(1)` for the compounding frequency

# Examples

```julia-repl
julia> y = Yields.Constant(0.05)
julia> discount(y,2)
0.9070294784580498     # 1 / (1.05) ^ 2
```
"""
struct Constant{T} <: AbstractYield
    rate::T
end

function Constant(rate::T) where {T<:Real}
    return Constant(Rate(rate, Periodic(1)))
end

rate(c::Constant) = c.rate
rate(c::Constant, time) = c.rate
discount(c::T, time) where {T<:Real} = discount(Constant(c), time)
discount(r::Constant, time) = 1 / accumulation(r, time)

accumulation(r::Constant, time) = accumulation(r.rate.compounding, r, time)
accumulation(c::T, time) where {T>:Real} = accumulation(Constant(c), time)
accumulation(::Continuous, r::Constant, time) = exp(rate(r.rate) * time)
accumulation(::Periodic, r::Constant, time) = (1 + rate(r.rate) / r.rate.compounding.frequency)^(r.rate.compounding.frequency * time)

"""
    Step(rates,times)

Create a yield curve object where the applicable rate is the effective rate of interest applicable until corresponding time. If `rates` is not a `Vector{Rate}`, will assume `Periodic(1)` type.

# Examples

```julia-repl
julia>y = Yields.Step([0.02,0.05], [1,2])

julia>rate(y,0.5)
0.02

julia>rate(y,1.5)
0.05

julia>rate(y,2.5)
0.05
```
"""
struct Step <: AbstractYield
    rates
    times
end

Step(rates) = Step(rates, collect(1:length(rates)))

function rate(y::Step, time)
    i = findfirst(t -> time <= t, y.times)
    if isnothing(i)
        return y.rates[end]
    else
        return y.rates[i]
    end
end

function discount(y::Step, time)
    v = 1 / (1 + y.rates[1])^min(y.times[1], time)

    if y.times[1] >= time
        return v
    end

    for i = 2:length(y.times)

        if y.times[i] >= time
            # take partial discount and break
            v /= (1 + y.rates[i])^(time - y.times[i-1])
            break
        else
            # take full discount and continue
            v /= (1 + y.rates[i])^(y.times[i] - y.times[i-1])
        end

    end

    return v
end


"""
    Zero(rates,maturities)

Construct a yield curve with given zero-coupon spot `rates` at the given `maturities`. If `rates` is not a `Vector{Rate}`, will assume `Periodic(1)` type.
"""
function Zero(rates::Vector{<:Rate}, maturities)
    # bump to a constant yield if only given one rate
    length(rates) == 1 && return Constant(first(rates))

    continuous_zeros = rate.(convert.(Continuous(), rates))
    return YieldCurve(
        rates,
        maturities,
        linear_interp([0.0; maturities], [first(continuous_zeros); continuous_zeros])
    )
end

#fallback if `rates` aren't `Rate`s. Assume `Periodic(1)` per Zero docstring
function Zero(rates, maturities)
    Zero(Periodic.(rates, 1), maturities)
end

function Zero(rates)
    # bump to a constant yield if only given one rate
    maturities = collect(1:length(rates))
    return Zero(rates, maturities)
end

"""
    Par(rate, maturity)

Construct a curve given a set of bond equivalent yields and the corresponding maturities. Assumes that maturities <= 1 year do not pay coupons and that after one year, pays coupons with frequency equal to the CompoundingFrequency of the corresponding rate.

# Examples

```julia-repl

julia> par = [6.,8.,9.5,10.5,11.0,11.25,11.38,11.44,11.48,11.5] ./ 100
julia> maturities = [t for t in 1:10]
julia> curve = Par(par,maturities);
julia> zero(curve,1)
Rate(0.06000000000000005, Periodic(1))

```
"""
function Par(rates::Vector{<:Rate}, maturities)
    # bump to a constant yield if only given one rate
    if length(rates) == 1
        return Constant(rate[1])
    end
    return YieldCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise semi-annual
        # per Hull 4.7
        bootstrap(rates, maturities, [m <= 1 ? nothing : 1 / r.compounding.frequency for (r, m) in zip(rates, maturities)])
    )
end

function Par(rates::Vector{T}, maturities) where {T<:Real}
    return Par(Rate.(rates), maturities)
end


"""
    Forward(rate_vector,maturities)

Takes a vector of 1-period forward rates and constructs a discount curve.
"""
function Forward(rates, maturities)
    # convert to zeros and pass to Zero
    disc_v = similar(rates)
    v = 1.0

    for i = 1:length(rates)
        Δt = maturities[i] - (i == 1 ? 0 : maturities[i-1])
        v *= discount(Constant(rates[i]), Δt)
        disc_v[i] = v
    end

    z = (1.0 ./ disc_v) .^ (1 ./ maturities) .- 1 # convert disc_v to zero
    return Zero(z, maturities)
end

Forward(rates) = Forward(rates, collect(1:length(rates)))

"""
Takes CMT yields (bond equivalent), and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual.
"""
function CMT(rates::Vector{T}, maturities) where {T<:Real}
    rs = map(zip(rates, maturities)) do (r, m)
        if m <= 1
            Rate(r, Periodic(1 / m))
        else
            Rate(r, Periodic(2))
        end
    end

    CMT(rs, maturities)
end

function CMT(rates::Vector{<:Rate}, maturities)
    return YieldCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise semi-annual
        # per Hull 4.7
        bootstrap(rates, maturities, [m <= 1 ? nothing : 0.5 for m in maturities])
    )
end


"""
    OIS(rates,maturities)
Takes Overnight Index Swap rates, and assumes that instruments <= one year maturity are settled once and other agreements are settled quarterly with a corresponding CompoundingFrequency

"""
function OIS(rates::Vector{T}, maturities) where {T<:Real}
    rs = map(zip(rates, maturities)) do (r, m)
        if m <= 1
            Rate(r, Periodic(1 / m))
        else
            Rate(r, Periodic(4))
        end
    end

    return OIS(rs, maturities)
end
function OIS(rates::Vector{<:Rate}, maturities)
    return YieldCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise quarterly
        # per Hull 4.7
        bootstrap(rates, maturities, [m <= 1 ? nothing : 1 / 4 for m in maturities])
    )
end

abstract type ObservableQuote end

"""
    ZeroCouponQuote(price, maturity)

Quote for a set of zero coupon bonds with given `price` and `maturity`. 

# Examples

```julia-repl
julia> prices = [1.3, 0.1, 4.5]
julia> maturities = [1.2, 2.5, 3.6]
julia> swq = Yields.ZeroCouponQuote.(prices, maturities)
```
"""
struct ZeroCouponQuote <: ObservableQuote
    price
    maturity
end

"""
    SwapQuote(yield, maturity, frequency)

Quote for a set of interest rate swaps with the given `yield` and `maturity` and a given payment `frequency`.

# Examples

```julia-repl
julia> maturities = [1.2, 2.5, 3.6]
julia> interests = [-0.02, 0.3, 0.04]
julia> prices = [1.3, 0.1, 4.5]
julia> frequencies = [2,1,2]
julia> swq = Yields.SwapQuote.(interests, maturities, frequencies)
```
"""
struct SwapQuote <: ObservableQuote
    yield
    maturity
    frequency
    function SwapQuote(yield, maturity, frequency)
        frequency <= 0 && throw(DomainError("Payment frequency must be positive"))
        return new(yield, maturity, frequency)
    end
end


"""
    BulletBondQuote(yield, price, maturity, frequency)

Quote for a set of fixed interest bullet bonds with given `yield`, `price`, `maturity` and a given payment frequency `frequency`.

Construct a vector of quotes for use with SmithWilson methods, e.g. by broadcasting over an array of inputs.

# Examples

```julia-repl
julia> maturities = [1.2, 2.5, 3.6]
julia> interests = [-0.02, 0.3, 0.04]
julia> prices = [1.3, 0.1, 4.5]
julia> frequencies = [2,1,2]
julia> bbq = Yields.BulletBondQuote.(interests, maturities, prices, frequencies)
```
"""
struct BulletBondQuote <: ObservableQuote
    yield
    price
    maturity
    frequency

    function BulletBondQuote(yield, maturity, price, frequency)
        frequency <= 0 && throw(DomainError("Payment frequency must be positive"))
        return new(yield, maturity, price, frequency)
    end
end


"""
    SmithWilson(zcq::Vector{ZeroCouponQuote}; ufr, α)
    SmithWilson(swq::Vector{SwapQuote}; ufr, α)
    SmithWilson(bbq::Vector{BulletBondQuote}; ufr, α)
    SmithWilson(times<:AbstractVector, cashflows<:AbstractMatrix, prices<:AbstractVector; ufr, α)
    SmithWilson(u, qb; ufr, α)
    
Create a yield curve object that implements the Smith-Wilson interpolation/extrapolation scheme.

Positional arguments to construct a curve:

- Quoted instrument as the first argument: either a `Vector` of `ZeroCouponQuote`s, `SwapQuote`s, or `BulletBondQuote`s, or 
- A set of `times`, `cashflows`, and `prices`, or
- A curve can be with `u` is the timepoints coming from the calibration, and `qb` is the internal parameterization of the curve that ensures that the calibration is correct. Users may prefer the other constructors but this mathematical constructor is also available.

Required keyword arguments:

- `ufr` is the Ultimate Forward Rate, the forward interest rate to which the yield curve tends, in continuous compounding convention. 
- `α` is the parameter that governs the speed of convergence towards the Ultimate Forward Rate. It can be typed with `\\alpha[TAB]`
"""
struct SmithWilson{TU<:AbstractVector,TQb<:AbstractVector} <: AbstractYield
    u::TU
    qb::TQb
    ufr
    α

    # Inner constructor ensures that vector lengths match
    function SmithWilson{TU,TQb}(u, qb; ufr, α) where {TU<:AbstractVector,TQb<:AbstractVector}
        if length(u) != length(qb)
            throw(DomainError("Vectors u and qb in SmithWilson must have equal length"))
        end
        return new(u, qb, ufr, α)
    end
end

SmithWilson(u::TU, qb::TQb; ufr, α) where {TU<:AbstractVector,TQb<:AbstractVector} = SmithWilson{TU,TQb}(u, qb; ufr = ufr, α = α)

"""
    H_ordered(α, t_min, t_max)

The Smith-Wilson H function with ordered arguments (for better performance than using min and max).
"""
function H_ordered(α, t_min, t_max)
    return α * t_min + exp(-α * t_max) * sinh(-α * t_min)
end

"""
    H(α, t1, t2)

The Smith-Wilson H function implemented in a faster way.
"""
function H(α, t1::T, t2::T) where {T}
    return t1 < t2 ? H_ordered(α, t1, t2) : H_ordered(α, t2, t1)
end

H(α, t1, t2) = H(α, promote(t1, t2)...)

H(α, t1vec::AbstractVector, t2) = [H(α, t1, t2) for t1 in t1vec]
H(α, t1vec::AbstractVector, t2vec::AbstractVector) = [H(α, t1, t2) for t1 in t1vec, t2 in t2vec]
# This can be optimized by going to H_ordered directly, but it might be a bit cumbersome 
H(α, tvec::AbstractVector) = H(α, tvec, tvec)


discount(sw::SmithWilson, t) = exp(-sw.ufr * t) * (1.0 + H(sw.α, sw.u, t) ⋅ sw.qb)
Base.zero(sw::SmithWilson, t) = Continuous(sw.ufr - log(1.0 + H(sw.α, sw.u, t) ⋅ sw.qb) / t)
Base.zero(sw::SmithWilson, t, cf::CompoundingFrequency) = convert(cf, zero(sw, t))

function SmithWilson(times::AbstractVector, cashflows::AbstractMatrix, prices::AbstractVector; ufr, α)
    Q = Diagonal(exp.(-ufr * times)) * cashflows
    q = vec(sum(Q, dims = 1))  # We want q to be a column vector
    QHQ = Q' * H(α, times) * Q
    b = QHQ \ (prices - q)
    Qb = Q * b
    return SmithWilson(times, Qb; ufr = ufr, α = α)
end

""" 
    timepoints(zcq::Vector{ZeroCouponQuote})
    timepoints(bbq::Vector{BulletBondQuote})

Return the times associated with the `cashflows` of the instruments.
"""
function timepoints(qs::Vector{Q}) where {Q<:ObservableQuote}
    frequency = maximum(q.frequency for q in qs)
    timestep = 1 / frequency
    maturity = maximum(q.maturity for q in qs)
    return [timestep:timestep:maturity...]
end


"""
    cashflows(interests, maturities, frequency)
    timepoints(zcq::Vector{ZeroCouponQuote})
    timepoints(bbq::Vector{BulletBondQuote})

Produce a cash flow matrix for a set of instruments with given `interests` and `maturities`
and a given payment frequency `frequency`. All instruments are assumed to have their first payment at time 1/`frequency`
and have their last payment at the largest multiple of 1/`frequency` less than or equal to the input maturity.
"""
function cashflows(interests, maturities, frequencies)
    frequency = lcm(frequencies)
    fq = inv.(frequencies)
    timestep = 1 / frequency
    floored_mats = floor.(maturities ./ timestep) .* timestep
    times = timestep:timestep:maximum(floored_mats)
    # we need to determine the coupons in relation to the payment date, not time zero
    time_adj = floored_mats .% fq

    cashflows = [
        # if on a coupon date and less than maturity, pay coupon
        ((((t + time_adj[instrument]) % fq[instrument] ≈ 0) && t <= floored_mats[instrument]) ? interests[instrument] / frequencies[instrument] : 0.0) +
        (t ≈ floored_mats[instrument] ? 1.0 : 0.0) # add maturity payment
        for t in times, instrument = 1:length(interests)
    ]

    return cashflows
end

function cashflows(qs::Vector{Q}) where {Q<:ObservableQuote}
    yield = [q.yield for q in qs]
    maturity = [q.maturity for q in qs]
    frequency = [q.frequency for q in qs]
    return cashflows(yield, maturity, frequency)
end

# Utility methods for calibrating Smith-Wilson directly from quotes
function SmithWilson(zcq::Vector{ZeroCouponQuote}; ufr, α)
    n = length(zcq)
    maturities = [q.maturity for q in zcq]
    prices = [q.price for q in zcq]
    return SmithWilson(maturities, Matrix{Float64}(I, n, n), prices; ufr = ufr, α = α)
end

function SmithWilson(swq::Vector{SwapQuote}; ufr, α)
    times = timepoints(swq)
    cfs = cashflows(swq)
    ones(length(swq))
    return SmithWilson(times, cfs, ones(length(swq)), ufr = ufr, α = α)
end

function SmithWilson(bbq::Vector{BulletBondQuote}; ufr, α)
    times = timepoints(bbq)
    cfs = cashflows(bbq)
    prices = [q.price for q in bbq]
    return SmithWilson(times, cfs, prices, ufr = ufr, α = α)
end

# https://github.com/dpsanders/hands_on_julia/blob/master/during_sessions/Fractale%20de%20Newton.ipynb
newton(f, f′, x) = x - f(x) / f′(x)
function solve(g, g′, x0, max_iterations = 100)
    x = x0

    tolerance = 2 * eps(x0)
    iteration = 0

    while (abs(g(x) - 0) > tolerance && iteration < max_iterations)
        x = newton(g, g′, x)
        iteration += 1
    end

    return x
end

function bootstrap(rates, maturities, settlement_frequency; interp_function = linear_interp)
    discount_vec = zeros(length(rates)) # construct a placeholder discount vector matching maturities
    # we have to take the first rate as the starting point
    discount_vec[1] = discount(Constant(rates[1]), maturities[1])

    for t = 2:length(maturities)
        if isnothing(settlement_frequency[t])
            # no settlement before maturity
            discount_vec[t] = discount(Constant(rates[t]), maturities[t])
        else
            # need to account for the interim cashflows settled
            times = settlement_frequency[t]:settlement_frequency[t]:maturities[t]
            cfs = [rate(rates[t]) * settlement_frequency[t] for s in times]
            cfs[end] += 1

            function pv(v_guess)
                v = interp_function([[0.0]; maturities[1:t]], vcat(1.0, discount_vec[1:t-1], v_guess...))
                return sum(v.(times) .* cfs)
            end
            target_pv = sum(map(t2 -> discount(Constant(rates[t]), t2), times) .* cfs)
            root_func(v_guess) = pv(v_guess) - target_pv
            root_func′(v_guess) = ForwardDiff.derivative(root_func, v_guess)
            discount_vec[t] = solve(root_func, root_func′, rate(rates[t]))
        end

    end
    zero_vec = -log.(discount_vec) ./ maturities
    return linear_interp([0.0; maturities], [first(zero_vec); zero_vec])
end

## Generic and Fallbacks
"""
    discount(rate,to)
    discount(rate,from,to)

The discount factor for the `rate` for times `from` through `to`. If rate is a `Real` number, will assume a `Constant` interest rate.
"""
discount(yc, time) = exp(-yc.zero(time) * time)
discount(rate::Rate{<:Real,<:CompoundingFrequency}, from, to) = discount(Constant(rate), from, to)
discount(rate::Rate{<:Real,<:CompoundingFrequency}, to) = discount(Constant(rate), to)



discount(yc, from, to) = discount(yc, to) / discount(yc, from)

"""
    forward(curve,from,to,CompoundingFrequency=Periodic(1))

The forward `Rate` implied by the curve between times `from` and `to`.
"""
function forward(yc, from, to)
    return forward(yc, from, to, Periodic(1))
end

function forward(yc, from, to, cf::T) where {T<:CompoundingFrequency}

    r = Periodic((accumulation(yc, to) / accumulation(yc, from))^(1 / (to - from)) - 1, 1)
    return convert(cf, r)
end

function forward(yc, from)
    to = from + 1
    return forward(yc, from, to)
end

"""
    accumulation(rate,from,to)

The accumulation factor for the `rate` for times `from` through `to`. If rate is a `Real` number, will assume a `Constant` interest rate.
"""
function accumulation(y::T, time) where {T<:AbstractYield}
    return 1 ./ discount(y, time)
end
accumulation(rate::Rate{<:Real,<:CompoundingFrequency}, to) = accumulation(Constant(rate), to)

function accumulation(y::T, from, to) where {T<:AbstractYield}
    return 1 ./ discount(y, from, to)
end
accumulation(rate::Rate{<:Real,<:CompoundingFrequency}, from, to) = accumulation(Constant(rate), from, to)

## Curve Manipulations
struct RateCombination{T,U,V} <: AbstractYield
    r1::T
    r2::U
    op::V
end

rate(rc::RateCombination, time) = rc.op(rate(rc.r1, time), rate(rc.r2, time))
function discount(rc::RateCombination, time)
    a1 = discount(rc.r1, time)^(-1 / time) - 1
    a2 = discount(rc.r2, time)^(-1 / time) - 1
    return 1 / (1 + rc.op(a1, a2))^time
end

"""
    Yields.AbstractYield + Yields.AbstractYield

The addition of two yields will create a `RateCombination`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will be added together.
"""
function Base.:+(a::AbstractYield, b::AbstractYield)
    return RateCombination(a, b, +)
end

function Base.:+(a::Constant, b::Constant)
    a_kind = rate(a).compounding
    rate_new_basis = rate(convert(a_kind, rate(b)))
    return Constant(
        Rate(
            rate(a.rate) + rate_new_basis,
            a_kind
        )
    )
end

function Base.:+(a::T, b) where {T<:AbstractYield}
    return a + Constant(b)
end

function Base.:+(a, b::T) where {T<:AbstractYield}
    return Constant(a) + b
end

"""
    Yields.AbstractYield - Yields.AbstractYield

The subtraction of two yields will create a `RateCombination`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the second curves will be subtracted from the first.
"""
function Base.:-(a::AbstractYield, b::AbstractYield)
    return RateCombination(a, b, -)
end

function Base.:-(a::Constant, b::Constant)
    a_kind = rate(a).compounding
    rate_new_basis = rate(convert(a_kind, rate(b)))
    return Constant(
        Rate(
            rate(a.rate) - rate_new_basis,
            a_kind
        )
    )
end

function Base.:-(a::T, b) where {T<:AbstractYield}
    return a - Constant(b)
end

function Base.:-(a, b::T) where {T<:AbstractYield}
    return Constant(a) - b
end

linear_interp(xs, ys) = Interpolations.extrapolate(
    Interpolations.interpolate((xs,), ys, Interpolations.Gridded(Interpolations.Linear())),
    Interpolations.Line()
)

function Base.show(io::IO, curve::T) where {T<:AbstractYield}
    r = zero(curve, 1)
    ylabel = isa(r.compounding, Continuous) ? "Continuous" : "Periodic($(r.compounding.frequency))"
    display(lineplot(t -> rate(zero(curve, t)), 0.01, 5, xlabel = "time", ylabel = ylabel, compact = true, name = "Zero rates"))
end

end
