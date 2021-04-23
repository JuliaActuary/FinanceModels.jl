# TODO: accumulalate -> accumulation
module Yields

import Interpolations
import ForwardDiff

# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCruve
export rate, discount, accumulation,forward, Yield, Rate, Continuous, Periodic, rate, spot

abstract type CompoundingFrequency end
Base.Broadcast.broadcastable(x::T) where{T<:CompoundingFrequency} = Ref(x) 

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
    Continuous(rate)

A convinience constructor for Rate(x,Continuous())

```julia-repl
julia> Continuous(0.01)
Rate(0.01, Continuous())
```

See also: [`Periodic`](@ref)
"""
Continuous(rate) = Rate(rate,Continuous())

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

A convinience constructor for Rate(rate,Periodic(frequency)).

# Examples

Creating a semi-annual bond equivalent yield:

```julia-repl
julia> Periodic(0.01,2)
Rate(0.01, Periodic(2))
```

See also: [`Continuous`](@ref)
"""
Periodic(x,freq) = Rate(x,Periodic(freq))

struct Rate
    value
    compounding::CompoundingFrequency
end

# Base.:==(r1::Rate,r2::Rate) = (r1.value == r2.value) && (r1.compounding == r2.compounding)

"""
    Rate(rate[,frequency=1])
    Rate(rate,frequency::CompoundingFrequency)

Rate is a type that encapsulates an interest `rate` along wtih its compounding `frequency`.

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
Rate(rate) = Rate(rate,Periodic(1))
Rate(x,frequency::T) where {T<:Real} = isinf(frequency) ? Rate(x,Continuous()) : Rate(x,Periodic(frequency))

"""
    convert(r::Rate,T::CompoundingFrequency)

Returns a `Rate` with an equivalent discount but represented with a different compounding frequency.

# Examples

```
julia> r = Rate(0.01,Periodic(12))
Rate(0.01, Periodic(12))

julia> convert(r,Periodic(1))
Rate(0.010045960887181016, Periodic(1))

julia> convert(r,Continuous())
Rate(0.009995835646701251, Continuous())
```
"""
Base.convert(r::Rate,T::CompoundingFrequency) = convert(r,r.compounding,T)
function Base.convert(r,from::Continuous,to::Continuous)
    return r
end

function Base.convert(r,from::Continuous,to::Periodic)
    return Rate(to.frequency * (exp(r.value/to.frequency) - 1),to)
end

function Base.convert(r,from::Periodic,to::Continuous)
    return Rate(from.frequency * log(1 + r.value / from.frequency),to)
end

function Base.convert(r,from::Periodic,to::Periodic)
    c = convert(r,from,Continuous())
    return convert(c,Continuous(),to)
end

rate(r::Rate) = r.value


"""
An AbstractYield is an object which can be used as an argument to:

- zero-coupon spot rates viea [`zero`](@ref)
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

struct YieldCurve <: AbstractYield
    rates
    maturities
    discount # discount function for time
end

"""
    zero(curve,time)
    zero(curve,time,CompoundingFrequency)

Return the zero rate for the curve at the given time. If not specified, will use `Periodic(1)` compounding.
"""
Base.zero(c::YieldCurve,time) = zero(c,time,Periodic(1))
function Base.zero(c::YieldCurve,time,cf::Periodic)
    d = discount(c,time)
    i = Rate(cf.frequency*(d^(-1/(time*cf.frequency))-1),cf)
    return i
end

function Base.zero(c::YieldCurve,time,cf::Continuous)
    d = discount(c,time)
    i = log(1/d)/time
    return Rate(i,cf)
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
struct Constant <: AbstractYield
    rate
end

function Constant(rate::T) where {T <: Real}
    return Constant(Rate(rate,Periodic(1)))
end

rate(c::Constant) = c.rate
rate(c::Constant,time) = c.rate
discount(c::T,time) where {T <: Real} = discount(Constant(c),time)
discount(r::Constant,time) = 1 / accumulation(r,time)

accumulation(r::Constant,time) = accumulation(r.rate.compounding,r,time)
accumulation(::Continuous,r::Constant,time) = exp(rate(r.rate) * time)
accumulation(::Periodic,r::Constant,time) = (1 + rate(r.rate) / r.rate.compounding.frequency) ^ (r.rate.compounding.frequency * time)

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

    for i in 2:length(y.times)

        if y.times[i] >= time
            # take partial discount and break
            v /= (1 + y.rates[i])^(time - y.times[i - 1])
            break
        else
            # take full discount and continue
            v /=  (1 + y.rates[i])^(y.times[i] - y.times[i - 1])
        end

    end

    return v
end


"""
    Zero(rates,maturities)

Construct a yield curve with given zero-coupon spot `rates` at the given `maturities`. If `rates` is not a `Vector{Rate}`, will assume `Periodic(1)` type.
"""
function Zero(rates, maturities)
    # bump to a constant yield if only given one rate
    length(rates) == 1 && return Constant(rate[1])
    return YieldCurve(
        rates,
        maturities,
        linear_interp([0.;maturities],[1.;discount.(Constant.(rates),maturities)])
    )
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
function Par(rates::Vector{Rate}, maturities)
    # bump to a constant yield if only given one rate
    if length(rates) == 1
        return Constant(rate[1])
    end
    return YieldCurve(
            rates,
            maturities,
            # assume that maturities less than or equal to 12 months are settled once, otherwise semi-annual
            # per Hull 4.7
            bootstrap(rates,maturities,[m <= 1 ? nothing : 1 / r.compounding.frequency for (r,m) in zip(rates,maturities)])
        )
end

function Par(rates::Vector{T},maturities) where {T <: Real}
    return Par(Rate.(rates),maturities)  
end


"""
    Forward(rate_vector,maturities)

Takes a vector of 1-period forward rates and constructs a discount curve.
"""
function Forward(rates, maturities)
    # convert to zeros and pass to Zero
    disc_v = similar(rates)
    v = 1.

    for i in 1:length(rates)
        Δt = maturities[i] - (i == 1 ? 0 : maturities[i-1])
        v *= discount(Constant(rates[i]),Δt)
        disc_v[i] = v
    end

    z = (1. ./ disc_v) .^ ( 1 ./ maturities) .- 1 # convert disc_v to zero
    return Zero(z,maturities)
end

Forward(rates) = Forward(rates,collect(1:length(rates)))

"""
Takes CMT yields (bond equivalent), and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual.
"""
function CMT(rates::Vector{T}, maturities) where {T<:Real}
    rs = map(zip(rates,maturities)) do (r,m)
        if m <= 1
            Rate(r,Periodic(1 / m))
        else
            Rate(r,Periodic(2))
        end
    end

    CMT(rs,maturities)
end

function CMT(rates::Vector{Rate}, maturities)
    return YieldCurve(
            rates,
            maturities,
            # assume that maturities less than or equal to 12 months are settled once, otherwise semi-annual
            # per Hull 4.7
            bootstrap(rates,maturities,[m <= 1 ? nothing : 0.5 for m in maturities])
        )
end
    

"""
    OIS(rates,maturities)
Takes Overnight Index Swap rates, and assumes that instruments <= one year maturity are settled once and other agreements are settled quarterly with a corresponding CompoundingFrequency

"""
function OIS(rates::Vector{T}, maturities) where {T<:Real}
    rs = map(zip(rates,maturities)) do (r,m)
        if m <= 1
            Rate(r,Periodic(1 / m))
        else
            Rate(r,Periodic(4))
        end
    end

    return OIS(rs,maturities)
end
function OIS(rates::Vector{Rate}, maturities)
    return YieldCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise quarterly
        # per Hull 4.7
        bootstrap(rates,maturities,[m <= 1 ? nothing : 1/4 for m in maturities])
    )
end


# https://github.com/dpsanders/hands_on_julia/blob/master/during_sessions/Fractale%20de%20Newton.ipynb
newton(f, f′, x) = x - f(x) / f′(x)
function solve(g, g′, x0, max_iterations=100)
    x = x0

    tolerance = 2*eps(x0)
    iteration = 0

    while (abs(g(x) - 0) > tolerance && iteration < max_iterations)
        x = newton(g, g′, x)        
        iteration += 1
    end

    return x
end

function bootstrap(rates,maturities,settlement_frequency;interp_function=linear_interp)
    settlement_frequency,maturities,rates
    discount_vec = zeros(length(rates)) # construct a placeholder discount vector matching maturities
    # we have to take the first rate as the starting point
    discount_vec[1] = discount(Constant(rates[1]),maturities[1])

    for t in 2:length(maturities)
        if isnothing(settlement_frequency[t]) 
            # no settlment before maturity
            discount_vec[t] = discount(Constant(rates[t]),maturities[t])
        else
            # need to account for the interim cashflows settled
            times = settlement_frequency[t]:settlement_frequency[t]:maturities[t]
            cfs = [rate(rates[t]) * settlement_frequency[t] for s in times]
            cfs[end] += 1
            
            function pv(v_guess)
                v = interp_function([[0.];maturities[1:t]],vcat(1.,discount_vec[1:t-1],v_guess...))
                return sum(v.(times) .* cfs)
            end
            target_pv = sum(map(t2->discount(Constant(rates[t]),t2),times) .* cfs)
            root_func(v_guess) = pv(v_guess) - target_pv
            root_func′(v_guess) = ForwardDiff.derivative(root_func,v_guess)
            discount_vec[t] = solve(root_func,root_func′,rate(rates[t]))
        end

    end
    return linear_interp([[0.];maturities],[[1.];discount_vec])
end

## Generic and Fallbacks
"""
    discount(rate,to)
    discount(rate,from,to)

The discount factor for the `rate` for times `from` through `to`. If rate is a `Real` number, will assume a `Constant` interest rate.
"""
discount(yc,time) = yc.discount(time)
discount(rate::Rate,from,to) = discount(Constant(rate),from,to)
discount(rate::Rate,to) = discount(Constant(rate),to)



discount(yc,from,to) = discount(yc, to) / discount(yc, from)

function forward(yc, from, to)
    return (accumulation(yc, to) / accumulation(yc, from))^(1 / (to - from)) - 1
end
function forward(yc, from)
    to = from + 1
    return forward(yc, from, to)
end

"""
    accumulation(rate,from,to)

The accumulation factor for the `rate` for times `from` through `to`. If rate is a `Real` number, will assume a `Constant` interest rate.
"""
function accumulation(y::T, time) where {T <: AbstractYield}
    return 1 / discount(y, time)
end
accumulation(rate::Rate,to) = accumulation(Constant(rate),to)

function accumulation(y::T,from,to) where {T <: AbstractYield}
    return 1 / discount(y,from,to)
end
accumulation(rate::Rate,from,to) = accumulation(Constant(rate),from,to)

## Curve Manipulations
struct RateCombination <: AbstractYield
    r1
    r2
    op
end

rate(rc::RateCombination,time) = rc.op(rate(rc.r1, time), rate(rc.r2, time))
function discount(rc::RateCombination, time) 
    a1 = discount(rc.r1,time)^(-1/time) - 1  
    a2 = discount(rc.r2,time)^(-1/time) - 1
    return 1 / (1 + rc.op(a1,a2)) ^ time
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
    rate_new_basis = rate(convert(rate(b),a_kind))
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
    rate_new_basis = rate(convert(rate(b),a_kind))
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

linear_interp(xs,ys) = Interpolations.extrapolate(
    Interpolations.interpolate((xs,), ys, Interpolations.Gridded(Interpolations.Linear())), 
    Interpolations.Line()
    ) 
end
