module Yields

import Interpolations
import ForwardDiff

# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCruve
export rate, discount, accumulation,forward, Yield, Rate
# USTreasury,  AbstractYield
# Zero,Constant, Forward
"""
An AbstractYield is an object which can be called with:

- `rate(yield,time)` for the spot rate at a given time
- `discount(yield,time)` for the spot discount rate at a given time

"""
abstract type AbstractYield end

abstract type CompoundingFrequency end

struct Continuous <: CompoundingFrequency end

struct Periodic <: CompoundingFrequency
    frequency::Int
end

struct Rate
    compounding
    value
end

Rate(x) = Rate(Periodic(1),x)
Base.convert(T::CompoundingFrequency,r::Rate) = Base.convert(r.compounding,T,r)
function Base.convert(from::Continuous,to::Continuous,r)
    return r.value
end
function Base.convert(from::Continuous,to::Periodic,r)
    to.frequency * (exp(r.value/to.frequency) - 1)
end
function Base.convert(from::Periodic,to::Continuous,r)
    from.frequency * log(1 + r.value / from.frequency)
end

# make interest curve broadcastable so that you can broadcast over multiple`time`s in `interest_rate`
Base.Broadcast.broadcastable(ic::AbstractYield) = Ref(ic) 

struct YieldCurve <: AbstractYield
    rates
    maturities
    discount
end

# Wrapping a a scalar value in this type allows for dispatch to operate as intended 
# (otherwise `Base.accumulate(<:Real,<:Real) tries to do something other than accumulate interest)
"""
    Constant(rate)

Construct a yield object where the spot rate is constant for all maturities.

# Examples

```julia-repl
julia> y = Yields.Constant(0.05)
julia> discount(y,2)
0.9070294784580498     # 1 / (1.05) ^ 2
```
"""
struct Constant <: AbstractYield
    compounding::CompoundingFrequency
    rate
end

function Constant(rate::T) where {T <: Real}
    return Constant(Periodic(1),rate)
end

rate(c::Constant) = c.rate
rate(c::Constant,time) = c.rate
discount(c::Constant,time) = 1 / (1 + rate(c, time))^time
discount(c::T,time) where {T <: Real} = discount(Constant(c),time)

discount(r::Constant,time) = 1 / accumulation(r,time)


accumulation(r::Constant,time) = accumulation(r.compounding,r,time)
accumulation(::Continuous,r::Constant,time) = exp(r.rate * time)
accumulation(::Periodic,r::Constant,time) = (1 + r.rate / r.compounding.frequency) ^ (r.compounding.frequency * time)


"""
    Step(rates,times)

Create a yield curve object where the applicable rate is the effective rate of interest applicable until corresponding time.

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
    compounding
    rates
    times
end

Step(rates) = Step(rates, collect(1:length(rates)))

Step(rates,times) = Step(Periodic(1),rates,times)
Step(cf::CompoundingFrequency, rates) = Step(cf,rates,collect(1:length(rates)))

function rate(y::Step, time)
    i = findfirst(t -> time <= t, y.times)
    if isnothing(i)
        return y.rates[end]
    else
        return y.rates[i]
    end
end
discount(y::Step,time) = 1 / accumulation(y::Step,time)

accumulation(y::Step,time) = accumulation(y.compounding,y,time)

function accumulation(::Periodic,y::Step, time)
    m = y.compounding.frequency
    v =  (1 + y.rates[1] / m )^ (m * min(y.times[1], time))

    if y.times[1] >= time
        return v
    end

    for i in 2:length(y.times)

        if y.times[i] >= time
            # take partial accumulation and break
            v *= (1 + y.rates[i]/m)^(m*(time - y.times[i - 1]))
            break
        else
            # take full accumulation and continue
            v *=  (1 + y.rates[i]/m)^(m*(y.times[i] - y.times[i - 1]))
        end

    end

    return v
end

function accumulation(::Continuous,y::Step, time)
    v = exp(y.rates[1] *min(y.times[1], time))

    if y.times[1] >= time
        return v
    end

    for i in 2:length(y.times)

        if y.times[i] >= time
            # take partial accumulation and break
            v *= exp(y.rates[i]*(time - y.times[i - 1]))
            break
        else
            # take full accumulation and continue
            v *=  exp(y.rates[i]*(y.times[i] - y.times[i - 1]))
        end

    end

    return v
end

function Zero(p::CompoundingFrequency,rates, maturities)
    # bump to a constant yield if only given one rate
    length(rates) == 1 && return Constant(p,rate[1])

    discounts = map(zip(rates,maturities)) do (r,m)
        discount(Constant(p,r),m)
    end

    # discount at time 0 should always be 1.0
    if ~iszero(first(maturities))
        maturities = [zero(eltype(maturities));maturities]
        discounts = [1.; discounts]
    end


    return YieldCurve(
        p,
        rates,
        maturities,
        linear_interp(maturities,discounts)
    )
end

Zero(rates,maturities) = Zero(Periodic(1),rates,maturities)

Zero(rates) = Zero(Periodic(1),rates)
function Zero(p::CompoundingFrequency,rates)
    maturities = collect(1:length(rates))
    return Zero(p,rates, maturities)
end

"""
    Par(rates,maturities)
    Par(rates)
    Par(::CompoundingFrequency,...)

Construct a curve given a set of bond yields priced at par. If no `CompoundingFrequency` pased, will assume once per period.
"""
function Par(p::Periodic,rates, maturities)
    # bump to a constant yield if only given one rate
    if length(rates) == 1
        return Constant(p,rates[1])
    end

    m = p.frequency

    discounts = similar(rates) 

    discounts[1] = discount(Constant(p,rates[1]),maturities[1])
    cashflows =  similar(rates,length(rates) * m)
    times =  collect(1/m:1/m:m*last(maturities))

    for i in 2:length(rates)
        # solve for the rates that gets you to the same yield for the cashflows
        for t in 1:i*m
            cashflows[t] = rates[i] / m
        end
        cashflows[i*m] += 1. # add the final par payment

        ytm_disc = map(t->discount(Constant(p,rates[i]),t),times[1:m*i])
        target_pv = @views sum(cashflows[1:m*i] .*  ytm_disc)
        
        function f(x) 
            
            df = linear_interp(times[1:i],[discounts[1:i-1];x])
            return sum(cashflows[t] * df(t) for t in 1:m*i) - target_pv
        end
        f′(x) = ForwardDiff.derivative(f,x)


        discounts[i] = newton(f,f′,rates[i])

    end

    # discount at time 0 should always be 1.0
    if ~iszero(first(maturities))
        maturities = [zero(eltype(maturities));maturities]
        discounts = [1.; discounts]
    end

    return YieldCurve(
        p,
        rates,
        maturities,
        linear_interp(maturities,discounts)
        )
end

function Par(p::Vector{Periodic},rates, maturities)
    # bump to a constant yield if only given one rate
    if length(rates) == 1
        return Constant(p,rates[1])
    end

    m = p.frequency

    discounts = similar(rates) 

    discounts[1] = discount(Constant(p,rates[1]),maturities[1])
    cashflows =  similar(rates,length(rates) * m)
    times =  collect(1/m:1/m:m*last(maturities))

    for i in 2:length(rates)
        # solve for the rates that gets you to the same yield for the cashflows
        for t in 1:i*m
            cashflows[t] = rates[i] / m
        end
        cashflows[i*m] += 1. # add the final par payment

        ytm_disc = map(t->discount(Constant(p,rates[i]),t),times[1:m*i])
        target_pv = @views sum(cashflows[1:m*i] .*  ytm_disc)
        
        function f(x) 
            
            df = linear_interp(times[1:i],[discounts[1:i-1];x])
            return sum(cashflows[t] * df(t) for t in 1:m*i) - target_pv
        end
        f′(x) = ForwardDiff.derivative(f,x)


        discounts[i] = newton(f,f′,rates[i])

    end

    # discount at time 0 should always be 1.0
    if ~iszero(first(maturities))
        maturities = [zero(eltype(maturities));maturities]
        discounts = [1.; discounts]
    end

    return YieldCurve(
        p,
        rates,
        maturities,
        linear_interp(maturities,discounts)
        )
end

# https://github.com/dpsanders/hands_on_julia/blob/master/during_sessions/Fractale%20de%20Newton.ipynb
newton(f, f′, x) = x - f(x) / f′(x)
function solve(g, g′, x0, max_iterations=100)  # valeur par defaut
    x = x0

    tolerance = 2*eps(x0)
    iteration = 0

    while (abs(g(x) - 0) > tolerance && iteration < max_iterations)
        x = newton(g, g′, x)        
        iteration += 1
    end
    
    return x
end

"""
    Forward(rate_vector)

Takes a vector of 1-period forward rates and constructs a discount curve.
"""
function Forward(rate_vector)
    zeros = similar(rate_vector)
    zeros[1] = rate_vector[1]
    for i in 2:length(rate_vector)
        zeros[i] = (prod(1 .+ rate_vector[1:i]))^(1 / i) - 1
    end
    return Zero(zeros, 1:length(rate_vector))
end

function Forward(rate_vector, times)
    disc_v = similar(rate_vector)
    disc_v[1] = 1 / (1 + rate_vector[1])^times[1]
    for i in 2:length(rate_vector)
        ∇t = times[i] - times[i - 1]
        disc_v[i] = disc_v[i - 1] / (1 + rate_vector[i])^∇t
    end

    return Zero(1 ./ disc_v.^(1 ./ times) .- 1, times)
end

"""
    USTreasury(rates,maturities)

Takes CMT yields (bond equivalent), and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual.
"""
function USTreasury(rates, maturities)
    z = zeros(length(rates))

    # use the discount rate for T-Bills with maturities <= 1 year
    for (i, (rate, mat)) in enumerate(zip(rates, maturities))
        
        if mat <= 1 
            z[i] = (1 + rate * mat) ^ (1/mat) -1
        else
            # uses interpolation b/c of common, but uneven maturities often present under 1 year.
            curve = linear_interp(maturities, z)
            pmts = [rate / 2 for t in 0.5:0.5:mat] # coupons only
            pmts[end] += 1 # plus principal

            discount =  1 ./ (1 .+ curve.(0.5:0.5:(mat - .5)))
            z[i] = ((1 - sum(discount .* pmts[1:end - 1])) / pmts[end])^- (1 / mat) - 1

        end




        
    end

    return YieldCurve(rates, maturities, linear_interp(maturities, z))


    return YieldCurve(
        rate,
        maturity,
        linear_interp(maturity,rate)
        )
end


## Generic and Fallbacks
"""
    rate(yield,time)

The annual effective spot rate at `time` for the given `yield`.
"""
rate(yc,time) = rate(yc.compounding,yc,time)
rate(::Continuous,yc,time) = log(accumulation(yc,time)) / time - 1
function rate(p::Periodic,yc,time) 
    m = p.frequency
    return (accumulation(yc,time) ^ (1 / (m * time)) - 1) * m
end

"""
    discount(yield,time)

The discount factor for the `yield` from time zero through `time`. If yield is a `Real` number, will assume a `Constant` interest rate.
"""
discount(yc,time) = yc.discount(time) 

"""
    discount(yield,from,to)

The discount factor for the `yield` from time `from` through `to`.
"""
discount(yc,from,to) = discount(yc, to) / discount(yc, from)

forward(yc,from,to) = forward(yc.compounding,yc,from,to)
function forward(p::Periodic,yc, from, to)
    x = (accumulation(yc, to) / accumulation(yc, from))
    m = p.frequency
    # convert to periodic before returning
    return m*(x ^(1 / (to - from) * m) - 1)
end

function forward(p::Continuous,yc, from, to)
    return log(accumulation(yc, to) / accumulation(yc, from))
end

function forward(yc, to)
    from = to - 1 
    return forward(yc, from, to)
end

"""
    accumulation(yield,time)

The accumulation factor for the `yield` from time zero through `time`.
"""
function accumulation(y::T, time) where {T <: AbstractYield}
    return 1 / discount(y, time)
end

function accumulation(y::T,from,to) where {T <: AbstractYield}
    return 1 / discount(y,from,to)
end

## Curve Manipulations
struct RateCombination <: AbstractYield
    r1
    r2
    op
end

rate(rc::RateCombination,time) = rc.op(rate(rc.r1, time), rate(rc.r2, time))
function discount(rc::RateCombination, time) 
    r = rc.op(rate(rc.r1, time), rate(rc.r2, time))
    return 1 / (1 + r)^time
end

"""
    Yields.AbstractYield + Yields.AbstractYield

The addition of two yields will create a `RateCombination`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will be added together.
"""
function Base.:+(a::AbstractYield, b::AbstractYield)
    return RateCombination(a, b, +) 
end

function Base.:+(a::Constant, b::Constant) where {T<:AbstractYield}
    return Constant(a.rate + b.rate)
end

function Base.:+(a::T, b) where {T<:AbstractYield}
    return a + Yield(b)
end

function Base.:+(a, b::T) where {T<:AbstractYield}
    return Yield(a) + b
end

"""
    Yields.AbstractYield - Yields.AbstractYield

The subtraction of two yields will create a `RateCombination`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the second curves will be subtracted from the first.
"""
function Base.:-(a::AbstractYield, b::AbstractYield)
    return RateCombination(a, b, -) 
end

function Base.:-(a::Constant, b::Constant)
    return Constant(a.rate - b.rate) 
end

function Base.:-(a::T, b) where {T<:AbstractYield}
    return a - Yield(b)
end

function Base.:-(a, b::T) where {T<:AbstractYield}
    return Yield(a) - b
end

""" 
    Yield(rate)
    Yield(forwards)

Yields provides a default, convienience construction for an AbstractYield.

"""
function Yield(i::T) where {T<:Real}
    return Constant(i)
end

function Yield(i::Vector{T}) where {T<:Real}
    return Forward(i)
end

linear_interp(xs,ys) = Interpolations.extrapolate(
    Interpolations.interpolate((xs,), ys, Interpolations.Gridded(Interpolations.Linear())), 
    Interpolations.Flat()
    ) 
end
