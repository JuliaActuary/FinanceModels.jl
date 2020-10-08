module Yields

using Dierckx

# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCruve
export rate, discount, forward, Yield
# USTreasury,  AbstractYield
# Zero,Constant, Forward
"""
An AbstractYield is an object which can be called with:

- `rate(yield,time)` for the spot rate at a given time
- `discount(yield,time)` for the spot discount rate at a given time

"""
abstract type AbstractYield end

# make interest curve broadcastable so that you can broadcast over multiple`time`s in `interest_rate`
Base.Broadcast.broadcastable(ic::AbstractYield) = Ref(ic) 

struct YieldCurve <: AbstractYield
    rates # spot rates
    maturities
    spline
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
    rate
end

rate(c::Constant,time) = c.rate
discount(c::Constant,time) = 1 / (1 + rate(c, time))^time

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

function Zero(rates, maturities)
    # bump to a constant yield if only given one rate
    length(rates) == 1 && return Constant(rate[1])

    return YieldCurve(
        rates,
        maturities,
        Spline1D(
            maturities,
            rates; 
            k=min(3, length(rates) - 1) # spline dim has to be less than number of given rates
            )
        )
end


function Zero(rates)
    # bump to a constant yield if only given one rate
    maturities = collect(1:length(rates))
    return Zero(rates, maturities)
end

"""
Construct a curve given a set of bond yields priced at par with a single coupon per period.
"""
function Par(rate, maturity;)
    # bump to a constant yield if only given one rate
    if length(rate) == 1
        return Constant(rate[1])
    end

    spot = similar(rate) 

    spot[1] = rate[1]

    for i in 2:length(rate)
        coupon_pv = sum(rate[i] / (1 + spot[j])^maturity[j] for j in 1:i - 1) # not including the one paid at maturity

        spot[i] = ((1 + rate[i]) / (1 - coupon_pv))^(1 / maturity[i]) - 1
    end



    return YieldCurve(
        rate,
        maturity,
        Spline1D(
            maturity,
            spot; 
            k=min(3, length(rate) - 1) # spline dim has to be less than number of given rates
            )
        )
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

function USTreasury(rates, maturities)
    z = zeros(length(rates))

    # use the discount rate for T-Bills with maturities <= 1 year
    for (i, (rate, mat)) in enumerate(zip(rates, maturities))
        
        if mat <= 1 
            z[i] = rate
        else
            # uses spline b/c of common, but uneven maturities often present under 1 year.
            curve = Spline1D(maturities, z)
            pmts = [rate / 2 for t in 0.5:0.5:mat] # coupons only
            pmts[end] += 1 # plus principal

            discount =  1 ./ (1 .+ curve.(0.5:0.5:(mat - .5)))
            z[i] = ((1 - sum(discount .* pmts[1:end - 1])) / pmts[end])^- (1 / mat) - 1

        end




        
    end

    return YieldCurve(rates, maturities, Spline1D(maturities, z))


    return YieldCurve(
        rate,
        maturity,
        Spline1D(
            maturity,
            spot; 
            k=min(3, length(rate) - 1) # spline dim has to be less than number of given rates
            )
        )
end

function ParYieldCurve(rates, maturities)

end


## Generic and Fallbacks
"""
    rate(yield,time)

The spot rate at `time` for the given `yield`.
"""
rate(yc,time) = yc.spline(time)

"""
    discount(yield,time)

The discount factor for the `yield` from time zero through `time`.
"""
discount(yc,time) = 1 / (1 + rate(yc, time))^time

"""
    discount(yield,from,to)

The discount factor for the `yield` from time `from` through `to`.
"""
discount(yc,from,to) = discount(yc, to) / discount(yc, from)

function forward(yc, from, to)
    return (accumulate(yc, to) / accumulate(yc, from))^(1 / (to - from)) - 1
end
function forward(yc, to)
    from = to - 1 
    return forward(yc, from, to)
end

"""
    accumulate(yield,time)

The accumulation factor for the `yield` from time zero through `time`.
"""
function Base.accumulate(y::T, time) where {T <: AbstractYield}
    return 1 / discount(y, time)
end

function Base.accumulate(y::T,from,to) where {T <: AbstractYield}
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
    yield(rate)
    yield(forwards)

Yields provides a default, convienience construction for an AbstractYield.

"""

function Yield(i::T) where {T<:Real}
    return Constant(i)
end

function Yield(i::Vector{T}) where {T<:Real}
    return Forward(i)
end

end
