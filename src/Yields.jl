module Yields

using Dierckx

# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCruve
export rate, disc, forward
# USTreasury,  AbstractYield
# Zero,Constant, Forward
"""
An AbstractInterestCurve is an object which can be called with:

- `rate` for the spot rate at a given time
- `disc` for the spot discount rate at a given time

"""
abstract type AbstractYield end

# make interest curve broadcastable so that you can broadcast over multiple`time`s in `interest_rate`
Base.Broadcast.broadcastable(ic::AbstractYield) = Ref(ic) 

struct YieldCurve <: AbstractYield
    rates
    maturities
    spline
end

struct Constant <: AbstractYield
    rate
end

rate(c::Constant,time) = c.rate
disc(c::Constant,time) = 1/ (1 + rate(c,time)) ^ time


function Zero(rates,maturities)
    # bump to a constant yield if only given one rate
    length(rates) == 1 && return Constant(rate[1])

    return YieldCurve(
        rates,
        maturities,
        Spline1D(
            maturities,
            rates; 
            k=min(3,length(rates)-1) # spline dim has to be less than number of given rates
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
        zeros[i] = (prod(1 .+ rate_vector[1:i]))  ^ (1/i) - 1
    end
    return Zero(zeros,1:length(rate_vector))
end

function USTreasury(rates,maturities)
    z = zeros(length(rates))

    # use the discount rate for T-Bills with maturities <= 1 year
    for (i,(rate,mat)) in enumerate(zip(rates,maturities))
        
        if mat <= 1 
            z[i] = rate
        else
            i, rate, mat
            curve = Spline1D(maturities,z)
            pmts = [rate / 2 for t in 0.5:0.5:mat] # coupons only
            pmts[end] += 1 # plus principal

            disc =  1 ./ (1 .+ curve.(0.5:0.5:(mat - .5)))
            z[i] = ((1 - sum(disc .* pmts[1:end-1])) / pmts[end]) ^ - (1/mat) - 1

        end




        
    end

    return YieldCurve(rates,maturities,Spline1D(maturities,z))


end

function ParYieldCurve(rates,maturities)

end

rate(yc,time) = yc.spline(time)

disc(yc,time) = 1 / (1 + rate(yc,time)) ^ time

function forward(yc,from,to) 
    (rate(yc,to) * to - rate(yc,from) * from) / (to - from)
end

struct RateCombination <: AbstractYield
    r1
    r2
    op
end

rate(rc::RateCombination,time) = rc.op(rate(rc.r1,time) ,rate(rc.r2,time))
function disc(rc::RateCombination,time) 
    r = rc.op(rate(rc.r1,time) ,rate(rc.r2,time))
    return 1 / (1+r) ^ time
end

### Curve Manipulations

function Base.:+(a::AbstractYield,b::AbstractYield)
   return RateCombination(a, b,+) 
end


function Base.:-(a::AbstractYield,b::AbstractYield)
    return RateCombination(a, b,-) 
end

end
