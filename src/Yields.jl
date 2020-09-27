module Yields

using Dierckx


export ZeroCurve,ConstantYield, ForwardYields,
rate,forward,TreasuryYieldCurve, disc, AbstractYieldCurve

"""
An AbstractInterestCurve is an object which can be called with:

- `rate` for the spot rate at a given time
- `disc` for the spot discount rate at a given time

"""
abstract type AbstractYieldCurve end

# make interest curve broadcastable so that you can broadcast over multiple`time`s in `interest_rate`
Base.Broadcast.broadcastable(ic::AbstractYieldCurve) = Ref(ic) 

struct YieldCurve <: AbstractYieldCurve
    rates
    maturities
    spline
end

struct ConstantYield <: AbstractYieldCurve
    rate
end

rate(c::ConstantYield,time) = c.rate
disc(c::ConstantYield,time) = 1/ (1 + c.rate) ^ time


function ZeroCurve(rates,maturities)
    return YieldCurve(rates,maturities,Spline1D(maturities,rates))
end

"""
    Forwards(rate_vector)

Takes a vector of 1-period forward rates and constructs a discount curve.
"""
function ForwardYields(rate_vector)
    zeros = similar(rate_vector)
    zeros[1] = rate_vector[1]
    for i in 2:length(rate_vector)
        zeros[i] = (prod(1 .+ rate_vector[1:i]))  ^ (1/i) - 1
    end
    return ZeroCurve(zeros,1:length(rate_vector))
end

function TreasuryYieldCurve(rates,maturities)
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

struct RateCombination <: AbstractYieldCurve
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

function Base.:+(a::AbstractYieldCurve,b::AbstractYieldCurve)
   return RateCombination(a, b,+) 
end


function Base.:-(a::AbstractYieldCurve,b::AbstractYieldCurve)
    return RateCombination(a, b,-) 
end

end
