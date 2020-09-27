module Yields

using Dierckx

export ZeroCurve, rate,forward,TreasuryYieldCurve

"""
An AbstractInterestCurve is an object which can be called with:

- `rate` for the spot rate at a given time
- `disc` for the spot discount rate at a given time

"""
abstract type AbstractInterestCurve end

struct InterestCurve <: AbstractInterestCurve
    rates
    maturities
    spline
end

struct ConstantCurve <: AbstractInterestCurve
    rate
end

rate(c::ConstantCurve,time) = c.rate
disc(c::ConstantCurve,time) = 1/ (1 + c.rate) ^ 2


function ZeroCurve(rates,maturities)
    return YieldCurve(rates,maturities,Spline1D(maturities,rates))
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

end
