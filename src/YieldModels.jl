module YieldModels

using Dierckx

export ZeroCurve, rate,forward,TreasuryYieldCurve


struct YieldCurve
    rates
    maturities
    spline
end

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
