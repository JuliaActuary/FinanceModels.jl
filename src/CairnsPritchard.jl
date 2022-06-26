struct CairnsPritchardFixedC{T} <: YieldCurveFitParameters
    c::Vector{T}
end

struct CairnsPritchardFitC <: YieldCurveFitParameters
    n::Int
end


struct CairnsPritchardCurve
    c
    b
    # function CairnsPritchardCurve(c,b)

end

rates =[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
mats = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]
m(x, params) = @. params[1] + params[2]*exp(-0.2*x) + params[3] * exp(-0.4*x) + params[4] * exp(-0.8*x)
Curvefit(rates,mats,m,ones(4),LBFGS())