abstract type CairnsPritchardFit <: YieldCurveFitParameters end
struct CairnsPritchardFixedC{T} <: CairnsPritchardFit
    c::Vector{T}
end

struct CairnsPritchardFitC <: CairnsPritchardFit
    n::Int
end


struct CairnsPritchardCurve <: ParametricModel
    b
    c
end
function __fit(fit::CairnsPritchardFixedC,func,yields,maturities)
    function m(x,params) 
        c = CairnsPritchardCurve(params,fit.c)
        return rate.(func.(c,x))
    end

    initial_params = ones(length(fit.c)+1) 
    A = DataInterpolations.Curvefit(yields,maturities,m,initial_params,DataInterpolations.LBFGS())
    b = A.pmin
    return CairnsPritchardCurve(b,fit.c)
end

function Zero(fit::T,yields, maturities=eachindex(yields)) where {T<:CairnsPritchardFit}
    __fit(fit,zero,yields,maturities)
end
function Par(fit::T,yields, maturities=eachindex(yields)) where {T<:CairnsPritchardFit}
    __fit(fit,par,yields,maturities)
end
function Forward(fit::T,yields, maturities=eachindex(yields)) where {T<:CairnsPritchardFit}
    __fit(fit,forward,yields,maturities)
end

function Base.zero(cpc::CairnsPritchardCurve, t)
    b = cpc.b
    c = cpc.c
    Continuous.(sum(b[i+1]*exp.(-cᵢ.*t) for (i,cᵢ) in enumerate(c)) .+ b[1])
end

function discount(cpc::CairnsPritchardCurve, t)
    discount.(zero.(cpc,t),t)
end





    
# rates =[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
# mats = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]
# m(x, params) = @. params[1] + params[2]*exp(-0.2*x) + params[3] * exp(-0.4*x) + params[4] * exp(-0.8*x)
# Curvefit(rates,mats,m,ones(4),LBFGS())