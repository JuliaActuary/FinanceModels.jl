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
    A = DataInterpolations.Curvefit(
        yields,
        maturities,
        m,
        initial_params,
        DataInterpolations.LBFGS(),
        true,
        __cairns_lb(fit),
        __cairns_ub(fit),
        )
    b = A.pmin
    return CairnsPritchardCurve(b,fit.c)
end

function __fit(fit::CairnsPritchardFitC,func,yields,maturities)
    n = fit.n
    # params = [b1, b2, ..., c1, c2, ...]
    function m(x,params) 
        c = CairnsPritchardCurve(first(params,n+1),last(params,n))
        return rate.(func.(c,x))
    end

    initial_params = ones(2*n+1) 
    A = DataInterpolations.Curvefit(
        yields,
        maturities,
        m,
        initial_params,
        DataInterpolations.LBFGS(),
        true,
        __cairns_lb(fit),
        __cairns_ub(fit),
        )
    p = A.pmin
    return CairnsPritchardCurve(first(p,n+1),last(p,n))
end

__cairns_lb(fit::CairnsPritchardFixedC) = fill(-10.,length(fit.c)+1)
__cairns_lb(fit::CairnsPritchardFitC) = fill(-10.,fit.n*2+1)
__cairns_ub(fit::CairnsPritchardFixedC) = fill(10.,length(fit.c)+1)
__cairns_ub(fit::CairnsPritchardFitC) = fill(10.,fit.n*2+1)

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