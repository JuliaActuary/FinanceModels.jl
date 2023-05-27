abstract type AbstractModel end

# a model for when you don't really need a model
# (e.g. determining nominal cashflows for fixed income contract)
struct NullModel <: AbstractModel end

module Yield
    import ..AbstractModel
    import ..FinanceCore
    using FinanceCore: Continuous, discount

    export discount, rate
    
    abstract type AbstractYieldModel <: AbstractModel end

    struct Constant{R} <: AbstractYieldModel
        rate::R
    end

    function Constant(rate::R) where {R<:Real}
        Constant(FinanceCore.Rate(rate))
    end

    FinanceCore.discount(c::Constant,t) = FinanceCore.discount(c.rate,t)
    FinanceCore.rate(c::Constant,t) = c.rate
end

function pv(model::M,c::AbstractContract) where {M<:Yield.AbstractYieldModel}
    p = Projection(model,CashflowProjection(),c)
    cfs = collect(p)
    return mapreduce(+,cfs) do cf
        FinanceCore.discount(model,cf.time) * cf.amount
    end
end
