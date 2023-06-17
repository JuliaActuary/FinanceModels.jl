abstract type AbstractModel end
Base.Broadcast.broadcastable(x::T) where {T<:AbstractModel} = Ref(x)

# a model for when you don't really need a model
# (e.g. determining nominal cashflows for fixed income contract)
struct NullModel <: AbstractModel end

# useful for round-tripping or iterating on quotes?
function Quote(m::M, c::C) where {M<:AbstractModel,C<:AbstractContract}
    return Quote(pv(m, c), c)
end

include("Spline.jl")
include("Yield.jl")
include("Volatility.jl")
include("Equity.jl")

function FinanceCore.present_value(model::M, c::AbstractContract; cur_time=0.0) where {M<:Yield.AbstractYieldModel}
    p = Projection(c, model, CashflowProjection())
    xf = p |> Filter(cf -> cf.time >= cur_time) |> Map(cf -> FinanceCore.discount(model, cur_time, cf.time) * cf.amount)
    foldxl(+, xf)
end