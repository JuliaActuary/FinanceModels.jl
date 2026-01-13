abstract type AbstractModel end
Base.Broadcast.broadcastable(x::T) where {T <: AbstractModel} = Ref(x)

"""
    NullModel()
A singleton type representing a placeholder model for when you don't really need a model. For example: determining nominal cashflows for fixed income contract.
"""
struct NullModel <: AbstractModel end

# useful for round-tripping or iterating on quotes?
function FinanceCore.Quote(m::M, c::C) where {M <: AbstractModel, C <: FinanceCore.AbstractContract}
    return FinanceCore.Quote(pv(m, c), c)
end

include("Spline.jl")
include("Yield.jl")
include("Volatility.jl")
include("Equity.jl")

"""
    present_value(model,contract,current_time=0.0)
    present_value(model,projection,current_time=0.0)
 
 Return the value of the contract as corresponding with the valuation assumptions embedded in the `model` for the given contract or projection with `CashflowProjection` kind.

# Examples

```julia
m = Equity.BlackScholesMerton(0.01, 0.02, 0.15)

a = Option.EuroCall(CommonEquity(), 1.0, 1.0)

pv(m, a) # ≈ 0.05410094201902403
```
"""
function FinanceCore.present_value(model, c::FinanceCore.AbstractContract, cur_time = 0.0)
    p = Projection(c, model, CashflowProjection())
    xf = p |> Filter(cf -> cf.time >= cur_time) |> Map(cf -> FinanceCore.discount(model, cur_time, cf.time) * cf.amount)
    return foldxl(+, xf)
end
