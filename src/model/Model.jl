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

# The closed-form value of a contract under the model that prices it (see `valuation_model`): the
# model files below add their formulas as methods.
function __closed_form end

include("Spline.jl")
include("Yield.jl")
include("Volatility.jl")
include("Equity.jl")
include("FX.jl")
include("Stochastic.jl")
