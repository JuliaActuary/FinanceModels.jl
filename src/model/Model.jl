abstract type AbstractModel end

# a model for when you don't really need a model
# (e.g. determining nominal cashflows for fixed income contract)
struct NullModel <: AbstractModel end

# useful for round-tripping or iterating on quotes?
function Quote(m::M, c::C) where {M<:AbstractModel,C<:AbstractContract}
    return Quote(pv(m, c), c)
end

include("Yield.jl")
include("Volatility.jl")
include("Equity.jl")