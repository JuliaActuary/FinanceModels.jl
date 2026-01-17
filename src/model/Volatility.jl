module Volatility
import ..AbstractModel

abstract type AbstractVolatilityModel <: AbstractModel end

"""
    Volatility.Constant(σ)

A constant volatility per period. If σ is not explicitly passed, then it is set to zero.
"""
struct Constant{T} <: AbstractVolatilityModel
    σ::T
end
Constant() = Constant(0.0)

# Define approximate comparison for Volatility.Constant with numbers
Base.isapprox(a::Constant, b::Real; kwargs...) = isapprox(a.σ, b; kwargs...)
Base.isapprox(a::Real, b::Constant; kwargs...) = isapprox(a, b.σ; kwargs...)
Base.isapprox(a::Constant, b::Constant; kwargs...) = isapprox(a.σ, b.σ; kwargs...)

end
