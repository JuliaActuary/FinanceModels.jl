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
end
