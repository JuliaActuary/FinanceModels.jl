
module Volatility
import ..AbstractModel

abstract type AbstractVolatilityModel <: AbstractModel end

struct Constant{T} <: AbstractVolatilityModel
    σ::T
end
Constant() = Constant(0.0)
end