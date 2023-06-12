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

FinanceCore.discount(c::Constant, t) = FinanceCore.discount(c.rate, t)
FinanceCore.rate(c::Constant, t) = c.rate
end

function value(model::M, c::AbstractContract; cur_time=0.0) where {M<:Yield.AbstractYieldModel}
    p = Projection(model, c, CashflowProjection())
    foldxl(+, Map(cf -> FinanceCore.discount(model, cf.time - cur_time) * cf.amount), p)
end

abstract type AbstractVolatilityModel <: AbstractModel end

struct ConstantVolatility{T} <: AbstractVolatilityModel
    σ::T
end

abstract type AbstractEquityModel <: AbstractModel end

struct BlackScholesMerton{T,U,V} <: AbstractEquityModel
    r::T # risk free rate
    q::U # dividend yield
    σ::V # roughly equivalent to the volatility in the usual lognormal model multiplied by F^{1-β}_{0}
end

function get_volatility(vol::ConstantVolatility, strike_ratio, time_to_maturity)
    return vol.σ
end

function value(model::M, c::Option.EuroCall{Equity,K,T}) where {M<:BlackScholesMerton,K,T}
    eurocall(; S=1.0, K=c.strike, τ=c.maturity, r=model.r, q=model.q, σ=model.σ)

end