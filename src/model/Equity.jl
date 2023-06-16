
module Equity
import ..AbstractModel

abstract type AbstractEquityModel <: AbstractModel end

struct BlackScholesMerton{T,U,V} <: AbstractEquityModel
    r::T # risk free rate
    q::U # dividend yield
    σ::V # roughly equivalent to the volatility in the usual lognormal model multiplied by F^{1-β}_{0}
end
end

function volatility(vol::Volatility.Constant, strike_ratio, time_to_maturity)
    return vol.σ
end

function FinanceCore.present_value(model::M, c::Option.EuroCall{CommonEquity,K,T}) where {M<:Equity.BlackScholesMerton,K,T}
    eurocall(; S=1.0, K=c.strike, τ=c.maturity, r=model.r, q=model.q, σ=model.σ)

end