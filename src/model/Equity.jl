"""
The `Equity` module provides equity-related model definitions.

See also: the [`Volatility`](@ref) module.
"""
module Equity
import ..AbstractModel

abstract type AbstractEquityModel <: AbstractModel end

"""
    BlackScholesMerton(r, q, σ) <: AbstractEquityModel

A struct representing the Black-Scholes-Merton model for equity prices.

# Arguments
- `r`: The risk-free rate.
- `q`: The dividend yield.
- `σ`: The volatility model of the underlying asset (see [`Volatility`](@ref) module) 

# Fields
- `r`: The risk-free rate.
- `q`: The dividend yield.
- `σ`: The volatility model of the underlying asset (see [`Volatility`](@ref) module)

When [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss)ting, the volatility will be solved-for; volatility itself is a sub-model that will be optimized with a default optimization bound of `0.0 .. 10.0`

# Examples
```julia-repl
julia> model = BlackScholesMerton(0.05, 0.02, 0.2)
BlackScholesMerton{Float64, Float64, Float64}(0.05, 0.02, 0.2)
```

Valuing an option:
```julia
m = Equity.BlackScholesMerton(0.01, 0.02, 0.15)

a = Option.EuroCall(CommonEquity(), 1.0, 1.0)

@test pv(m, a) ≈ 0.05410094201902403
```

Fitting a set of option prices:

```julia
qs = [
    Quote(0.0541, a),
    Quote(0.072636, b),
]
m = Equity.BlackScholesMerton(0.01, 0.02, Volatility.Constant())
fit(m, qs)
@test fit(m, qs).σ ≈ 0.15 atol = 1e-4

```
"""
struct BlackScholesMerton{T,U,V} <: AbstractEquityModel
    r::T # risk free rate
    q::U # dividend yield
    σ::V # roughly equivalent to the volatility in the usual lognormal model multiplied by F^{1-β}_{0}
end

end

"""
    volatility(volatiltiy_model,strike_ratio,time_to_maturity)

Returns the volatility associated with the moneyness (strike/price ratio) and time to maturity.
"""
function volatility(vol::Volatility.Constant, strike_ratio, time_to_maturity)
    return vol.σ
end

function FinanceCore.present_value(model::M, c::Option.EuroCall{CommonEquity,K,T}) where {M<:Equity.BlackScholesMerton,K,T}
    eurocall(; S=1.0, K=c.strike, τ=c.maturity, r=model.r, q=model.q, σ=model.σ)

end