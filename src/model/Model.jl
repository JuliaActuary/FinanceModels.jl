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

"""
    closed_form(model, contract)

The value of `contract` under `model` by a formula that depends on the model: Black–Scholes for a
European option under `Equity.BlackScholesMerton`, Jamshidian's decomposition for a swaption under
`ShortRate.Vasicek` or `ShortRate.HullWhite`.

`present_value(ctx, contract)` is the extension point for valuing a contract, and `closed_form` is
the kernel it calls when the formula depends on the type of the pricing model. FinanceModels'
options, caps, floors and swaptions are valued by
`present_value(ctx, c) = closed_form(valuation_model(ctx), c)`, so a model prices them by adding a
`closed_form` method. A custom contract opts in the same way:

```julia
# pays S_T - K at T on a unit stock
struct EquityForward <: FinanceCore.AbstractContract
    strike::Float64
    maturity::Float64
end
FinanceCore.present_value(ctx, c::EquityForward) = FinanceModels.closed_form(valuation_model(ctx), c)
FinanceModels.closed_form(m::Equity.BlackScholesMerton, c::EquityForward) =
    exp(-m.q * c.maturity) - c.strike * exp(-m.r * c.maturity)
```

It then values under the model, under [`Models`](@ref), in a portfolio and in a `Composite`:

```julia
m = Equity.BlackScholesMerton(0.03, 0.01, 0.2)
fwd = EquityForward(1.0, 2.0)
call = Option.EuroCall(CommonEquity(), 1.0, 2.0)
present_value(m, fwd)
present_value(Models(m; index = Yield.Constant(0.03)), fwd)   # the same value
present_value(m, [fwd, call])                                 # the sum of the two values
present_value(m, Composite(fwd, call))                        # the same sum
```

A context whose model has no method for the contract throws a `MethodError`. A formula that reads
only discount factors and models by key (`discount(ctx, t)`, `ctx[key]`) needs no model kernel: it
is written on `present_value` directly.
"""
function closed_form end

include("Spline.jl")
include("Yield.jl")
include("Volatility.jl")
include("Equity.jl")
include("FX.jl")
include("Stochastic.jl")
