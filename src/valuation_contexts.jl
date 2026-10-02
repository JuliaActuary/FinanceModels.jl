"""
    Models(model, store)
    Models(model; index)

A valuation context that holds the models a contract observes as well as the model that values it.
`model` discounts (`discount(ctx, t)`) and prices (`valuation_model(ctx)`); the models a contract
reads by key, such as the index curve of a `Bond.Floating` or the FX model of an `FX.Converted` leg,
come from `store` (`ctx[key]`), any collection indexed by key (a `Dict` or a `NamedTuple`).
`Models(model; index)` reads `index` under every key.

Any model is a valuation context of its own: `present_value(curve, bond)` discounts on `curve`, and
`present_value(bsm, call)` prices on `bsm`. `Models` adds the observed models:

```julia
curve = Yield.Constant(0.04)
swap = InterestRateSwap(curve, 5.0; frequency = 1)
present_value(Models(curve; index = curve), swap)   # approximately zero
present_value(Models(ois, Dict("SOFR" => sofr, "EURUSD" => fx)), portfolio)
collect(Projection(swap, Models(curve; index = curve)))   # its cashflows
```
"""
struct Models{M, S} <: AbstractModel
    model::M
    store::S
end
Models(model; index) = Models(model, __EveryKey(index))

# The store of `Models(model; index)`: `index` under every key.
struct __EveryKey{I}
    index::I
end
Base.getindex(s::__EveryKey, key) = s.index

Base.getindex(m::Models, key) = m.store[key]
FinanceCore.discount(m::Models, t...) = FinanceCore.discount(m.model, t...)

"""
    valuation_model(ctx)

The model that values contracts in the valuation context `ctx`: the `model` of a [`Models`](@ref)
context, and any other model itself. A contract with a closed-form value reads the model it needs
through it, so that the formula applies alike to `present_value(model, c)` and to
`present_value(Models(model, store), c)`:

```julia
FinanceCore.present_value(ctx, c::MyOption) = my_formula(valuation_model(ctx), c)
my_formula(m::MyModel, c::MyOption) = ...
```
A context whose model does not price the contract fails in the formula, with a `MethodError`.
"""
valuation_model(m::Models) = m.model
valuation_model(m) = m

# FinanceModels' closed-form contracts, valued by the context's model (`__closed_form`).
const __ClosedFormContract = Union{
    Option.EuroCall, Option.EuroPut, Option.ZCBCall, Option.ZCBPut, Option.Cap, Option.Floor, Option.Swaption,
}
FinanceCore.present_value(ctx, c::__ClosedFormContract) = __closed_form(valuation_model(ctx), c)
