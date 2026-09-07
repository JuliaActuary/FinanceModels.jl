"""
    model_requirements(contract) -> Tuple

Return the `key => model_type` requirements for a cashflow projection. Fixed
cashflows have no requirements; floating bonds require a yield model, and
`FX.Converted` wrappers require an FX model as well as their inner requirements.
Composites, forwards, portfolios, and Transducers eductions retain the requirements
of their underlying contracts. A key may occur more than once: its model must
satisfy every occurrence.

Custom projectable contracts should implement this function, returning `()` only
when they require no model. There is no model-independent fallback for unknown
contracts. Explicit `Projection(contract, model_store)` remains available without
implementing this convenience protocol.
"""
function model_requirements end

model_requirements(::Cashflow) = ()
model_requirements(::Bond.Fixed) = ()
model_requirements(::FX.BasisSwapLeg) = ()
model_requirements(c::Bond.Floating) = (c.key => Yield.AbstractYieldModel,)
model_requirements(c::FinanceCore.Composite) = (model_requirements(c.a)..., model_requirements(c.b)...)
model_requirements(c::Forward) = model_requirements(c.instrument)
model_requirements(c::Transducers.Eduction) = model_requirements(c.coll)
model_requirements(c::FX.Converted) = (c.key => FX.AbstractFXModel, model_requirements(c.contract)...)
model_requirements(cs::AbstractArray) = Tuple(req for c in cs for req in model_requirements(c))

"""
    Projection(contract; index)

Project a contract or portfolio using `index` for every required model key.
The model must satisfy all [`model_requirements`](@ref). For contracts with
different index curves or an FX conversion, pass an explicit model store instead:
`Projection(contract, Dict("SOFR" => sofr, "EURUSD" => fx))`.

This form is useful inside a valuation closure: rebuilding the projection with a
bumped `index` recomputes floating coupons, including transformed swap legs.
Omitting `index` retains the default model-free projection.

```julia
curve = Yield.Constant(0.04)
swap = InterestRateSwap(curve, 5.0)
value(index, discount_curve) = present_value(discount_curve, Projection(swap; index))
value(curve, curve) # approximately zero
```
"""
function Projection(c; index = nothing)
    isnothing(index) && return Projection(c, NullModel(), CashflowProjection())
    requirements = model_requirements(c)
    for (key, model_type) in requirements
        index isa model_type || throw(
            ArgumentError(
                "Projection: key $(repr(key)) requires $model_type, but index is $(typeof(index)). " *
                    "Pass an explicit model store with Projection(contract, models)."
            )
        )
    end
    return Projection(c, Dict(key => index for (key, _) in requirements), CashflowProjection())
end
