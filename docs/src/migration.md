# Migration Guide

## Within v5

### `TransformedYield` renamed to `TenorShift`; new `ProjectedShift`

`Yield.TransformedYield` has been renamed to [`Yield.TenorShift`](@ref FinanceModels.Yield.TenorShift) to sit alongside the new [`Yield.ProjectedShift`](@ref FinanceModels.Yield.ProjectedShift), which adds a second time axis (projection / as-of time) to the shift rule. Both are concrete subtypes of the new [`Yield.AbstractYieldShift`](@ref FinanceModels.Yield.AbstractYieldShift).

`Yield.TransformedYield` is retained as a deprecation alias (`const TransformedYield = TenorShift`) and will be removed one minor release after introduction. The field name `transform` has been renamed to `rule`; if downstream code accesses `.transform` directly it must switch to `.rule`. The `+` operator semantics (`curve + (z, t) -> ...`) are unchanged — only the returned struct's name changes.

Use [`ProjectedShift`](@ref FinanceModels.Yield.ProjectedShift) for shifts whose shape evolves across a projection horizon (BMA SBA phase-ins, IFRS17 macro scenarios, EV runoffs). See the Yield Shifts section in [Available Models - Yields](@ref) for usage.

## v4 to v5

### Yield curve `+` and `-` now operate in continuous zero-rate space

In v4, `curve_a + curve_b` added rates in whatever compounding convention the curves happened to use. In v5, `+` and `-` always work in **continuous zero-rate (CZR) space**, which is equivalent to multiplying/dividing discount factors:

```julia
# v5 behavior:
combined = curve_a + curve_b
discount(combined, t) == discount(curve_a, t) * discount(curve_b, t)
```

This is the economically correct way to combine deflators — see [Yield Curve Arithmetic](@ref) for a full explanation.

**What to check when upgrading:** If your v4 code added curves whose rates were expressed in `Periodic` conventions, the combined discount factors will now differ by the cross-term. For small rates and short horizons the difference is minor, but it compounds over long projections (e.g. 10 bps/year for a 5% base + 2% spread).

### `ForwardYields` renamed to `ForwardYield`

The plural `ForwardYields` has been renamed to `ForwardYield` for consistency with other singular type names (`Yield.Constant`, `ZCBYield`, etc.).

## v3 to v4

### Yields.jl is now FinanceModels.jl

This re-write accomplishes three primary things:

- Provide a composable set of **contracts** and **`Quotes`**
- Those contracts, when combined with a **model** produce a **`Cashflow`** via a flexibly defined `Projection`
- **models** can be `fit` with a new unified API: `fit(model_type,quotes,fit_method)`

### Migrating Code

#### Update Dependencies

You should remove `Yields` from your project's dependencies and add `FinanceModels` instead. ([link to Pkg documentation on how to do this](https://pkgdocs.julialang.org/v1/managing-packages/))

#### API Changes

Previously, the API pattern was, e.g.:

```julia
model = Yields.Par(SmitWilson(...), rates,timepoints)
```

Now, follow the pattern of:

1. Define the quotes you want to fit the model to
2. `fit` the model to those quotes

Example:

```julia
quotes = ParYield.(rates,timepoints)
model = fit(SmithWilson(),quotes)
```

#### Details of changes

Previously the kind of contract, the implied quotes, the type of model, and how the fitting process worked were all combined into a single call (`Yields.Par`). This minimized the amount of code needed to construct a yield curve, but left it fairly cumbersome to extend the package. For example, for every new yield curve model, methods for `Par`, `CMT`, `OIS`, `Zero`, ... had to be defined. Additionally, all of the inputs needed to be yields - specifying a price was not available as an argument to fit.

With the new design of the package, creating a completely new model is much easier, as only the model itself and the valuation primitives need to be defined. For example, defining a new yield curve type that works to value contracts instrument quotes only requires defining the `discount` method. To allow the model to be `fit` requires only defining a default set of parameters to optimize with `__default_optic`:

```julia
 using FinanceModels, FinanceCore
 using AccessibleModels 
 using IntervalSets
 
struct ABDiscountLine{A} <: FinanceModels.Yield.AbstractYieldModel
    a::A
    b::A
end

# define the default constructor for convenience
ABDiscountLine() = ABDiscountLine(0.,0.)

function FinanceCore.discount(m::ABDiscountLine,t)
    #discount rate is approximated by a straight lined, floored at 0.0 and capped at 1.0
    clamp(m.a*t + m.b, 0.0,1.0) 
end


# `@optic` indicates what in our model variables needs to be updated (from AccessibleModels.jl)
# `-1.0 .. 1.0` says to bound the search from negative to positive one (from IntervalSets.jl)
FinanceModels.__default_optic(m::ABDiscountLine) = (
    @optic(_.a) => -1.0 .. 1.0,
    @optic(_.b) => -1.0 .. 1.0,
)

quotes = ZCBPrice([0.9, 0.8, 0.7,0.6])

m = fit(ABDiscountLine(),quotes)
```
