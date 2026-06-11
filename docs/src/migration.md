# Migration Guide

## v5.x to v6

- **Continuous zero rates are the curve primitive.** Curve composition and shift arithmetic (`+`, `-`, `*`, `/`, `TenorShift`, `ProjectedShift`) operate in continuous-zero-rate space, which is equivalent to multiplying/dividing/exponentiating discount factors. See [Yield Curve Arithmetic](@ref).
- **[Breaking] `Spline.Linear()`, `Spline.Quadratic()`, and `Spline.Cubic()` now default to local polynomial interpolants** instead of B-splines. Use `Spline.BSpline(order)` to retain the previous (global) B-spline behavior.
- **`ZeroRateCurve` eagerly builds its interpolation at construction** rather than on first evaluation.
- **Bootstrap `fit` (`Fit.Bootstrap()`) is now an exact per-knot root-solve.** Results match the previous optimizer-based bootstrap within the old optimizer's tolerance, but each quote is now repriced exactly.
- **`TransformedYield` is deprecated — use `Yield.TenorShift`.** The old name remains available as a `Base.@deprecate_binding` alias but will be removed in a future release.
- **`MonotoneConvex` forward rates are now continuous at and beyond the last knot** (extrapolation is anchored at the boundary instantaneous forward; see `Yield.instantaneous_forward`), and the Hagan-West positivity collar has been corrected (it is now generalized to handle negative discrete forwards).
- **`par` now throws an informative `ArgumentError`** when the requested maturity implies a stub period that cannot be represented with the given coupon frequency.
- **The Makie plotting extension now requires Makie ≥ 0.24.**
- **With FinanceCore v3, `irr` / `internal_rate_of_return` return `Periodic(NaN, 1)` instead of `nothing`** when no root is found. Replace `isnothing(irr(x))` checks with `isnan(rate(irr(x)))`.

## v5.4 to v5.5

### `TransformedYield` renamed to `TenorShift`; new `ProjectedShift`

`Yield.TransformedYield` has been renamed to [`Yield.TenorShift`](@ref FinanceModels.Yield.TenorShift) to sit alongside the new [`Yield.ProjectedShift`](@ref FinanceModels.Yield.ProjectedShift), which adds a second time axis (projection / as-of time) to the shift rule. Both are concrete subtypes of the new [`Yield.AbstractYieldShift`](@ref FinanceModels.Yield.AbstractYieldShift).

Use [`ProjectedShift`](@ref FinanceModels.Yield.ProjectedShift) for shifts whose shape evolves across a projection horizon (BMA SBA phase-ins, IFRS17 macro scenarios, EV runoffs). See the Yield Shifts section in [Available Models - Yields](@ref) for usage.

!!! warning "Breaking changes shipped under a minor bump"
    This release is tagged minor (5.4 → 5.5) but contains two breaking behavior changes
    that downstream code may need to react to:

    1. **Field rename: `.transform` → `.rule`.** Direct field access on
       `TransformedYield` instances (e.g., `ty.transform`) will fail. The
       `TransformedYield` type name itself is preserved via
       `const TransformedYield = TenorShift`, so constructor and `+`-operator
       call sites continue to work unchanged.
    2. **Strict `Rate` return contract.** `Base.zero` on `TenorShift` /
       `ProjectedShift` now type-asserts the rule's return value as
       `FinanceCore.Rate`. Rules that previously returned a plain `Real`
       (silently coerced to `Continuous`) will now raise a `TypeError` at
       call time. Replace `(z, t) -> z.continuous_value + 0.01` with
       `(z, t) -> Continuous(z.continuous_value + 0.01)`, or more
       idiomatically `(z, t) -> z + Continuous(0.01)` and let `Rate`
       arithmetic carry compounding convention.

    The `TransformedYield` alias is slated for removal one minor release after
    introduction. The `+` operator semantics (`curve + (z, t) -> Rate`) are
    unchanged — only the returned struct's name changes.

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
model = Yields.Par(SmithWilson(...), rates,timepoints)
```

Now, follow the pattern of:

1. Define the quotes you want to fit the model to
2. `fit` the model to those quotes

Example:

```julia
quotes = ParYield.(rates,timepoints)
model = fit(Yield.SmithWilson(ufr=0.03, α=0.1), quotes)
```

Note that `SmithWilson` is not exported at the top level (qualify it as `Yield.SmithWilson`) and that the `ufr` and `α` keyword arguments are required: they are model hyperparameters that are not solved for in the fit.

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
