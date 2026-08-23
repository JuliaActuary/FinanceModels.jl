# Migration Guide

## v6.4 to v6.5

!!! warning "New errors from `ZeroRateCurve` construction"
    `ZeroRateCurve` now owns its data and validates its inputs. Inputs that previously produced silently wrong curves (`NaN` discount factors from duplicate tenors, misordered knots with `MonotoneConvex`, stale interpolation caches after `@set`) are now loud `ArgumentError`s or cannot happen at all.

- **`ZeroRateCurve` copies `rates`/`tenors` into owned storage** and promotes them to one concrete floating-point type (`Int` tenors become `Float64`; `Float32`+`BigFloat` → `BigFloat`; `Float64`+`ForwardDiff.Dual` → `Dual`; ranges and tuples are accepted). Mutating the vectors you passed in no longer affects the curve.
- **`zrc.rates` and `zrc.tenors` are read-only.** Indexed assignment (`zrc.rates[1] = …`, `.=`, `sort!`, writes through `view`) throws an `ArgumentError`. Use `copy(zrc.rates)` for a mutable copy, and `Accessors.@set zrc.rates[i] = x` / `@set zrc.tenors = grid` / `@set zrc.spline = …` to derive a modified curve — these now re-validate and **rebuild the interpolation cache** (previously `@set` silently kept the stale interpolant, so two curves could be `==` yet price differently).
- **The interpolation cache is internal**: `propertynames(zrc)` is `(:rates, :tenors, :spline)`, `zrc._model` throws, and `@set zrc._model = …` throws. The undocumented 4-argument `ZeroRateCurve(rates, tenors, spline, model)` constructor is removed; `ConstructionBase.setproperties`/`constructorof` reconstruct correctly through the 3-argument form.
- **Validation** — `ArgumentError` when `rates` and `tenors` differ in length or are empty; when any rate or tenor is non-finite (`NaN`, `±Inf`); when any tenor is negative; when tenors are not strictly increasing (unsorted or duplicated); or when there are fewer knots than the interpolant needs — **`Spline.PCHIP()` and `Spline.Akima()` 3; `Spline.Linear()`/`Quadratic()`/`Cubic()`/`BSpline(n)` 2; `Spline.MonotoneConvex()` 1**. `Spline.BSpline(d)`/`Spline.PolynomialSpline(d)` now require `d ≥ 1` at construction. Negative rates and a `t = 0` knot in the direct form remain valid.
- **The direct form no longer sorts unsorted tenors** (it throws); the `ZeroRateCurve(curve::AbstractYieldModel, tenors)` sampling form still sorts its tenor grid, still requires `t > 0`, and now samples through `zero(curve, t)` instead of `-log(discount(curve, t))/t`, which is numerically stable at very small and very large tenors.
- `isequal(::ZeroRateCurve, ::ZeroRateCurve)` is now fieldwise (previously it fell back to `==`), so `isequal(a, b)` implies `hash(a) == hash(b)`; `ZeroRateCurve([-0.0], t)` and `ZeroRateCurve([0.0], t)` are `==` but not `isequal`, exactly like the underlying arrays.
- `fit(zrc::ZeroRateCurve, quotes)` is now supported.

!!! warning "`Yield.Spline` and `Yield.MonotoneConvex` share the same storage and validation"
    The two curves `ZeroRateCurve` is built on had the same aliasing and stale-cache exposure. Both now go through the shared knot-grid construction, so direct construction, `ZeroRateCurve`, `Yield.build_model`, `fit` and bootstrap raise identical `ArgumentError`s for the same invalid grids.

- **`Yield.Spline(spline, tenors, rates)` copies its inputs** (DataInterpolations previously kept the caller's `Vector`s for `Spline.Linear()`/`Quadratic()`/`Cubic()`, so mutating them later changed the curve while its cached coefficients went stale) and applies the validation list above (length, empty, non-finite, negative/unsorted/duplicate tenors, minimum knots).
- **`Yield.MonotoneConvex(rates, times)` accepts any iterable of reals** (ranges, tuples, `Int`s — previously `Vector` only), copies and promotes them, and applies the same validation. **`c.rates`, `c.times`, `c.f` and `c.fᵈ` are read-only** (indexed assignment throws; `copy(c.rates)` for a mutable copy) so the cached forwards can no longer be desynchronised from the knots. `propertynames(c)` lists `(:rates, :times)`; the derived `f`/`fᵈ` are still readable and still listed by `propertynames(c, true)`. `@set c.rates[i] = …` / `@set c.times = …` re-validate and recompute the forwards; `@set c.f = …` / `@set c.fᵈ = …` throw instead of being silently discarded.
- **Fitting knot curves is one rebuild per optimizer candidate.** `__default_optic` for `Yield.MonotoneConvex` and `ZeroRateCurve` is now a single batch `FinanceModels.KnotRatesOptic()` over all knot rates (`Accessors.getall` → tuple of rates, `setall` → one reconstruction) instead of one `@optic(_.rates[i])` per knot, which rebuilt the curve once per knot per candidate — O(n²) in the number of knots (measured 5×/50×/500× less allocation per candidate at 5/50/500 knots). If you passed your own `variables = (@optic(_.rates[i]) => …, …)` to `fit`, that still works; `length(FinanceModels.__default_optic(curve))` is now `1`.
- **`fit` validates the knot grid before optimising**: `fit(::Spline.SplineCurve, quotes, Fit.Loss(…))`, `fit(Yield.MonotoneConvex(), quotes)` and `Fit.Bootstrap()` throw the construction `ArgumentError` for duplicate/non-positive maturities or too few quotes for the interpolant (e.g. two quotes with `Spline.PCHIP()`) up front, rather than failing inside the solver. Optimizer trial curves are built without validation (a non-finite candidate is a bad loss, not an exception); the returned curve is validated, so a fit that diverged to non-finite rates now throws instead of returning a curve of `NaN`s.

## v6.0 to v6.1

!!! warning "Changed numbers and new errors"
    Several items below change computed values (curve extrapolation, fitted bootstrap curves where the prior optimizer had not fully converged) or convert previously-silent mispricing into loud errors. Review each against your pipelines before upgrading.

- **`MonotoneConvex` (the default `ZeroRateCurve` interpolant) — two value-changing corrections:**
  - *Forward rates are now continuous at and beyond the last knot*: extrapolation is anchored at the boundary instantaneous forward `f(tₙ)` instead of the last discrete forward (see `Yield.instantaneous_forward`). **Extrapolated zero rates change** — on a typical upward-sloping curve with a 10y last knot, the 20y zero moves on the order of +10bp (about −2% PV for a 20y cashflow). For steeply inverted/humped curves the boundary forward can be collared to 0, giving a 0% forward tail beyond the last knot — extend your knot grid past your longest cashflow if you discount far beyond it.
  - *The Hagan-West positivity collar was corrected* (it previously clamped the wrong nodes and left one node unclamped, so the guaranteed-positive-forwards property could fail). Fitted/interpolated values change only where a clamp binds (sharply non-monotone forward curves); the collar is also generalized to negative discrete forwards.
  - The module-local `Yield.forward(mc::MonotoneConvex, t)` (instantaneous forward) was renamed `Yield.instantaneous_forward(mc, t)`. `Yield.forward` now refers to `FinanceCore.forward`, so the *same call* returns the discrete one-period forward as a `Rate` — update qualified callers.
- **Bootstrap `fit` (`Fit.Bootstrap()`) is now an exact per-knot root-solve** instead of a per-knot optimizer pass. For zero-coupon quotes (any interpolant) and for coupon quotes with *local* interpolants (`Spline.Linear/Quadratic/Cubic`), every quote is repriced to root-finder precision; with *global* interpolants (`Spline.BSpline`) later knots still reshape earlier segments, so earlier coupon quotes reprice approximately (comparable to the previous behavior). Quotes are now sorted by maturity internally; duplicate maturities are an error.
- **`ZeroRateCurve` eagerly builds its interpolation at construction** rather than on first evaluation, and **`discount(zrc, t)` for `t < 0` now throws a `DomainError`** (it previously returned `1.0` silently — a misprice for anything that actually discounted at negative times). Notably, a `Bond.Floating` whose maturity is not an integer multiple of the coupon period generates a stub first coupon that references `forward(model, t - 1/freq, t)` with a *negative* start time: on a `ZeroRateCurve` this was previously a silent half-sized stub forward and is now a loud error. Align floater maturities/resets to the coupon period.
- **`par` now throws an informative `ArgumentError`** when the requested maturity implies a stub period that cannot be represented with the given coupon frequency (previously a bare `InexactError`).
- **`TransformedYield` is deprecated — use `Yield.TenorShift`.** The old name remains available as a `Base.@deprecate_binding` alias but will be removed in a future release.
- **The Makie plotting extension now targets Makie ≥ 0.24 directly** (the previous MakieCore-based extension stopped loading when Makie 0.24 absorbed MakieCore, so plot recipes had been silently unavailable). Makie < 0.24 is no longer supported. The UnicodePlots extension now renders only for rich (`text/plain` MIME) display; `print`/string interpolation of curves no longer embeds a chart.
- **With FinanceCore v3, `irr` / `internal_rate_of_return` return `Periodic(NaN, 1)` instead of `nothing`** when no root is found. Replace `isnothing(irr(x))` checks with `isnan(rate(irr(x)))`.

## v5.x to v6

- **Continuous zero rates are the curve primitive.** Curve composition and shift arithmetic (`+`, `-`, `*`, `/`, `TenorShift`, `ProjectedShift`) operate in continuous-zero-rate space, which is equivalent to multiplying/dividing/exponentiating discount factors. See [Yield Curve Arithmetic](@ref).

## v5.4 to v5.5

### `TransformedYield` renamed to `TenorShift`; new `ProjectedShift`

`Yield.TransformedYield` has been renamed to [`Yield.TenorShift`](@ref FinanceModels.Yield.TenorShift) to sit alongside the new [`Yield.ProjectedShift`](@ref FinanceModels.Yield.ProjectedShift), which adds a second time axis (projection / as-of time) to the shift rule. Both are concrete subtypes of the new [`Yield.AbstractYieldShift`](@ref FinanceModels.Yield.AbstractYieldShift).

Use [`ProjectedShift`](@ref FinanceModels.Yield.ProjectedShift) for shifts whose shape evolves across a projection horizon (BMA SBA phase-ins, IFRS17 macro scenarios, EV runoffs). See the Yield Shifts section in [Available Models - Yields](@ref) for usage.

!!! warning "Breaking changes shipped under a minor bump"
    This release is tagged minor (5.4 → 5.5) but contains two breaking behavior changes
    that downstream code may need to react to:

    1. **Field rename: `.transform` → `.rule`.** Direct field access on
       `TransformedYield` instances (e.g., `ty.transform`) will fail. The
       `TransformedYield` type name itself was preserved in v5.5 via
       `const TransformedYield = TenorShift`, so constructor and `+`-operator
       call sites continue to work unchanged. (As of v6 the alias emits a
       deprecation warning — see the v5.x → v6 section above.)
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
