# Migration Guide

## v6.x to v7.0

- **`Fit.Bootstrap()` supports only `Spline.Linear()`** (equivalently
  `Spline.PolynomialSpline(1)` or `Spline.BSpline(1)`). Bootstrapping solves one
  quote at a time and requires that each new knot leave earlier curve segments
  unchanged. With quadratic, cubic, higher-order B-spline, PCHIP, or Akima
  interpolation a later knot reshapes earlier segments, so earlier coupon quotes
  silently stopped repricing: on uneven par quotes the residuals reached 0.40%
  (quadratic) and 0.11% (cubic). These strategies now throw an `ArgumentError`.
  **Migration:** use `Spline.Linear()` with `Fit.Bootstrap()`, or fit the smoother
  strategy to all quotes at once with `Fit.Loss(x -> x^2)`; fit a monotone convex
  curve with `fit(Spline.MonotoneConvex(), quotes)` (bootstrap previously switched
  to this loss fit silently). Zero-coupon quote sets were exact with every strategy,
  so their fitted curves change only if you switch to a loss fit.
- Bootstrap validates its inputs up front: empty quote sets, non-finite or
  non-positive maturities, and duplicate maturities throw an `ArgumentError`.
  After solving, every quote is repriced on the returned curve, and a residual
  beyond root-finder precision throws.
- Full-curve `Fit.Loss` spline fits start from slightly sloped rates near 5%
  instead of a flat 5%, which lets PCHIP and Akima fits converge. Converged fits
  of other strategies move by at most about 1e-9 in zero rate.

### `ZeroRateCurve` returns the curve it builds

`ZeroRateCurve` is now a construction function rather than a type. It returns a
`Yield.MonotoneConvex` for `Spline.MonotoneConvex()` (the default) and a `Yield.Spline` for every
other method. Both are subtypes of the new `Yield.AbstractInterpolatedZeroCurve`, and so are the
curves returned by spline `fit`s and `Fit.Bootstrap()`. All of them share one interface:

| v6 | v7 |
|:---|:---|
| `zrc isa ZeroRateCurve`, `f(z::ZeroRateCurve)` | `zrc isa Yield.AbstractInterpolatedZeroCurve` (or `Yield.AbstractYieldModel` when only discounting) |
| `zrc.rates`, `zrc.tenors`, `mc.times` | `knot_rates(curve)`, `knot_tenors(curve)` (read-only vectors) |
| `@set zrc.rates[2] = 0.031` | `reconstruct(curve; rates = new_rates)` (`@set` on `rates`, `tenors`, `spline` or `extrapolation` still works and calls `reconstruct`) |
| `ZeroRateCurve(dual_rates, zrc.tenors, zrc.spline)` in a gradient | `reconstruct(curve; rates = dual_rates)` |
| `Yield.build_model(spline, tenors, rates; extrapolation)` | `ZeroRateCurve(rates, tenors, spline; extrapolation)` (note the argument order) |
| `fit(Yield.MonotoneConvex(), quotes)` | `fit(Spline.MonotoneConvex(), quotes)` |
| `mc.f`, `mc.fᵈ` | `Yield.instantaneous_forward(mc, t)`; the node forwards are internal |
| `Yield.Spline(fn)` for a callable zero-rate function | `Yield.Constant(0.0) + ((z, t) -> Continuous(fn(t)))` |

- **`Yield.build_model` is removed**, and so are the `Yield.MonotoneConvex()` placeholder and
  its callable form (#272). `Spline.MonotoneConvex()` is the only monotone convex selector: it
  works with `ZeroRateCurve`, loss `fit`s (whose default optimizer for it is `LBFGS()`), and
  `FX.Forwards`, and every interpolation method is built by the same code wherever the curve
  comes from. `Yield.Spline(Spline.MonotoneConvex(), …)` throws an `ArgumentError` pointing to
  `ZeroRateCurve`.
- **Every knot curve owns and validates its data.** Inputs are copied (mutating the vectors you
  passed in no longer changes the curve) and promoted to one concrete floating-point type per
  vector (`Int` → `Float64`; `Float32` + `BigFloat` → `BigFloat`; `Float64` +
  `ForwardDiff.Dual` → `Dual`; ranges and tuples are accepted). Construction throws an
  `ArgumentError` for a length mismatch, empty or non-finite inputs, negative, unsorted or
  duplicate tenors, or too few knots: `Spline.PCHIP()` and `Spline.Akima()` need 3; every other
  method accepts a single knot, which gives a flat curve. Direct construction, `reconstruct`,
  `fit`, and bootstrap raise the same errors.
- **Knot curves change only through `reconstruct`.** The knot vectors are read-only (indexed
  assignment, `.=`, `sort!`, and writes through `view` throw). `reconstruct` and `Accessors.@set`
  revalidate and rebuild every derived cache, so two curves can no longer compare `==` yet price
  differently. Setting a derived cache (`_f`, `_fn`, `_tail`, …) throws.
- **Equality is structural for every knot curve.** `==`, `isequal`, and `hash` compare the knot
  rates, knot tenors, interpolation method, and extrapolation policy, so independently built
  curves from equal inputs are equal and work as `Dict` keys. `isequal` keeps the signed-zero
  distinction of the underlying rates. Curves with different methods are unequal even when they
  are numerically identical (`Spline.Linear()` and `Spline.BSpline(1)`).
- **Negative times throw a `DomainError`** for every knot curve, as they already did for
  `ZeroRateCurve`.
- **`show` prints the construction call**, for example
  `ZeroRateCurve([0.02, 0.03], [1.0, 2.0], PolynomialSpline(1); extrapolation = :flat_forward)`,
  which rebuilds an equal curve.
- **Spline descriptors validate their order**: `Spline.PolynomialSpline(order)` accepts only
  orders 1, 2 and 3 (previously any order above 3 silently built a cubic spline), and
  `Spline.BSpline(d)` requires `d ≥ 1`.
- **The sampling form** `ZeroRateCurve(curve::AbstractYieldModel, tenors)` still sorts its tenor
  grid and requires `t > 0`. It now checks the method's minimum knot count before evaluating the
  source curve, and samples through `zero(curve, t)` instead of `-log(discount(curve, t))/t`,
  which is numerically stable at very small and very large tenors.
- **Fitting a knot curve is one rebuild per optimizer candidate.** `fit(curve, quotes)` works for
  any knot curve and varies all knot rates through a single batch `FinanceModels.KnotRatesOptic()`
  (previously one `@optic(_.rates[i])` per knot rebuilt the curve once per knot per candidate,
  O(n²) in the number of knots). Custom `variables` still work.
- **`fit` validates the knot grid before optimising**: loss fits and `Fit.Bootstrap()` throw the
  construction `ArgumentError` for duplicate or non-positive maturities, too few quotes for the
  interpolant, or `extrapolation = :extension` with `Spline.MonotoneConvex()`, before any solver
  runs. Optimizer trial curves are not validated (a non-finite candidate is a bad loss, not an
  exception). The returned curve is, so a fit that diverged to non-finite rates throws instead of
  returning a curve of `NaN`s.
- **Unsuccessful optimizer fits throw `FitConvergenceError`** carrying the solver's `retcode`,
  instead of returning the unfitted starting model. Catch it to retry with a different seed model
  or optimizer.

### Flat short end

!!! warning "Changed numbers: zero rates before the first knot"
    DataInterpolations-backed curves (`Spline.Linear()`, `Quadratic()`, `Cubic()`, `PCHIP()`,
    `Akima()`, `BSpline(d)`) now hold the zero rate flat at the first knot's rate between `t = 0`
    and the first knot, instead of extending the first interpolation piece back to `t = 0`. Values
    at and beyond the first knot do not change, and `Spline.MonotoneConvex()`, whose first interval
    from `t = 0` is part of the Hagan-West construction, does not change.

    For zero rates `[0.02, 0.025, 0.03, 0.035, 0.04]` at `[1, 2, 5, 10, 30]` with `Spline.Linear()`,
    the zero rate at 0.25 years moves from 1.625% (the 1-to-2-year line extended) to 2.0%, and at
    `t = 0` from 1.5% to 2.0%. A loss fit is affected only through cash flows before its earliest
    quote maturity, such as early coupons of a longer bond when the first knot is later than
    them.

- **`Fit.Bootstrap()` knots are exactly the quote maturities.** The returned curve used to carry
  an extra knot at `t = 0` tied to the first zero rate; the flat short end now gives the same
  curve, so bootstrapped values do not change, but `knot_tenors(curve)` has one entry per quote
  and `reconstruct(curve; rates = r)` takes one rate per quote.

### Long-end extrapolation

!!! warning "Changed numbers: zero rates beyond the last knot"
    Every DataInterpolations-backed curve and fit now extrapolates with `:flat_forward` by default instead of continuing its final interpolation piece. This changes zero rates, discount factors and present values **beyond the last knot** for `Yield.Spline`, `ZeroRateCurve` and `fit` (loss fits and `Fit.Bootstrap()`) with `Spline.Linear()`, `Quadratic()`, `Cubic()`, `PCHIP()`, `Akima()` and `BSpline(d)`. Values at and before the last knot do not change. `Spline.MonotoneConvex()`, the `ZeroRateCurve` default, keeps its v6.1 tail and does not change.

    For example, a linear bootstrap of `ZCBYield.([0.02, 0.025, 0.031, 0.036], [1, 2, 5, 10])` gives these continuously compounded zero rates:

    | Maturity | v6 (final piece continued) | v7 default | Change |
    |---------:|---------------------------:|-----------:|-------:|
    | 10 (last knot) | 3.5367% | 3.5367% | none |
    | 15 | 4.0205% | 3.6980% | −32 bp |
    | 20 | 4.5043% | 3.7786% | −73 bp |
    | 30 | 5.4719% | 3.8592% | −161 bp |

    At 30 years the discount factor moves from 0.1937 to 0.3142, so the present value of a 30-year cashflow rises by 62%. The v7 tail holds the last discrete forward, 4.0205% (the average forward from 5 to 10 years); the v6 linear continuation kept raising the zero rate by about 10 bp a year.

    **To keep the previous values**, pass `extrapolation = :extension` wherever the curve is built: `Yield.Spline(spline, tenors, rates; extrapolation = :extension)`, `ZeroRateCurve(rates, tenors, spline; extrapolation = :extension)`, or `fit(spline, quotes, method; extrapolation = :extension)`.

- **Knot curves take an `extrapolation` policy for the long end.** `ZeroRateCurve`, `Yield.Spline`, `Yield.MonotoneConvex`, `reconstruct` and the spline `fit` methods accept `extrapolation = :flat_forward` (default), `:flat_zero`, `:linear`, `:extension`, or `Yield.FlatForwardAt(rate)`. `:extension` restores the v6 continuation of the final DataInterpolations piece and is unavailable for `Spline.MonotoneConvex()`.
- **The `:flat_forward` anchor.** Beyond the last knot `(tₙ, zₙ)` the instantaneous forward is held at a constant `fₙ`, so `z(t) = fₙ + (zₙ - fₙ)tₙ/t` and discount factors are continuous at `tₙ`. DataInterpolations-backed curves use the last discrete forward `fₙ = (zₙtₙ - zₙ₋₁tₙ₋₁)/(tₙ - tₙ₋₁)` (for a single knot, `z₁`), so the forward can jump at `tₙ`; the tail no longer depends on the interpolant's end piece, whose endpoint forward can be extreme or negative (about −6.8% for `Spline.BSpline(3)` on an upward-sloping 2% to 4% curve). `Spline.MonotoneConvex()` keeps its boundary instantaneous forward, so its forward curve stays continuous. See [Interpolation Methods](interpolation.md).
- **`Yield.FlatForwardAt(rate)` supplies an independent terminal forward.** It requires a `FinanceCore.Rate` such as `Continuous(0.035)` or `Periodic(0.035, 1)`; a bare number throws an `ArgumentError`, because `Yield.Constant(0.035)` reads a bare number as annual effective. The rate is stored continuously compounded and stays fixed through fitting and Accessors updates. It preserves the final-knot discount factor and usually introduces a forward jump.
- **The policy is part of the curve.** It is preserved by `fit`, `reconstruct`, and `Accessors.@set`, compared by `==`, and readable as `curve.extrapolation`. The type parameters of `Yield.Spline` and `Yield.MonotoneConvex` are internal: dispatch on the type names or `Yield.AbstractInterpolatedZeroCurve`.
- **MonotoneConvex supports the tail policies natively.** Direct construction and loss-fitting `Spline.MonotoneConvex()` return a native `Yield.MonotoneConvex` for `:flat_forward`, `:flat_zero`, `:linear`, and `Yield.FlatForwardAt(rate)`; `Fit.Bootstrap()` still rejects the descriptor. At the final knot, `instantaneous_forward` reports the interior (left) value.

### Quote conventions

- **`OISYield` pays annually beyond one year.** Maturities over one year now build
  annual-pay par swaps, matching SOFR, €STR, and SONIA overnight index swaps; they
  were quarterly. Maturities of one year or less still settle once. Bootstrapped
  OIS curves change slightly at maturities over one year.
- **`ParSwapYield` requires `frequency`.** The quarterly default was removed because
  fixed-leg conventions differ by market. Write, for example,
  `ParSwapYield(r, t; frequency = 1)` for OIS-style annual fixed legs or
  `frequency = 2` for semiannual legs.
- **`InterestRateSwap` requires `frequency`** for the same reason; both legs use it.
  `InterestRateSwap(curve, 10)` becomes `InterestRateSwap(curve, 10; frequency = 4)`
  to keep the former quarterly legs.
- **`ParYield` rejects a conflicting `frequency`.** A `Periodic` rate sets its own
  frequency; passing a different `frequency` now throws an `ArgumentError` instead
  of being silently ignored. Convert the rate first, e.g. `Periodic(1)(r)`.

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
