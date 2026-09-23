# FinanceModels.jl release notes

## v7.0.0 (unreleased)

### `implied_quote`

`implied_quote(curve, family, maturity)` returns the quote at which a quote
constructor (`CMTYield`, `OISYield`, `ZCBYield`, `ZCBPrice`, or a closure such as
`(r, t) -> ParYield(r, t; frequency = 1)`) reprices on a curve. Its first-order
ForwardDiff derivatives with respect to curve parameters are exact: they come
from the implicit function theorem at the solution rather than from solver
iterations. FinanceModels now depends on ForwardDiff directly.

### Bootstrap requires linear interpolation

`Fit.Bootstrap()` now accepts only `Spline.Linear()` and throws an `ArgumentError`
for other strategies. Smoother interpolants let later knots reshape earlier segments,
so bootstrapped coupon quotes silently stopped repricing (residuals up to 0.40% of
par). Fit smoother curves to all quotes at once with `Fit.Loss`. Bootstrap also
validates its inputs and reprices every quote on the returned curve. See the
migration guide.

### `ZeroRateCurve` returns the curve it builds; one knot-curve interface

`ZeroRateCurve` is now a function returning a `Yield.MonotoneConvex` (for the default
`Spline.MonotoneConvex()`) or a `Yield.Spline`. Both subtype the new
`Yield.AbstractInterpolatedZeroCurve`, as do spline `fit` and bootstrap results. Read knots with
`knot_rates`/`knot_tenors` and build a changed curve with `reconstruct(curve; rates, tenors,
spline, extrapolation)`, which also accepts dual numbers for knot-rate gradients. Knot curves copy
and validate their inputs through one shared construction (unsorted, duplicate, negative or
non-finite tenors, non-finite rates, and too few knots throw an `ArgumentError`), keep read-only
knot vectors, compare structurally (`==`/`isequal`/`hash` on knots, method and extrapolation
policy), throw a `DomainError` for negative times, and print as the `ZeroRateCurve(...)` call that
rebuilds them. Fitting a knot curve rebuilds it once per optimizer candidate instead of once per
knot. `Yield.build_model` and the `Yield.MonotoneConvex()` fit placeholder are removed:
`Spline.MonotoneConvex()` is the one monotone convex selector (#272), fitted through the same path
as every other interpolation method. `Spline.PolynomialSpline` accepts only orders 1 to 3 and
`Spline.BSpline` requires degree 1 or more.

### Flat zero rate before the first knot (changed numbers)

DataInterpolations-backed curves now hold the first knot's zero rate between `t = 0` and the
first knot instead of extending the first interpolation piece. **Zero rates before the first knot
change**; values from the first knot on do not, and `Spline.MonotoneConvex()` does not change.
Bootstrapped curves are numerically unchanged, but their knots are now exactly the quote
maturities (the extra `t = 0` knot is gone), so a curve has one knot rate per quote.

### Failed optimizer fits throw `FitConvergenceError`

Every optimizer-backed `fit` now checks the solver's return code and throws a
`FitConvergenceError` (with the solver's `retcode`) instead of returning the
unfitted starting model.

### Flat-forward extrapolation beyond the last knot (changed numbers)

DataInterpolations-backed curves (`Spline.Linear()`, `Quadratic()`, `Cubic()`,
`PCHIP()`, `Akima()`, `BSpline(d)`) no longer continue their final interpolation
piece beyond the last knot. By default they now hold the last discrete forward,
`(zₙtₙ - zₙ₋₁tₙ₋₁)/(tₙ - tₙ₋₁)`, constant there, which keeps discount factors
continuous at the last knot and stops the end polynomial from driving the tail
(continuing it, or even holding its end slope's forward, can give extreme or
negative long-end forwards).
**Zero rates and present values beyond the last knot change** for every such
curve and fit; values up to the last knot do not. For a linear bootstrap of
`ZCBYield.([0.02, 0.025, 0.031, 0.036], [1, 2, 5, 10])` the 30-year continuous
zero rate moves from 5.47% to 3.86% (−161 bp; the 30-year discount factor rises
from 0.194 to 0.314). `Spline.MonotoneConvex()`, the `ZeroRateCurve` default,
keeps its boundary forward and is unchanged.

Pass `extrapolation = :extension` to `Yield.Spline`, `ZeroRateCurve` or `fit` to
keep the previous values. The other policies are
`:flat_zero`, `:linear` and `Yield.FlatForwardAt(rate)`, which takes a
`FinanceCore.Rate` such as `Continuous(0.035)` (a bare number throws). See the
migration guide.

### Quote conventions

- `OISYield` builds annual-pay par swaps beyond one year (was quarterly), matching
  SOFR, €STR, and SONIA overnight index swaps.
- `ParSwapYield` and `InterestRateSwap` require an explicit `frequency`.
- `ParYield` throws an `ArgumentError` when an explicit `frequency` conflicts with a
  `Periodic` rate's own compounding; previously the keyword was ignored.

See the migration guide for upgrade steps.

## v6.4.0

### `Spline.BSpline` fitted values changed on non-uniform tenor grids

The only change in this release is the move from DataInterpolations 8 to 9 (#293).
That bump silently changed what `Spline.BSpline(d)` produces whenever the knot
tenors are not equally spaced, which includes every standard CMT/Treasury grid.
Fits on equally spaced tenors are bit-identical to v6.3.0.

**Mechanism.** DataInterpolations 8 evaluated `BSplineInterpolation` as
`S(p(t))`: the maturity `t` was first mapped piecewise-linearly onto an equally
spaced parameter grid (one parameter value per knot), and the degree-`d` spline
was built in that parameter space. FinanceModels passed the `:Uniform`
parameter-vector option. In effect, `BSpline(d)` was a degree-`d` spline in
*knot index*, warped back to maturity. Where two neighbouring tenor intervals
differ in width the warp has a slope break, so the zero curve had a derivative
kink and the instantaneous forward rate jumped at that knot. On the CMT grid
below the v6.3.0 `BSpline(4)` forward curve jumps by roughly −63 bps at the 10y
knot; `BSpline(3)` by −71 bps.

DataInterpolations 9 removed the parameter vector. The spline is now built and
evaluated directly in maturity with de Boor averaged knots, so `BSpline(d)` is
what its name says: a C^(d−1) piecewise polynomial in `t`. Forwards are
continuous at every knot. The old curve cannot be reproduced; the concept no
longer exists upstream.

**Consequence.** A true degree-4 interpolant across tenors spanning 1 month to
30 years oscillates between the long knots. With the representative CMT quotes

```julia
tenors = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]
quotes = [0.0530,0.0528,0.0525,0.0515,0.0490,0.0450,0.0430,0.0420,0.0425,0.0435,0.0465,0.0455]
c = fit(Spline.BSpline(4), CMTYield.(quotes, tenors), Fit.Bootstrap())
```

the fitted curve reprices its own 20y quote with an error of −20.8 bps
(v6.3.0: +1.0 bps) and the instantaneous forward at the 20y knot is negative.
`BSpline(3)` on the same grid improves, from 1.0 bps to 0.2 bps at 20y.
`BSpline(2)`, `Linear`, `Cubic`, `PCHIP`, `Akima` and `MonotoneConvex` are
unaffected or slightly better.

**What to do.**

- If you pinned curve values from a `BSpline` fit on unequally spaced tenors,
  re-baseline; the change is expected.
- Prefer `Spline.BSpline(3)` (or `Cubic`, `PCHIP`, `MonotoneConvex`) for
  Treasury-style grids. Degree 4 and above on strongly non-uniform tenors will
  oscillate; that is inherent to a global interpolating spline and not a bug.
