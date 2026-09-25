# Sensitivities Through Calibration

A fitted curve is a function of the market quotes it was fitted to. Two questions follow:

- how does a valuation change when the **quotes** change (differentiate through `fit`); and
- what quote does a curve **imply**, and how does it change when the curve changes
  ([`implied_quote`](@ref FinanceModels.Yield.implied_quote)).

Both are answered with [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), and both
give exact first-order derivatives: they come from the implicit function theorem at the solved
point, not from differentiating a solver's iterations. The values are always the primal
calculation, bitwise.

## The contract at a glance

- **First order only through solves.** Derivatives through `fit` and `implied_quote` are exact
  first derivatives. Nested dual numbers (a Hessian, or convexity through a calibration) throw.
  `reconstruct` solves nothing: a curve rebuilt from dual knot rates has derivatives of any
  order, except at an interpolation kink, where nested dual numbers throw.
- **With respect to what you differentiate.** Differentiating through `fit` gives risk to the
  quotes passed to it: the original market inputs. Quotes implied from a fitted curve in another
  family (for example `implied_quote` at the curve's knots) are a synthetic family. Risk to them
  equals risk to the market quotes only when the curve was fitted to that family at those tenors.
- **Exact fits only.** A differentiated fit must reprice its quotes. Bootstrap does; a loss fit
  is accepted only when one Newton correction of its knot rates is at most `1e-6` in every
  component (see "Accuracy" below).
- **Kinks.** Linear, quadratic, cubic, and B-spline interpolation are differentiable everywhere.
  At a kink of `Spline.MonotoneConvex()` (a flat stretch of the curve, or two equal adjacent
  forwards), each partial is the centered response to a bump in its own direction. **On a flat
  stretch these responses need not add up to the response to a parallel shift, so they do not
  aggregate like a gradient.** `Spline.PCHIP()`, `Spline.Akima()`, and fits that sit on a kink
  throw. See [Kinks](@ref calibration-kinks).

## Differentiating through `fit`

Pass dual numbers in the quotes and the fitted curve carries their derivatives:

```julia
using FinanceModels, ForwardDiff

tenors = [1.0, 2.0, 3.0, 5.0, 10.0]
rates = [0.03, 0.032, 0.035, 0.037, 0.04]
liability = Cashflow.([4.0, 4.0, 4.0, 4.0, 104.0], [1.0, 2.0, 3.0, 4.0, 5.0])

value(r) = pv(fit(Spline.Linear(), OISYield.(r, tenors), Fit.Bootstrap()), liability)

value(rates)                          # the valuation
ForwardDiff.gradient(value, rates)    # its sensitivity to each quoted OIS rate
```

The gradient is the change in value per unit change in each quoted rate, with every other quote
held fixed and the curve refitted. The same works for a loss fit with any interpolation method
(`fit(Spline.MonotoneConvex(), quotes)`, `fit(Spline.Cubic(), quotes)`, …), for pipelines that
build several curves, and for `FX.Forwards` fitted with a spline foreign curve (spot, the
domestic curve, and the quotes may all carry dual numbers). The shape-preserving interpolations
have kinks where no derivative exists; see [Kinks](@ref calibration-kinks) below.

### How it works

A spline fit places one knot at each quote maturity, so the fitted knot rates `z` solve the
square repricing system

```math
R_i(z, p) = \operatorname{pv}(\text{curve}(z), q_i) - \text{price}_i = 0 .
```

`fit` solves it on primal copies of the quotes, then attaches
``\mathrm{d}z = -(\partial R/\partial z)^{-1} (\partial R/\partial p)\,\mathrm{d}p`` to the
knot rates. The Jacobian ``\partial R/\partial z`` is computed on primal values, and your dual
numbers enter only through one evaluation of ``R`` with ordinary arithmetic, so FinanceModels'
internal derivatives never mix with yours.

### What is supported

| Calibration | Differentiable | Notes |
|:------------|:---------------|:------|
| `fit(spline, quotes)` (`Fit.Loss`), every interpolation method | yes | the fit must reprice its quotes; see below |
| `fit(Spline.Linear(), quotes, Fit.Bootstrap())` | yes | exact to root-finder precision |
| `fit(FX.Forwards(pair, spot, domestic, spline), quotes, …)` | yes | through the implied foreign quotes |
| `fit(Yield.SmithWilson(…), quotes)` | quote prices only | closed form; dual coupon amounts are not yet supported |
| other models (`Yield.NelsonSiegel`, `Yield.Constant`, short-rate models, refitting an existing curve) | no | `fit` throws an `ArgumentError` |

Dual numbers may appear in quote prices, in the rates of `Bond.Fixed` and `Bond.Floating`, in
the amounts of `Cashflow`, `Composite`, and `FX.BasisSwapLeg` instruments (which covers
`ZCBPrice`, `ZCBYield`, `ParYield`, `CMTYield`, `OISYield`, `ParSwapYield`, and FX quotes), and
in a `Yield.FlatForwardAt` extrapolation forward.

### Accuracy

The derivative is exact for the exactly repricing curve. Bootstrap reprices to root-finder
precision. A loss fit stops at its optimizer's tolerance, and `fit` refuses to differentiate one
unless the largest absolute component of one Newton correction of its knot rates towards the
exact fit is at most `1e-6`. The correction is a local estimate of the fit's error, not a
guaranteed distance to the exact solution; like the conditioning check, it does not depend on
the quotes' notionals. To tighten a loss fit, pass solver settings through `fit`'s
`solve_kwargs`, for example `fit(Spline.MonotoneConvex(), quotes; solve_kwargs = (; g_tol = 1e-12))`.
Refitting with bumped quotes and taking finite differences is a much noisier check: optimizer
noise of `1e-11` in the fitted rates becomes an error of order `1e-4` in a difference quotient.

### Errors

`fit` throws an `ArgumentError`, rather than return a curve with missing or wrong derivatives,
when:

- a quote maturity or cashflow time carries a dual number (derivatives with respect to maturities
  are not supported);
- the dual numbers are nested (second-order derivatives) or come from two different ForwardDiff
  calls;
- a dual number sits inside a contract type `fit` does not know how to strip;
- a loss fit does not reprice its quotes (for example with a loss function whose minimum is not
  at zero residual);
- the quote prices do not determine the knot rates (a singular or ill-conditioned repricing
  Jacobian);
- the fitted curve lies on, or within the fit's precision of, a kink of its interpolation (for
  example, flat quotes fitted with `Spline.MonotoneConvex()`, `Spline.PCHIP()`, or
  `Spline.Akima()`); or
- the model is not a spline fit (see the table).

Reverse-mode AD is not supported. Derivatives are first order only: second derivatives through
a fit (convexity with respect to the quotes) throw the nested-dual error above.

## [Kinks in shape-preserving interpolations](@id calibration-kinks)

`Spline.Linear()`, `Spline.Cubic()`, `Spline.Quadratic()`, and `Spline.BSpline(n)` are linear in
their knot rates, so their derivatives exist everywhere. `Spline.MonotoneConvex()`,
`Spline.PCHIP()`, and `Spline.Akima()` preserve shape by switching formula as the knots change,
so they are only piecewise smooth in their knot rates. The switches are special configurations
of the knot *values*, not of the tenors:

| Interpolation | Kinks |
|:--------------|:------|
| `Spline.MonotoneConvex()` | two adjacent discrete forwards equal (every flat stretch of the curve); a node forward exactly at its positivity bound |
| `Spline.PCHIP()` | two adjacent knot rates equal (a flat segment) |
| `Spline.Akima()` | three or more knots on a straight line; a knot whose slope weight is exactly `1e-9` of the largest, where the curve's value jumps |

Flat curves sit on a kink, and repeated or rounded quotes can put a fitted curve on one.
Away from a kink, every derivative on this page is exact. Close to one it is still exact, but a
bump of one basis point may cross the kink and move the value differently.

At a kink, the derivative depends on the direction of the bump:

- **`Spline.MonotoneConvex()`** reports, for each partial, the limit of a centered bump in that
  partial's direction (the average of the up and down derivatives). Where two adjacent forwards
  are equal on an otherwise sloped curve, this limit is linear in the direction, so single-knot
  sensitivities add up to the parallel one. On a flat stretch it is not: for a flat curve, the
  sensitivities to each knot need not sum to the sensitivity to a parallel shift. A parallel
  shift of a flat curve keeps it flat, and its derivative is exact. Second derivatives do not
  exist at a kink and throw, as does a node forward at its bound when a discrete forward is
  exactly zero (a 0% curve).
- **`Spline.PCHIP()`** and **`Spline.Akima()`** throw an `ArgumentError` when dual knot rates move
  a kink, since their implementation returns a one-sided value, `NaN`, or the derivative of a
  fallback formula there, and Akima's value can jump. A parallel shift of a flat PCHIP or Akima curve is fine.
- **`fit`** throws when the fitted curve lies on a kink, or within the fit's precision of one:
  the refitted curve has no derivative with respect to the quotes there.

For key-rate risk on any of these curves, shifts added on top of the curve (such as
ActuaryUtilities' `KeyRates`) are smooth: they do not re-interpolate the knots.

## Implied quotes

[`implied_quote`](@ref FinanceModels.Yield.implied_quote) inverts a quote constructor on a curve:
it returns the rate (or price) at which the quote reprices, in the constructor's own convention.

```julia
curve = fit(Spline.Linear(), OISYield.(rates, tenors), Fit.Bootstrap())

implied_quote(curve, OISYield, 7.0)                                   # 7-year OIS rate
implied_quote(curve, (r, t) -> ParYield(r, t; frequency = 2), 7.0)   # semiannual par yield
implied_quote.(curve, OISYield, tenors) ≈ rates                       # the curve's own quotes
```

Its derivatives with respect to the curve are exact, so for example a curve bumped with a dual
spread gives the sensitivity of a quote to that spread:

```julia
bumped(s) = curve + Yield.Constant(Continuous(s))
ForwardDiff.derivative(s -> implied_quote(bumped(s), OISYield, 7.0), 0.0)
```

`implied_quote` agrees with [`par`](@ref FinanceModels.Yield.par) for par yields:
`rate(par(curve, t; frequency = f)) ≈ implied_quote(curve, (r, T) -> ParYield(r, T; frequency = f), t)`.

## Knot-rate sensitivities without refitting

To differentiate with respect to a fitted curve's own knot rates, rather than the quotes behind
them, rebuild the curve with dual rates using [`reconstruct`](@ref FinanceModels.Yield.reconstruct)
(see [Kinks](@ref calibration-kinks) for the shape-preserving interpolations):

```julia
z = collect(knot_rates(curve))
ForwardDiff.gradient(z -> pv(reconstruct(curve; rates = z), liability), z)
```
