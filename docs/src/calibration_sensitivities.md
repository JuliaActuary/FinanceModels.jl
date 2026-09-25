# Sensitivities Through Calibration

A fitted curve is a function of the market quotes it was fitted to. Two questions follow:

- how does a valuation change when the **quotes** change (differentiate through `fit`); and
- what quote does a curve **imply**, and how does it change when the curve changes
  ([`implied_quote`](@ref FinanceModels.Yield.implied_quote)).

Both are answered with [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), and both
give exact first-order derivatives: they come from the implicit function theorem at the solved
point, not from differentiating a solver's iterations. The values are always the primal
calculation, bitwise.

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
domestic curve, and the quotes may all carry dual numbers).

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
whose largest relative repricing residual exceeds `1e-6`. Refitting with bumped quotes and
taking finite differences is a much noisier check: optimizer noise of `1e-11` in the fitted rates
becomes an error of order `1e-4` in a difference quotient.

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
  Jacobian); or
- the model is not a spline fit (see the table).

Reverse-mode AD is not supported.

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
them, rebuild the curve with dual rates using [`reconstruct`](@ref FinanceModels.Yield.reconstruct):

```julia
z = collect(knot_rates(curve))
ForwardDiff.gradient(z -> pv(reconstruct(curve; rates = z), liability), z)
```
