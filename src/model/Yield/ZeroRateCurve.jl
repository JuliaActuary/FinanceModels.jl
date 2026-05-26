"""
    ZeroRateCurve(rates, tenors, [spline])
    ZeroRateCurve(curve::AbstractYieldModel, tenors; spline=Spline.MonotoneConvex())

Construct an `AbstractYieldModel` from continuously-compounded zero rates at
the given tenors, interpolated by `spline` (default: `Spline.MonotoneConvex()`).

`ZeroRateCurve` is a factory function that returns the appropriate underlying
curve type — `Yield.MonotoneConvex` for the default spline, `Yield.Spline` for
all other splines. The returned object supports the full `AbstractYieldModel`
interface (`discount`, `zero`, `forward`, callable `(t)` syntax, `pv`, etc.),
so downstream code that programs against the abstract interface continues to
work without modification.

## Constructing from rates + tenors

```julia
using FinanceModels

rates = [0.02, 0.03, 0.035, 0.04]
tenors = [1.0, 2.0, 5.0, 10.0]

zrc = ZeroRateCurve(rates, tenors)                              # default: MonotoneConvex
zrc_pchip = ZeroRateCurve(rates, tenors, Spline.PCHIP())        # PCHIP
zrc_lin = ZeroRateCurve(rates, tenors, Spline.Linear())          # linear
zrc_cubic = ZeroRateCurve(rates, tenors, Spline.Cubic())         # cubic

discount(zrc, 1.0)   # exp(-0.02 * 1.0)
discount(zrc, 3.5)   # interpolated rate at t=3.5
zero(zrc, 5.0)       # Continuous(0.035)
```

## Resampling another yield model

The second form samples zero rates from any existing `AbstractYieldModel`
(e.g. `Yield.Constant`, `Yield.NelsonSiegel`, a fitted spline curve, etc.) at
the specified `tenors`, then re-interpolates with the chosen spline:

```julia
ns = Yield.NelsonSiegel(1.0, 0.04, -0.02, 0.01)
zrc_ns = ZeroRateCurve(ns, [1.0, 2.0, 5.0, 10.0, 20.0])
```

All tenors must be positive (`t > 0`); the zero rate is undefined at `t = 0`.

## Why a factory function

Earlier versions of FinanceModels exposed `ZeroRateCurve` as a concrete struct
holding `(rates, tenors, spline)` with a per-call lazy rebuild of the
interpolation model. That design existed to support an AD pathway that tagged
`zrc.rates` with `ForwardDiff` duals; that pathway has been superseded by
`Yield.TenorShift` bumps in ActuaryUtilities. The factory function form avoids
the per-call rebuild and removes a redundant concrete type from the curve
hierarchy.

If you need access to the spline-resampled rates+tenors data that the old
struct stored, compute them once at construction and keep them as local
variables alongside the curve.

## Forward curve smoothness

The default `Spline.MonotoneConvex()` guarantees positive continuous forward
rates and produces C1-smooth forward curves ([Hagan & West, 2006](https://doi.org/10.1080/13504860600829233)).
For C2 smoothness, use `Spline.Cubic()`. `Spline.Linear()` produces kinks in
the forward curve at tenor points.
"""
ZeroRateCurve(rates, tenors, spline = Sp.MonotoneConvex()) =
    Yield.build_model(spline, tenors, rates)

function ZeroRateCurve(curve::AbstractYieldModel, tenors; spline=Sp.MonotoneConvex())
    all(t -> t > zero(t), tenors) || throw(ArgumentError(
        "All tenors must be positive (t > 0). The zero rate is undefined at t = 0."))
    tenors_f = collect(float.(tenors))
    rates = [-log(FinanceCore.discount(curve, t)) / t for t in tenors_f]
    perm = sortperm(tenors_f)
    return Yield.build_model(spline, tenors_f[perm], rates[perm])
end
