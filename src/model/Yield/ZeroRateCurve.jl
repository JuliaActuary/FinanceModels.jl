"""
    ZeroRateCurve(rates, tenors, [spline])
    ZeroRateCurve(curve::AbstractYieldModel, tenors; spline=Spline.MonotoneConvex())

A yield curve defined by continuously-compounded zero rates at specified tenors,
with interpolation between tenors via the existing `Spline` infrastructure.

The `spline` argument is a `Spline.SplineCurve` object (e.g. `Spline.MonotoneConvex()`,
`Spline.PCHIP()`, `Spline.Linear()`, `Spline.Cubic()`). Defaults to `Spline.MonotoneConvex()`.

The curve stores the raw `rates` vector, making it compatible with ForwardDiff:
construct `ZeroRateCurve(dual_rates, tenors, spline)` inside an AD closure and
the interpolation will propagate dual numbers.

## Constructing from another yield model

The second form samples zero rates from any `AbstractYieldModel` (e.g. `Yield.Constant`,
`Yield.NelsonSiegel`, a fitted spline curve, etc.) at the specified `tenors`, producing
a `ZeroRateCurve` suitable for key rate analysis with ActuaryUtilities.jl. All tenors
must be positive (`t > 0`).

# Examples

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

# From a NelsonSiegel model:
ns = Yield.NelsonSiegel(1.0, 0.04, -0.02, 0.01)
zrc_ns = ZeroRateCurve(ns, [1.0, 2.0, 5.0, 10.0, 20.0])
```

## Performance note

`discount` reconstructs the interpolation model on each call for AD compatibility
(dual numbers must flow through the interpolation). For non-AD usage where `discount`
is called many times on the same curve, construct `Yield.build_model(zrc.spline,
zrc.tenors, zrc.rates)` once and use that instead. The AD pathway in ActuaryUtilities
builds the model once per gradient step, avoiding this per-call overhead.

## Forward curve smoothness

The default `Spline.MonotoneConvex()` guarantees positive continuous forward rates
and produces C1-smooth forward curves ([Hagan & West, 2006](https://doi.org/10.1080/13504860600829233)).
For C2 smoothness, use `Spline.Cubic()`. `Spline.Linear()` produces kinks in the
forward curve at tenor points.
"""
struct ZeroRateCurve{R<:AbstractVector, T<:AbstractVector, S<:Sp.SplineCurve} <: AbstractYieldModel
    rates::R      # continuously-compounded zero rates
    tenors::T     # time points
    spline::S     # e.g., Spline.Linear(), Spline.Cubic()
end

ZeroRateCurve(rates, tenors) = ZeroRateCurve(rates, tenors, Sp.MonotoneConvex())

function ZeroRateCurve(curve::AbstractYieldModel, tenors; spline=Sp.MonotoneConvex())
    all(t -> t > zero(t), tenors) || throw(ArgumentError(
        "All tenors must be positive (t > 0). The zero rate is undefined at t = 0."))
    tenors_f = collect(float.(tenors))
    rates = [-log(FinanceCore.discount(curve, t)) / t for t in tenors_f]
    perm = sortperm(tenors_f)
    return ZeroRateCurve(rates[perm], tenors_f[perm], spline)
end

function FinanceCore.discount(zrc::ZeroRateCurve, t)
    if t <= zero(t)
        return one(eltype(zrc.rates))
    end
    curve = Yield.build_model(zrc.spline, zrc.tenors, zrc.rates)
    return discount(curve, t)
end

(zrc::ZeroRateCurve)(t) = FinanceCore.discount(zrc, t)
