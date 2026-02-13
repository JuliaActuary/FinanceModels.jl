"""
    ZeroRateCurve(rates, tenors, [spline])

A yield curve defined by continuously-compounded zero rates at specified tenors,
with interpolation between tenors via the existing `Spline` infrastructure.

The `spline` argument is a `Spline.SplineCurve` object (e.g. `Spline.Linear()`,
`Spline.Cubic()`). Defaults to `Spline.Linear()`.

The curve stores the raw `rates` vector, making it compatible with ForwardDiff:
construct `ZeroRateCurve(dual_rates, tenors, spline)` inside an AD closure and
the interpolation will propagate dual numbers.

# Examples

```julia
using FinanceModels

rates = [0.02, 0.03, 0.035, 0.04]
tenors = [1.0, 2.0, 5.0, 10.0]

zrc = ZeroRateCurve(rates, tenors)                    # linear interpolation
zrc_cubic = ZeroRateCurve(rates, tenors, Spline.Cubic())  # cubic

discount(zrc, 1.0)   # exp(-0.02 * 1.0)
discount(zrc, 3.5)   # interpolated rate at t=3.5
zero(zrc, 5.0)       # Continuous(0.035)
```
"""
struct ZeroRateCurve{R<:AbstractVector, T<:AbstractVector, S<:Sp.SplineCurve} <: AbstractYieldModel
    rates::R      # continuously-compounded zero rates
    tenors::T     # time points
    spline::S     # e.g., Spline.Linear(), Spline.Cubic()
end

ZeroRateCurve(rates, tenors) = ZeroRateCurve(rates, tenors, Sp.BSpline(1))

function FinanceCore.discount(zrc::ZeroRateCurve, t)
    if t <= zero(t)
        return one(eltype(zrc.rates))
    end
    curve = Yield.Spline(zrc.spline, zrc.tenors, zrc.rates)
    return discount(curve, t)
end

(zrc::ZeroRateCurve)(t) = FinanceCore.discount(zrc, t)
