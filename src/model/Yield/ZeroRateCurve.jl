"""
    ZeroRateCurve(rates, tenors, spline=Spline.MonotoneConvex(); extrapolation=:flat_forward)
    ZeroRateCurve(rates, tenors; spline=Spline.MonotoneConvex(), extrapolation=:flat_forward)
    ZeroRateCurve(curve::AbstractYieldModel, tenors;
        spline=Spline.MonotoneConvex(), extrapolation=:flat_forward)

Build a yield curve that interpolates continuously-compounded zero `rates` at `tenors` with the
interpolation method `spline`: `Spline.MonotoneConvex()` (the default), `Spline.PCHIP()`,
`Spline.Akima()`, `Spline.Linear()`, `Spline.Quadratic()`, `Spline.Cubic()`, or
`Spline.BSpline(n)`.

The result is a [`Yield.AbstractInterpolatedZeroCurve`](@ref): a [`Yield.MonotoneConvex`](@ref)
for `Spline.MonotoneConvex()`, and a [`Yield.Spline`](@ref) otherwise. `ZeroRateCurve` is a
construction function, not a type. Dispatch on `Yield.AbstractInterpolatedZeroCurve` when code
needs a curve's knots, or on `Yield.AbstractYieldModel` when it only discounts. Read the knots
with [`knot_rates`](@ref) and [`knot_tenors`](@ref), and build a changed curve with
[`reconstruct`](@ref).

Before the first tenor, the DataInterpolations-backed methods hold the zero rate flat at the
first rate; `Spline.MonotoneConvex()` interpolates from `t = 0` as part of its construction. The
default `extrapolation=:flat_forward` holds the instantaneous forward rate constant beyond the
last tenor while preserving the last-tenor discount factor. With `Spline.MonotoneConvex()` the
constant is its boundary instantaneous forward, so the forward curve stays continuous. With the
DataInterpolations-backed splines it is the last discrete forward,
`(zₙtₙ - zₙ₋₁tₙ₋₁)/(tₙ - tₙ₋₁)`, and the forward can jump at the last tenor (see
[`Yield.Spline`](@ref FinanceModels.Yield.Spline)). `:flat_zero` instead holds the last zero
rate; `:linear` extends the zero rate at its boundary slope; and `:extension` continues the
final DataInterpolations polynomial piece (unavailable with `Spline.MonotoneConvex()`).
[`FlatForwardAt`](@ref) holds a user-specified forward rate fixed, including when knot rates are
fitted or changed.

## Constructing from another yield model

The third form samples zero rates from any `AbstractYieldModel` (e.g. `Yield.Constant`,
`Yield.NelsonSiegel`, a fitted curve) at the given `tenors`. All tenors must be positive
(`t > 0`); they are sorted before sampling, and the grid is validated before the source curve is
evaluated.

# Examples

```julia
using FinanceModels

rates = [0.02, 0.03, 0.035, 0.04]
tenors = [1.0, 2.0, 5.0, 10.0]

zrc = ZeroRateCurve(rates, tenors)                              # Yield.MonotoneConvex
zrc_pchip = ZeroRateCurve(rates, tenors, Spline.PCHIP())        # Yield.Spline (PCHIP)
zrc_lin = ZeroRateCurve(rates, tenors, Spline.Linear())         # Yield.Spline (linear)
zrc_flat_zero = ZeroRateCurve(rates, tenors, Spline.Cubic(); extrapolation=:flat_zero)

discount(zrc, 1.0)   # exp(-0.02 * 1.0)
zero(zrc, 5.0)       # Continuous(0.035)
knot_rates(zrc)      # the rates, read-only
reconstruct(zrc; rates = rates .+ 0.001)   # a new curve, 10bp higher at every knot

# From a NelsonSiegel model:
ns = Yield.NelsonSiegel(1.0, 0.04, -0.02, 0.01)
zrc_ns = ZeroRateCurve(ns, [1.0, 2.0, 5.0, 10.0, 20.0])
```

## Validation

`rates` and `tenors` are copied (later mutation of the vectors you passed in does not affect the
curve) and each promoted to a single concrete floating-point element type (`Int` → `Float64`,
`Float32` + `BigFloat` → `BigFloat`, `Float64` + `ForwardDiff.Dual` → `Dual`), so
`ZeroRateCurve(dual_rates, tenors, spline)` inside an AD closure propagates derivatives. Ranges
and tuples are accepted. Construction throws an `ArgumentError` when:

- `rates` and `tenors` differ in length, or either is empty;
- any rate or tenor is not finite (`NaN`, `±Inf`);
- any tenor is negative, or the tenors are not strictly increasing (unsorted or duplicated);
- `extrapolation` is neither `FlatForwardAt(forward)` nor one of `:flat_forward`,
  `:flat_zero`, `:linear`, or `:extension`, or `:extension` is requested with
  `Spline.MonotoneConvex()`;
- there are fewer knots than the interpolant needs: `Spline.PCHIP()` and `Spline.Akima()`
  need 3; the other methods accept a single knot, which gives a flat curve.

A tenor of `0` is allowed in the direct form (you supply the instantaneous rate `r(0)`
explicitly); negative rates are allowed.

## Forward curve smoothness

With the default extrapolation and inputs implying positive discrete forwards,
`Spline.MonotoneConvex()` guarantees positive continuous forward rates ([Hagan & West, 2006](https://doi.org/10.1080/13504860600829233)).
For C2 zero-rate smoothness between knots, use `Spline.Cubic()`. `Spline.Linear()` produces
kinks in the forward curve at tenor points. At the last tenor, the default `:flat_forward`
policy keeps the instantaneous forward continuous for `Spline.MonotoneConvex()`; for the other
interpolants, and under the other policies, the forward can jump there.
"""
ZeroRateCurve(rates, tenors, spline::Sp.SplineCurve; extrapolation = :flat_forward) =
    __build_public(spline, KnotGrid(rates, tenors, spline; who = "ZeroRateCurve"); extrapolation)

ZeroRateCurve(rates, tenors; spline::Sp.SplineCurve = Sp.MonotoneConvex(), extrapolation = :flat_forward) =
    ZeroRateCurve(rates, tenors, spline; extrapolation)

# Sampling form. The grid is normalised ONCE (a stateful iterator must not be consumed
# twice) and validated BEFORE the source curve is touched (so `Inf` never reaches
# `discount(curve, Inf)`), then sorted, checked for duplicates and the method's minimum knot
# count, and sampled.
function ZeroRateCurve(
        curve::AbstractYieldModel, tenors;
        spline::Sp.SplineCurve = Sp.MonotoneConvex(), extrapolation = :flat_forward
    )
    t = sort!(__owned_float_vector(tenors, "tenors", "ZeroRateCurve"))
    all(isfinite, t) || throw(ArgumentError("ZeroRateCurve: all tenors must be finite (got $(t))."))
    first(t) > zero(eltype(t)) || throw(
        ArgumentError(
            "All tenors must be positive (t > 0). The zero rate is undefined at t = 0."
        )
    )
    # same adjacent strict comparison as `KnotGrid` (`allunique` would treat Duals with equal
    # primals but different partials as distinct)
    for i in 2:length(t)
        t[i] > t[i - 1] || throw(ArgumentError("ZeroRateCurve: tenors must be distinct (got $(t))."))
    end
    __check_min_knots(spline, length(t), "ZeroRateCurve")
    # Sample through the zero-rate interface rather than `-log(discount)/t`: that round-trip
    # is numerically unstable at extreme tenors (for a flat 5% curve it gives -0.0 at
    # t = 1e-20 and Inf at t = 2e4) and would trip the finite-rate validation on curves
    # that are mathematically fine.
    rates = [FinanceCore.rate(convert(Continuous(), Base.zero(curve, tᵢ))) for tᵢ in t]
    return ZeroRateCurve(rates, t, spline; extrapolation)
end
