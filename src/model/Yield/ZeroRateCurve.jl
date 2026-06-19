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

The interpolation model is built once at construction (`Yield.build_model(spline,
tenors, rates)`) and stored on the struct, so `discount` is a direct dispatch to
the prebuilt model rather than a rebuild per call. AD usage that creates a fresh
`ZeroRateCurve` per gradient step (the documented pattern) pays the build cost
once per step — same total cost as before, just front-loaded from `discount`-time
to construction-time. The `_model` field is internal; do not rely on it.

## Forward curve smoothness

The default `Spline.MonotoneConvex()` guarantees positive continuous forward rates
and produces C1-smooth forward curves ([Hagan & West, 2006](https://doi.org/10.1080/13504860600829233)).
For C2 smoothness, use `Spline.Cubic()`. `Spline.Linear()` produces kinks in the
forward curve at tenor points.
"""
struct ZeroRateCurve{R<:AbstractVector, T<:AbstractVector, S<:Sp.SplineCurve, M} <: AbstractYieldModel
    rates::R      # continuously-compounded zero rates
    tenors::T     # time points
    spline::S     # e.g., Spline.Linear(), Spline.Cubic()
    _model::M     # prebuilt interpolation; do not access directly
end

# Three-arg constructor eagerly builds the interpolation. The 4-arg form
# (passing a prebuilt `_model`) is internal — used to thread an already-built
# model through wrapping constructors without rebuilding.
function ZeroRateCurve(rates, tenors, spline)
    model = Yield.build_model(spline, tenors, rates)
    return ZeroRateCurve(rates, tenors, spline, model)
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
    if iszero(t)
        return one(eltype(zrc.rates))
    end
    # negative times previously returned 1 silently; the curve has no defined
    # behavior before time zero, so error rather than misprice
    t < zero(t) && throw(DomainError(t, "ZeroRateCurve discount is only defined for t ≥ 0"))
    return discount(zrc._model, t)
end
# The callable `zrc(t) ≡ discount(zrc, t)` comes from the generic `AbstractYieldModel` fallback.

# Structural equality on the value-carrying fields. The `_model` field is a
# deterministic function of (rates, tenors, spline) and may not implement `==`
# on its underlying interpolation; ignoring it here keeps `a == b` meaningful
# whenever two curves are built from equal inputs.
Base.:(==)(a::ZeroRateCurve, b::ZeroRateCurve) =
    a.rates == b.rates && a.tenors == b.tenors && a.spline == b.spline
Base.hash(z::ZeroRateCurve, h::UInt) =
    hash(z.rates, hash(z.tenors, hash(z.spline, h)))
