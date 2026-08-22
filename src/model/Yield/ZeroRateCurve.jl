"""
    ZeroRateCurve(rates, tenors, [spline])
    ZeroRateCurve(curve::AbstractYieldModel, tenors; spline=Spline.MonotoneConvex())

A yield curve defined by continuously-compounded zero rates at specified tenors,
with interpolation between tenors via the existing `Spline` infrastructure.

The `spline` argument is a `Spline.SplineCurve` object (e.g. `Spline.MonotoneConvex()`,
`Spline.PCHIP()`, `Spline.Linear()`, `Spline.Cubic()`). Defaults to `Spline.MonotoneConvex()`.

The curve stores the `rates` vector with its original element type, making it compatible
with ForwardDiff: construct `ZeroRateCurve(dual_rates, tenors, spline)` inside an AD closure
and the interpolation will propagate dual numbers.

## Constructing from another yield model

The second form samples zero rates from any `AbstractYieldModel` (e.g. `Yield.Constant`,
`Yield.NelsonSiegel`, a fitted spline curve, etc.) at the specified `tenors`, producing
a `ZeroRateCurve` suitable for key rate analysis with ActuaryUtilities.jl. All tenors
must be positive (`t > 0`); they are sorted before sampling.

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

## Storage and validation

The curve **owns** its data: `rates` and `tenors` are copied at construction (so later
mutation of the vectors you passed in does not affect the curve) and promoted to a single
concrete floating-point element type (`Int` → `Float64`, `Float32` + `BigFloat` → `BigFloat`,
`Float64` + `ForwardDiff.Dual` → `Dual`). Ranges and tuples are accepted. The public
properties are `rates`, `tenors` and `spline`; the first two are exposed as **read-only**
vectors — indexed assignment (`zrc.rates[1] = …`, `.=`, `sort!`, …) throws. Use
`copy(zrc.rates)` for a mutable copy. The interpolation cache is internal and not accessible.

Construction throws an `ArgumentError` when:

- `rates` and `tenors` differ in length, or either is empty;
- any rate or tenor is not finite (`NaN`, `±Inf`);
- any tenor is negative, or the tenors are not strictly increasing (unsorted or duplicated);
- there are fewer knots than the interpolant needs: `Spline.PCHIP()` and `Spline.Akima()`
  need 3, `Spline.Linear()`/`Quadratic()`/`Cubic()`/`BSpline(n)` need 2,
  `Spline.MonotoneConvex()` needs 1.

A tenor of `0` is allowed in the direct form (you supply the instantaneous rate `r(0)`
explicitly); negative rates are allowed.

To change a curve, use `Accessors.@set`, which re-validates and rebuilds the interpolation:

```julia
using Accessors
zrc2 = @set zrc.rates[2] = 0.031          # one knot rate
zrc3 = @set zrc.tenors = [1.0, 3.0, 6.0, 12.0]  # whole grid (must stay strictly increasing)
zrc4 = @set zrc.spline = Spline.Linear()  # different interpolant over the same knots
```

`fit(zrc, quotes)` is supported, with one bounded optic per knot rate (as for
`Yield.MonotoneConvex`).

## Performance note

The interpolation model is built once at construction (`Yield.build_model(spline,
tenors, rates)`) and stored on the struct, so `discount` is a direct dispatch to
the prebuilt model rather than a rebuild per call. AD usage that creates a fresh
`ZeroRateCurve` per gradient step (the documented pattern) pays the build cost
once per step. The cache is a pure function of `(rates, tenors, spline)` and is ignored
by `==`, `isequal` and `hash`.

## Forward curve smoothness

The default `Spline.MonotoneConvex()` guarantees positive continuous forward rates
and produces C1-smooth forward curves ([Hagan & West, 2006](https://doi.org/10.1080/13504860600829233)).
For C2 smoothness, use `Spline.Cubic()`. `Spline.Linear()` produces kinks in the
forward curve at tenor points.
"""
struct ZeroRateCurve{R, T, S <: Sp.SplineCurve, M} <: AbstractYieldModel
    rates::ReadOnlyVector{R}   # continuously-compounded zero rates (finite); read-only
    tenors::ReadOnlyVector{T}  # finite, ≥ 0, strictly increasing; read-only
    spline::S                  # e.g., Spline.Linear(), Spline.MonotoneConvex()
    _model::M                  # pure function of (rates, tenors, spline); internal, not settable

    # The only constructor. Copies and promotes the inputs, validates them, builds the
    # interpolation cache, and wraps the storage read-only — so `_model` can never disagree
    # with `(rates, tenors, spline)`. There is deliberately no 4-arg constructor:
    # Accessors/ConstructionBase reconstruct through `setproperties`/`constructorof`
    # (see src/fit.jl), both of which route back here.
    function ZeroRateCurve(rates, tenors, spline::Sp.SplineCurve)
        r = __owned_float_vector(rates, "rates")
        t = __owned_float_vector(tenors, "tenors")
        length(r) == length(t) || throw(ArgumentError(
            "ZeroRateCurve: `rates` and `tenors` must have the same length (got $(length(r)) and $(length(t)))."))
        all(isfinite, r) || throw(ArgumentError("ZeroRateCurve: all rates must be finite (got $(r))."))
        all(isfinite, t) || throw(ArgumentError("ZeroRateCurve: all tenors must be finite (got $(t))."))
        first(t) >= zero(eltype(t)) || throw(ArgumentError("ZeroRateCurve: tenors must be ≥ 0 (got $(t))."))
        # finiteness is checked first so NaN cannot defeat the ordering comparison
        for i in 2:length(t)
            t[i] > t[i - 1] || throw(ArgumentError(
                "ZeroRateCurve: tenors must be strictly increasing (sorted, no duplicates); got $(t). " *
                "Sort the (rate, tenor) pairs before constructing."))
        end
        k = __min_knots(spline)
        length(t) >= k || throw(ArgumentError(
            "ZeroRateCurve: $(spline) requires at least $k knots (got $(length(t)))."))
        model = Yield.build_model(spline, t, r)   # receives the raw owned Vectors
        return new{eltype(r), eltype(t), typeof(spline), typeof(model)}(
            ReadOnlyVector(r), ReadOnlyVector(t), spline, model)
    end
end

ZeroRateCurve(rates, tenors) = ZeroRateCurve(rates, tenors, Sp.MonotoneConvex())

# Sampling form. The grid is normalised ONCE (a stateful iterator must not be consumed
# twice) and validated BEFORE the source curve is touched (so `Inf` never reaches
# `discount(curve, Inf)`), then sorted, checked for duplicates, and sampled.
function ZeroRateCurve(curve::AbstractYieldModel, tenors; spline = Sp.MonotoneConvex())
    t = sort!(__owned_float_vector(tenors, "tenors"))
    all(isfinite, t) || throw(ArgumentError("ZeroRateCurve: all tenors must be finite (got $(t))."))
    first(t) > zero(eltype(t)) || throw(ArgumentError(
        "All tenors must be positive (t > 0). The zero rate is undefined at t = 0."))
    # same adjacent strict comparison as the inner constructor (`allunique` would treat
    # Duals with equal primals but different partials as distinct)
    for i in 2:length(t)
        t[i] > t[i - 1] || throw(ArgumentError("ZeroRateCurve: tenors must be distinct (got $(t))."))
    end
    # Sample through the zero-rate interface rather than `-log(discount)/t`: that round-trip
    # is numerically unstable at extreme tenors (for a flat 5% curve it gives -0.0 at
    # t = 1e-20 and Inf at t = 2e4) and would trip the finite-rate validation on curves
    # that are mathematically fine.
    rates = [FinanceCore.rate(convert(Continuous(), Base.zero(curve, tᵢ))) for tᵢ in t]
    return ZeroRateCurve(rates, t, spline)
end

# Minimum knot counts that the backing interpolants can handle (verified against
# DataInterpolations 9; `Yield.Spline` already reduces a polynomial/B-spline order to n-1).
__min_knots(::Sp.SplineCurve) = 2     # Linear / Quadratic / Cubic / BSpline
__min_knots(::Sp.PCHIP) = 3
__min_knots(::Sp.Akima) = 3           # 2 knots reach an UndefRefError inside DataInterpolations
__min_knots(::Sp.MonotoneConvex) = 1  # single knot is a flat curve

# Owned, concretely-typed float `Vector` from any iterable of reals. Always copies (even
# when handed a `Vector`), so the curve never aliases caller-owned memory; never narrows
# AD types (`float(::Dual)` is a `Dual`).
function __owned_float_vector(x, what)
    v = x isa AbstractVector ? x : collect(x)
    isempty(v) && throw(ArgumentError("ZeroRateCurve: `$what` must not be empty."))
    T = isconcretetype(eltype(v)) ? eltype(v) : mapreduce(typeof, promote_type, v)
    T <: Real || throw(ArgumentError("ZeroRateCurve: `$what` must be real numbers (got element type $T)."))
    F = try
        float(T)
    catch
        T   # `float(::Type{Real})` etc. has no method; rejected just below
    end
    (F <: Real && isconcretetype(F)) || throw(ArgumentError(
        "ZeroRateCurve: `$what` must promote to one concrete floating-point type (got $F)."))
    return try
        Vector{F}(v)   # Array-from-AbstractArray always allocates a fresh copy
    catch e
        throw(ArgumentError("ZeroRateCurve: could not convert `$what` to Vector{$F}: $(sprint(showerror, e))"))
    end
end

# `_model` is internal: hidden from property access. Internal code uses `getfield`.
Base.propertynames(::ZeroRateCurve, ::Bool = false) = (:rates, :tenors, :spline)
function Base.getproperty(z::ZeroRateCurve, s::Symbol)
    s === :_model && throw(ArgumentError(
        "ZeroRateCurve: `_model` is internal (derived from rates, tenors and spline) and not accessible; use `discount`/`zero`."))
    return getfield(z, s)
end

function FinanceCore.discount(zrc::ZeroRateCurve, t)
    if iszero(t)
        return one(eltype(zrc.rates))
    end
    # negative times previously returned 1 silently; the curve has no defined
    # behavior before time zero, so error rather than misprice
    t < zero(t) && throw(DomainError(t, "ZeroRateCurve discount is only defined for t ≥ 0"))
    return discount(getfield(zrc, :_model), t)
end
# The callable `zrc(t) ≡ discount(zrc, t)` comes from the generic `AbstractYieldModel` fallback.

# Structural equality on the value-carrying fields. `_model` is a pure function of
# `(rates, tenors, spline)` — the storage is owned and read-only and the cache cannot be
# patched — so ignoring it is sound. `isequal` is defined fieldwise (not via `==`) so that
# `isequal(a, b)` implies `hash(a) == hash(b)` with Julia's own array semantics
# (`[-0.0] == [0.0]` but `!isequal([-0.0], [0.0])`).
Base.:(==)(a::ZeroRateCurve, b::ZeroRateCurve) =
    a.rates == b.rates && a.tenors == b.tenors && a.spline == b.spline
Base.isequal(a::ZeroRateCurve, b::ZeroRateCurve) =
    isequal(a.rates, b.rates) && isequal(a.tenors, b.tenors) && isequal(a.spline, b.spline)
Base.hash(z::ZeroRateCurve, h::UInt) =
    hash(z.rates, hash(z.tenors, hash(z.spline, h)))
