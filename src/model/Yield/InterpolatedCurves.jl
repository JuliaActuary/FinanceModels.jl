# ── Interpolated zero curves ───────────────────────────────────────────────────────────

"""
    Yield.AbstractInterpolatedZeroCurve <: Yield.AbstractYieldModel

Supertype of the yield curves that interpolate continuously compounded zero rates over a knot
grid: [`Yield.Spline`](@ref) (the DataInterpolations-backed methods) and
[`Yield.MonotoneConvex`](@ref). [`ZeroRateCurve`](@ref), [`reconstruct`](@ref), and spline
[`fit`](@ref FinanceModels.fit)s return one of these.

Every such curve:

- owns validated copies of its knots, read with [`knot_rates`](@ref) and
  [`knot_tenors`](@ref) (read-only vectors);
- changes only by building a new curve with [`reconstruct`](@ref) (or `Accessors.@set` on its
  `rates`, `tenors`, `spline`, or `extrapolation` properties, which calls `reconstruct`);
- follows its `extrapolation` policy beyond the last knot;
- throws a `DomainError` for negative times; and
- compares structurally: two curves are `==` (and `isequal`, with equal hashes) when their knot
  rates, knot tenors, interpolation method, and extrapolation policy are.

Before the first knot, the DataInterpolations-backed curves hold the zero rate flat at the first
knot's rate. `Yield.MonotoneConvex` instead interpolates from `t = 0` as part of the Hagan-West
construction. Other fields of the concrete types are internal.
"""
abstract type AbstractInterpolatedZeroCurve <: AbstractYieldModel end

"""
    knot_rates(curve::Yield.AbstractInterpolatedZeroCurve)

The continuously compounded zero rates at the curve's knots, in tenor order, as a read-only
vector. Use `copy` or `collect` for a mutable `Vector`.

See also [`knot_tenors`](@ref), [`reconstruct`](@ref).
"""
knot_rates(c::AbstractInterpolatedZeroCurve) = c.rates

"""
    knot_tenors(curve::Yield.AbstractInterpolatedZeroCurve)

The curve's knot tenors (strictly increasing, in years) as a read-only vector. Use `copy` or
`collect` for a mutable `Vector`.

See also [`knot_rates`](@ref), [`reconstruct`](@ref).
"""
knot_tenors(c::AbstractInterpolatedZeroCurve) = c.tenors

"""
    reconstruct(curve::Yield.AbstractInterpolatedZeroCurve;
        rates = knot_rates(curve), tenors = knot_tenors(curve),
        spline = <curve's method>, extrapolation = <curve's policy>)

Build a new curve from `curve` with any of its knot `rates`, knot `tenors`, interpolation
method (`spline`), or `extrapolation` policy replaced; omitted arguments keep the curve's own.
The result goes through the same validation as [`ZeroRateCurve`](@ref) and is built once.
Replacing `spline` can change the concrete type, for example to `Yield.MonotoneConvex`.

`reconstruct(curve)` is equal to `curve` and prices identically. `rates` may hold ForwardDiff
dual numbers, so this is how to differentiate a valuation with respect to a curve's knot rates:

```julia
curve = ZeroRateCurve([0.03, 0.035, 0.04], [1.0, 5.0, 10.0])
ForwardDiff.gradient(z -> pv(reconstruct(curve; rates = z), cfs), collect(knot_rates(curve)))
```

The derivatives are first order. At a kink of `Spline.MonotoneConvex()` (for example a flat
stretch of the curve), each partial is the centered response to a bump of that knot, and on a
flat stretch these need not add up to the parallel-shift derivative. `Spline.PCHIP()` and
`Spline.Akima()` throw at their kinks. See [Sensitivities Through Calibration](@ref).
"""
reconstruct(
    c::AbstractInterpolatedZeroCurve;
    rates = knot_rates(c), tenors = knot_tenors(c), spline = c.spline, extrapolation = c.extrapolation
) = ZeroRateCurve(rates, tenors, spline; extrapolation)

function __check_time(t, who)
    t < zero(t) && __throw_negative_time(t, who)
    return nothing
end
# Out of line: building the message inline makes every discount body too large to inline.
@noinline __throw_negative_time(t, who) = throw(DomainError(t, "$who is only defined for t ≥ 0"))

# Structural equality on the construction inputs. Every derived cache is rebuilt from these
# by the only constructors, and the stored knots are read-only. `isequal` is fieldwise (not via
# `==`) so that `isequal(a, b)` implies `hash(a) == hash(b)` with Julia's own array semantics
# (`[-0.0] == [0.0]` but `!isequal([-0.0], [0.0])`).
Base.:(==)(a::AbstractInterpolatedZeroCurve, b::AbstractInterpolatedZeroCurve) =
    a.spline == b.spline && a.extrapolation == b.extrapolation && a.tenors == b.tenors &&
    a.rates == b.rates
Base.isequal(a::AbstractInterpolatedZeroCurve, b::AbstractInterpolatedZeroCurve) =
    isequal(a.spline, b.spline) && isequal(a.extrapolation, b.extrapolation) &&
    isequal(a.tenors, b.tenors) && isequal(a.rates, b.rates)
Base.hash(c::AbstractInterpolatedZeroCurve, h::UInt) =
    hash(c.rates, hash(c.tenors, hash(c.extrapolation, hash(c.spline, hash(:AbstractInterpolatedZeroCurve, h)))))

# Base's convention: the derived caches are listed only with `propertynames(c, true)`.
Base.propertynames(c::AbstractInterpolatedZeroCurve, private::Bool = false) =
    private ? fieldnames(typeof(c)) : (:spline, :rates, :tenors, :extrapolation)

# ConstructionBase protocol, whose derived caches must follow their inputs:
#  * the properties are (spline, rates, tenors, extrapolation); `getproperties` excludes caches;
#  * `setproperties` rebuilds through `reconstruct`, whose keywords are exactly those properties:
#    a patch to a cache (or any unknown key) is a MethodError rather than silently ignored, and
#    `setproperties(c, getproperties(c)) == c`;
#  * `constructorof` discards the cache arguments and rebuilds, so
#    `constructorof(typeof(c))(getfields(c)...) == c`.
# Every `@set`/`setall`/`fit` reconstruction therefore revalidates and rebuilds the caches.
Accessors.ConstructionBase.getproperties(c::AbstractInterpolatedZeroCurve) =
    (spline = c.spline, rates = c.rates, tenors = c.tenors, extrapolation = c.extrapolation)
Accessors.ConstructionBase.setproperties(c::AbstractInterpolatedZeroCurve, patch::NamedTuple) =
    reconstruct(c; patch...)
Accessors.ConstructionBase.constructorof(::Type{<:AbstractInterpolatedZeroCurve}) =
    (spline, rates, tenors, extrapolation, _...) -> ZeroRateCurve(rates, tenors, spline; extrapolation)

# Print the construction call, which rebuilds an equal curve.
function Base.show(io::IO, c::AbstractInterpolatedZeroCurve)
    print(io, "ZeroRateCurve(")
    show(io, collect(c.rates))
    print(io, ", ")
    show(io, collect(c.tenors))
    print(io, ", ")
    show(io, c.spline)
    print(io, "; extrapolation = ")
    show(io, c.extrapolation)
    return print(io, ")")
end
