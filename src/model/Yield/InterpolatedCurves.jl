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

# Every knot curve's time functions share one skeleton: a negative time throws; beyond the last knot
# the curve's tail, `__tail(c)` (a `CurveTail`); otherwise its interior, `__interior_zero`,
# `__interior_log_discount` and `__interior_forward`. A curve that continues its interpolant past the
# last knot (`:extension`, `__extends(c)`) has no tail. The zero rate and L are continuous at the last
# knot, so the tail starts after it; the instantaneous forward is the right-hand derivative of L, the
# tail's from the last knot on.
@inline __in_tail(c::AbstractInterpolatedZeroCurve, t) = !__extends(c) && t > __tail(c).last_tenor
function Base.zero(c::AbstractInterpolatedZeroCurve, t)
    __check_time(t, "zero")
    return Continuous(__in_tail(c, t) ? __tail(c)(t) : __interior_zero(c, t))
end
@inline __knot_log_discount(c::AbstractInterpolatedZeroCurve, t) =
    __in_tail(c, t) ? __tail(c)(t) * t : __interior_log_discount(c, t)

"""
    instantaneous_forward(curve::AbstractInterpolatedZeroCurve, t)

The instantaneous (continuously compounded) forward rate of a knot curve at `t`: the right-hand
derivative of `-log(discount(curve, t))`. At a knot where the interpolant's pieces meet it is the
forward of the piece that starts there (`Spline.MonotoneConvex()`'s forward is continuous at its
interior knots). From the last knot on it is the tail's forward, which follows `curve.extrapolation`:
the default `:flat_forward` holds it constant, `:flat_zero` holds it at the last zero rate,
`:linear` derives it from the linearly extended zero rate, and `FlatForwardAt(f)` holds it at `f`;
at `t = Inf` it is the tail's limit. Under `:extension` the last piece continues, and `t = Inf` is
unsupported.

Note this is distinct from `forward(curve, from, to)`, which is the *discrete* forward `Rate`
between two times and is defined for every yield model.
"""
function instantaneous_forward(c::AbstractInterpolatedZeroCurve, t)
    __check_time(t, "instantaneous_forward")
    (__extends(c) || t < __tail(c).last_tenor) || return __tail_forward(__tail(c), t)
    return __interior_forward(c, t)
end

function __log_tail(c::AbstractInterpolatedZeroCurve)
    __extends(c) && __throw_extension_tail()
    return __log_tail(__tail(c))
end
@noinline __throw_extension_tail() = throw(
    DomainError(
        Inf, "discount(curve, Inf) is unsupported with extrapolation = :extension, whose " *
            "polynomial continuation of the last piece is evaluated at finite times only. Use a " *
            "finite time, or a policy with a tail limit (:flat_forward, :flat_zero, :linear, FlatForwardAt)."
    )
)

# Every knot curve evaluates alike: a negative time throws, t = Inf is the limit of the curve's tail
# (`__log_tail`, out of line so that the finite path stays small enough to inline), and a finite
# time is the curve's own L, `__knot_log_discount(c, t)`.
function FinanceCore.discount(c::AbstractInterpolatedZeroCurve, t)
    __check_time(t, "discount")
    isinf(t) && return __discount_at_infinity(c)
    return exp(-__knot_log_discount(c, t))
end
function __log_discount(c::AbstractInterpolatedZeroCurve, t)
    __check_time(t, "discount")
    isinf(t) && return __log_discount_at_infinity(c)
    return __knot_log_discount(c, t)
end
__log_native(::AbstractInterpolatedZeroCurve) = true

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
