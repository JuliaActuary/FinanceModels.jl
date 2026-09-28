module Yield
import ..AbstractModel
import ..FinanceCore
import ..Spline as Sp
import ..ReadOnlyVector
import ..DataInterpolations
import ..Bond: coupon_times, __regular_schedule, __par_coupon
import ..__implicit_root, ..__primal, ..__ad_depth
import ..ForwardDiff

using ..FinanceCore: Continuous, Periodic, discount, accumulation, forward, pv, AbstractContract

export discount, zero, forward, par, implied_quote, pv, instantaneous_forward, knot_rates, knot_tenors, reconstruct

abstract type AbstractYieldModel <: AbstractModel end

# ── Internal curve capabilities ────────────────────────────────────────────────────────
#
# A curve needs only `discount(c, t)`. Each method below has a default that is safe for any such
# curve; a type adds its own where it knows more.
#
# | method                     | meaning                                        | default                        |
# |----------------------------|------------------------------------------------|--------------------------------|
# | `__log_discount(c, t)`     | L(t) = −log D(t)                               | `-log(discount(c, t))`         |
# | `__log_native(c)`          | L is direct, L(0) = 0, and D = exp(−L)         | `false`: the ratio D(b)/D(a)   |
# | `__log_interval(c, a, b)`  | L(b) − L(a), real where the interval factor is positive | from `__log_native`   |
# | `__log_tail(c)`            | the coefficients of L as t → ∞ (`__LogTail`)   | from `zero(c, Inf)`, below     |
#
# `__log_native` says that the log-difference route is valid for intervals: the curve's L is
# computed directly, with L(0) = 0, and its discount factor is exp(-L), so positive. It doesn't just
# mean that a type defines `__log_discount`: Smith–Wilson, `ForwardStarting`, Hull–White and the
# composite and scaled wrappers define L but take their intervals from `__log_interval`.
#
# - `__log_discount` + `__log_native`: `Constant`, `Spline`, `MonotoneConvex`, the zero-native curves
#   (`NelsonSiegel(Svensson)`, `CairnsPritchard(Extended)`, the yield shifts), `Vasicek`, `RatePath`.
# - `__log_discount` + `__log_interval`: `SmithWilson` (signed factors: from its interval ratio),
#   `ForwardStarting` (from the base curve's interval), `HullWhite` (its initial curve's),
#   `CompositeYield` and `ScaledYield` (their components' intervals, combined).
# - `__log_tail`, exact: `Constant`, `Spline`, `MonotoneConvex`, and the wrappers `CompositeYield`,
#   `ScaledYield`, `ForwardStarting`, `HullWhite`, which combine their components' tails. Every other
#   curve uses the fallback, which knows only the growth rate its zero rate at infinity implies.
# - `__PrimalCurve` forwards all four without partials.

# Cumulative log-discount L(t) = −log D(t): the continuously compounded force accumulated from
# valuation time 0 to `t`. `forward` and the generic `zero` are built from it. The fallback needs
# only `discount(c, t)`; built-in curves compute L directly.
__log_discount(c, t) = -log(FinanceCore.discount(c, t))

# Whether a curve's discount factor is exp(-L) for a directly computed L, and hence positive. Such
# curves take interval factors as exp(L(from) − L(to)): a difference of log-discounts stays finite
# where both discount factors underflow, and an interval from 0 reproduces `discount(c, t)` bit for
# bit. A curve that defines only `discount` may have D(0) ≠ 1 or a nonpositive discount factor, so
# its interval factors stay the ratio D(to)/D(from). The answer depends only on the curve's type.
__log_native(c) = false

# L(to) − L(from), the log-discount accumulated over an interval; `ForwardStarting` rebases a curve
# with it. A curve that isn't log-native takes the log of its own interval factor, which is real
# wherever that factor is positive (for a signed Smith-Wilson fit, also where both discount factors
# are negative).
__log_interval(c, from, to) = __log_native(c) ? __log_native_interval(c, from, to) :
    -log(FinanceCore.discount(c, from, to))

# L(to) − L(from) of a log-native curve. L(0) = 0, so an interval from 0 (every present value) skips
# it and is L(to) exactly. Under ForwardDiff 1.x `iszero` also requires zero partials, so a time
# derivative at 0 still evaluates L(from).
function __log_native_interval(c, from, to)
    L_to = __log_discount(c, to)
    return L_to - (iszero(from) ? zero(L_to) : __log_discount(c, from))
end

# The tail of L (see `__LogTail`) for a curve without a closed form: the leading slope is its zero
# rate at infinity and the intercept is not known. When those rates are finite, a composite of such
# curves has the same limit as the sum of their zero rates, and an exact cancellation of the slopes
# is NaN (0·Inf). An infinite zero rate means L grows faster than t, at an order the curve doesn't
# say; the forward slot then holds ±Inf, which the limit reads as that unknown growth.
function __log_tail(c)
    z = FinanceCore.rate(Base.zero(c, Inf))
    return __LogTail(zero(z), z, oftype(z, NaN))
end

# Wrappers evaluate L at t = +Inf from their combined tail; a negative time still reaches the
# components, which reject it. The limits are out of line so the finite-time path stays small.
__at_infinity(t) = isinf(t) && t > 0
@noinline __discount_at_infinity(c::AbstractYieldModel) = __discount_at_infinity(__log_tail(c))
@noinline __log_discount_at_infinity(c::AbstractYieldModel) = __log_discount_at_infinity(__log_tail(c))

# L of a zero-native curve: z(t)·t with z the continuous zero rate, which is exactly the exponent
# `discount(zero(c, t), t)` uses. At an exact zero time L is 0 whatever the zero rate there: a
# wrapper over a discount-native curve has a removable 0/0 zero rate at t = 0. Under ForwardDiff
# 1.x `iszero` also requires zero partials, so a time derivative at 0 still flows through z(t)·t.
function __zero_log_discount(c, t)
    L = FinanceCore.rate(Base.zero(c, t)) * t
    return iszero(t) ? zero(L) : L
end

# Discount factor of the curves defined by their zero rate (NelsonSiegel(Svensson),
# CairnsPritchard, MonotoneConvex, the yield shifts), shared through the one-line `discount`
# stubs at each curve definition. `exp(-L)` gives DF(0) = 1 exactly and keeps the promoted
# curve/time numeric type.
_discount_from_zero(c, t) = exp(-__zero_log_discount(c, t))

# Generic callable fallback: `curve(t) ≡ discount(curve, t)`. Covers every
# AbstractYieldModel subtype (Constant, Spline, CompositeYield, ScaledYield,
# TenorShift, ProjectedShift, NelsonSiegel, MonotoneConvex, …);
# each one routes through its own `discount`, so no per-type callable is needed.
(yc::AbstractYieldModel)(t) = FinanceCore.discount(yc, t)

"""
    Constant(rate)

A yield curve representing a flat term structure. `rate` can be a [`Rate`](@ref) object or a `Real` object.


If [`fit`](@ref FinanceModels.fit)ing with the default FinanceModels.jl settings, the solver will attempt to fit a discount rate with the range of: `-1.0 .. 1.0`
"""
struct Constant{R} <: AbstractYieldModel
    rate::R
end

function Constant(rate::R) where {R <: Real}
    return Constant(FinanceCore.Rate(rate))
end

Constant() = Constant(0.0)

FinanceCore.discount(c::Constant, t) = FinanceCore.discount(c.rate, t)

# The continuous zero rate of a flat curve is its (continuous) rate at every tenor,
# including t=0. Defining `zero` directly avoids the generic `-log(discount)/t`
# round-trip — which is `0/0 → NaN` at t=0 — and lets curves composed from a
# `Constant` stay in zero-rate space (see `CompositeYield`/`ScaledYield`).
Base.zero(c::Constant, t) = convert(Continuous(), c.rate)
__log_discount(c::Constant, t) = __zero_log_discount(c, t)
__log_native(::Constant) = true
function __log_tail(c::Constant)
    z = FinanceCore.rate(Base.zero(c, Inf))
    return __LogTail(zero(z), z, zero(z))
end

# ── Shared knot-grid construction ──────────────────────────────────────────────────────
#
# Every curve that interpolates zero rates over a knot grid (`Yield.Spline` and
# `Yield.MonotoneConvex`, however built: `ZeroRateCurve`, `reconstruct`, `fit`) obtains its
# knot data through `KnotGrid`, so all of them copy their inputs, promote to one concrete float
# type, and raise the same `ArgumentError`s for the same invalid grids.

# Minimum knot counts. Polynomial and B-spline orders reduce to `n - 1` on short grids, so one
# knot is a flat curve; PCHIP and Akima need three (two reach an UndefRefError inside
# DataInterpolations); a one-knot MonotoneConvex curve is flat.
__min_knots(::Sp.SplineCurve) = 1
__min_knots(::Sp.PCHIP) = 3
__min_knots(::Sp.Akima) = 3

__check_min_knots(::Nothing, n, who) = n >= 1 || throw(ArgumentError("$who: at least one knot is required (got $n)."))
function __check_min_knots(spline::Sp.SplineCurve, n, who)
    k = __min_knots(spline)
    n >= k || throw(ArgumentError("$who: $(spline) requires at least $k knots (got $n)."))
    return nothing
end


# Owned, concretely-typed float `Vector` from any iterable of reals. Always copies (even
# when handed a `Vector`), so the curve never aliases caller-owned memory; never narrows
# AD types (`float(::Dual)` is a `Dual`). Inputs that are not real numbers fail where they are
# converted (`float(String)`, `Float64(::Dual)` for mixed dual tags, …). An untyped empty input
# has no values to promote: `Bool`, which every real type absorbs, makes it `Float64[]`, so the
# grid's own length check reports it.
function __owned_float_vector(x)
    v = x isa AbstractVector ? x : collect(x)
    T = isconcretetype(eltype(v)) ? eltype(v) : mapreduce(typeof, promote_type, v; init = Bool)
    return Vector{float(T)}(v)   # Array-from-AbstractArray always allocates a fresh copy
end

# Marker for the internal, unvalidated `KnotGrid` construction used by optimizer trial curves.
struct Unchecked end

"""
    KnotGrid(rates, tenors[, spline::Spline.SplineCurve]; who = "KnotGrid")

Internal. The owned, validated knot data behind every `AbstractInterpolatedZeroCurve`
(`Yield.Spline`, `Yield.MonotoneConvex`). `rates` and `tenors` are each
copied from any iterable of reals (a `Vector` is copied too, so no curve aliases caller-owned
memory) and promoted **independently** to one concrete floating-point element type
(`Int` → `Float64`, `Float32` + `BigFloat` → `BigFloat`, `Float64` + `ForwardDiff.Dual` →
`Dual`, never narrowed).

Throws an `ArgumentError` (prefixed with `who`, the public constructor's name) when:

- `rates` and `tenors` differ in length, or either is empty;
- any rate or tenor is not finite (`NaN`, `±Inf`);
- any tenor is negative, or the tenors are not strictly increasing (unsorted or duplicated);
- there are fewer knots than `spline` needs (`__min_knots`). Without `spline`, one knot is
  enough; the curve constructors that accept a grid check their own minimum.

Inputs that are not real numbers fail where they are converted to floats (a `MethodError` such as
`float(::Type{String})`).

A curve built from a grid takes ownership of the grid's vectors. The grid is a construction
input, not a container to keep and modify.

`KnotGrid(Yield.Unchecked(), rates, tenors)` is internal: it neither copies nor validates. It
exists only for optimizer trial curves (`fit`, bootstrap) whose grid was validated up front
and whose candidate rates may be non-finite mid-search. Fitted results are always rebuilt
through the validating form.
"""
struct KnotGrid{R, T}
    rates::Vector{R}   # continuously-compounded zero rates (finite)
    tenors::Vector{T}  # finite, ≥ 0, strictly increasing
    KnotGrid(::Unchecked, rates::Vector{R}, tenors::Vector{T}) where {R, T} = new{R, T}(rates, tenors)
end
# unchecked trial form over non-`Vector` buffers (an optimizer may hand over a view or similar)
KnotGrid(u::Unchecked, rates::AbstractVector, tenors::AbstractVector) =
    KnotGrid(u, convert(Vector, rates), convert(Vector, tenors))

function KnotGrid(rates, tenors, spline::Union{Nothing, Sp.SplineCurve} = nothing; who = "KnotGrid")
    r = __owned_float_vector(rates)
    t = __owned_float_vector(tenors)
    length(r) == length(t) || throw(
        ArgumentError(
            "$who: `rates` and `tenors` must have the same length (got $(length(r)) and $(length(t)))."
        )
    )
    # the knot count first, so an empty grid reports it; then the tenors, which the rates of a
    # sampled grid depend on
    __check_min_knots(spline, length(t), who)
    all(isfinite, t) || throw(ArgumentError("$who: all tenors must be finite (got $(t))."))
    # primal values throughout: ForwardDiff orders Duals with equal values by their partials, so
    # `Dual(0, -1) >= 0` is false and `Dual(1, 1) > Dual(1, 0)` is true
    __primal(first(t)) >= 0 || throw(ArgumentError("$who: tenors must be ≥ 0 (got $(t))."))
    # finiteness is checked first so NaN cannot defeat the ordering comparison
    for i in 2:length(t)
        __primal(t[i]) > __primal(t[i - 1]) || throw(
            ArgumentError(
                "$who: tenors must be strictly increasing (sorted, no duplicates); got $(t). " *
                    "Sort the (rate, tenor) pairs before constructing."
            )
        )
    end
    all(isfinite, r) || throw(ArgumentError("$who: all rates must be finite (got $(r))."))
    return KnotGrid(Unchecked(), r, t)
end

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

# ── DataInterpolations-backed curve ────────────────────────────────────────────────────

"""
    Yield.Spline(spline::Spline.SplineCurve, tenors, rates; extrapolation=:flat_forward)

A yield curve that interpolates continuously-compounded zero rates over a knot grid with a
DataInterpolations interpolant chosen by the `spline` descriptor (`Spline.Linear()`,
`Spline.Cubic()`, `Spline.PCHIP()`, `Spline.BSpline(3)`, …). Polynomial and B-spline orders are
reduced to `length(tenors) - 1` when the grid is short, so a single knot gives a flat curve.
`Spline.MonotoneConvex()` is not a DataInterpolations method: build it with
[`ZeroRateCurve`](@ref) or [`Yield.MonotoneConvex`](@ref). Note the argument order
`(spline, tenors, rates)`, unlike `ZeroRateCurve(rates, tenors, spline)`, which returns the
same curve.

Before the first knot the zero rate is held flat at the first knot's rate.

The default `extrapolation=:flat_forward` holds the instantaneous forward rate constant beyond
the last knot at the **last discrete forward**, the average continuously compounded forward over
the last knot interval. For the last two knots `(tₙ₋₁, zₙ₋₁)` and `(tₙ, zₙ)`,

    fₙ = (zₙtₙ - zₙ₋₁tₙ₋₁) / (tₙ - tₙ₋₁),    z(t) = fₙ + (zₙ - fₙ)tₙ/t  for t > tₙ.

This preserves the last-knot discount factor (discount factors are continuous at `tₙ`) and does
not depend on the interpolant, so a quadratic or cubic end piece cannot drive the tail. The
instantaneous forward generally jumps at `tₙ`, from the interpolant's endpoint forward to `fₙ`.
Set `extrapolation` to `:flat_zero`, `:linear`, or `:extension` to instead hold the zero rate
constant, extend the zero rate at its left-hand boundary slope, or continue the final
interpolation piece, respectively.
[`FlatForwardAt`](@ref) holds a supplied forward rate fixed beyond the last knot.

Inputs are **copied** (later mutation of the vectors you passed in does not affect the curve)
and promoted to one concrete float element type each. An `ArgumentError` is raised for a length
mismatch, empty or non-finite inputs, negative/unsorted/duplicate tenors, or fewer knots than
the interpolant needs. `Yield.Spline` is a [`Yield.AbstractInterpolatedZeroCurve`](@ref): read
its knots with [`knot_rates`](@ref)/[`knot_tenors`](@ref) and change it with
[`reconstruct`](@ref).
"""
struct Spline{S <: Sp.SplineCurve, R, T, E, F} <: AbstractInterpolatedZeroCurve
    spline::S                  # the requested method (not the order reduced for short grids)
    rates::ReadOnlyVector{R}   # continuously-compounded zero rates at the knots
    tenors::ReadOnlyVector{T}  # finite, ≥ 0, strictly increasing
    extrapolation::E           # validated long-end policy
    _fn::F                     # t -> continuous zero rate (interpolant, short end, and tail)
    # Internal: the builders below pass an owned grid and a validated policy.
    Spline(::Unchecked, spline::S, g::KnotGrid{R, T}, extrapolation::E, fn::F) where {S, R, T, E, F} =
        new{S, R, T, E, F}(spline, ReadOnlyVector(g.rates), ReadOnlyVector(g.tenors), extrapolation, fn)
end

function __check_extension_at_infinity(c::Spline, t)
    c._fn.extend && throw(
        DomainError(
            t, "discount(curve, Inf) is unsupported with extrapolation = :extension, whose " *
                "polynomial continuation of the last piece is evaluated at finite times only. Use a " *
                "finite time, or a policy with a tail limit (:flat_forward, :flat_zero, :linear, FlatForwardAt)."
        )
    )
    return nothing
end

function FinanceCore.discount(c::Spline, t)
    __check_time(t, "discount")
    isinf(t) && return __spline_discount_at_infinity(c, t)
    return exp(-c._fn(t) * t)
end

function __log_discount(c::Spline, t)
    __check_time(t, "discount")
    isinf(t) && return __spline_log_discount_at_infinity(c, t)
    return c._fn(t) * t
end
# The limits are out of line so the finite-time path stays small enough to inline.
@noinline function __spline_discount_at_infinity(c::Spline, t)
    __check_extension_at_infinity(c, t)
    return __discount_at_infinity(c._fn.tail)
end
@noinline function __spline_log_discount_at_infinity(c::Spline, t)
    __check_extension_at_infinity(c, t)
    return __log_discount_at_infinity(c._fn.tail)
end
__log_native(::Spline) = true
function __log_tail(c::Spline)
    __check_extension_at_infinity(c, Inf)
    return __log_tail(c._fn.tail)
end

# `_fn(t)` is the continuous zero rate, so `zero` is a direct read of the interpolant. This
# avoids the generic `-log(discount)/t` round-trip (and its `0/0 → NaN` at t=0).
function Base.zero(c::Spline, t)
    __check_time(t, "zero")
    return Continuous(c._fn(t))
end

include("Yield/Extrapolation.jl")

# Public, validating form: every direct construction copies and checks its inputs.
Spline(spline::Sp.SplineCurve, tenors, rates; extrapolation = :flat_forward) =
    __build_public(spline, KnotGrid(rates, tenors, spline; who = "Yield.Spline"); extrapolation)
Spline(::Sp.MonotoneConvex, tenors, rates; extrapolation = :flat_forward) = throw(
    ArgumentError(
        "Yield.Spline implements the DataInterpolations methods only. Build a monotone convex " *
            "curve with ZeroRateCurve(rates, tenors, Spline.MonotoneConvex()) or Yield.MonotoneConvex(rates, tenors)."
    )
)

# The one builder per interpolation method, over an owned grid: validated, or an internal
# optimizer trial grid. `ZeroRateCurve`, `Yield.Spline`, `reconstruct`, and `fit` all end here,
# so the descriptor → curve-type mapping lives in one place.
function __build(s::Sp.SplineCurve, g::KnotGrid; extrapolation = :flat_forward)
    __check_min_knots(s, length(g.tenors), "Yield.Spline")
    interpolant = __interpolant(s, __interpolation_grid(g), __interpolation_extrapolation(extrapolation))
    return Spline(Unchecked(), s, g, extrapolation, __extrapolate(interpolant, g, extrapolation))
end

# DataInterpolations needs two points. A one-knot curve is flat, so its interpolant runs a flat
# segment from the knot to one year beyond it; the curve's own knots, tail, and short end are
# still those of the single knot. Returning a grid of the same type keeps construction inferable.
function __interpolation_grid(g::KnotGrid)
    length(g.tenors) == 1 || return g
    z, t = only(g.rates), only(g.tenors)
    return KnotGrid(Unchecked(), [z, z], [t, t + one(t)])
end

function __interpolant(b::Sp.BSpline, g::KnotGrid, extrapolation_kwargs)
    xs, ys = g.tenors, g.rates
    # documented short-grid rule (`Spline.BSpline`): the degree reduces to `k - 1` on `k` knots
    order = min(length(xs) - 1, b.order)
    knot_type = length(xs) < 3 ? :Uniform : :Average
    return DataInterpolations.BSplineInterpolation(ys, xs, order, knot_type; extrapolation_kwargs...)
end

__interpolant(::Sp.PCHIP, g::KnotGrid, extrapolation_kwargs) =
    DataInterpolations.PCHIPInterpolation(g.rates, g.tenors; extrapolation_kwargs...)

__interpolant(::Sp.Akima, g::KnotGrid, extrapolation_kwargs) =
    DataInterpolations.AkimaInterpolation(g.rates, g.tenors; extrapolation_kwargs...)

include("Yield/Kinks.jl")

function __interpolant(b::Sp.PolynomialSpline, g::KnotGrid, extrapolation_kwargs)
    xs, ys = g.tenors, g.rates
    # documented short-grid rule (`Spline.PolynomialSpline`): the order reduces to `k - 1` on `k` knots
    order = min(length(xs) - 1, b.order)
    # `cache_parameters = true` precomputes per-segment parameters at construction so that evaluation is
    # read-only and therefore thread-safe — notably for `QuadraticSpline`, which is B-spline-based and
    # otherwise overwrites a shared internal coefficient buffer on every call. It also avoids recomputing
    # parameters on each evaluation. (Curves here are built immutably, so the "do not mutate u/t" caveat
    # of cached parameters does not apply.)
    return if order == 1
        DataInterpolations.LinearInterpolation(ys, xs; extrapolation_kwargs..., cache_parameters = true)
    elseif order == 2
        DataInterpolations.QuadraticSpline(ys, xs; extrapolation_kwargs..., cache_parameters = true)
    else
        DataInterpolations.CubicSpline(ys, xs; extrapolation_kwargs..., cache_parameters = true)
    end
end

include("Yield/SmithWilson.jl")
include("Yield/NelsonSiegelSvensson.jl")
include("Yield/CairnsPritchard.jl")
include("Yield/MonotoneConvex.jl")

__build(::Sp.MonotoneConvex, g::KnotGrid; extrapolation = :flat_forward) = MonotoneConvex(g; extrapolation)

include("Yield/ZeroRateCurve.jl")
include("Yield/ImpliedQuote.jl")


## Generic and Fallbacks
"""
    discount(yc, to)
    discount(yc, from,to)

The discount factor for the yield curve `yc` for times `from` through `to`.
"""
function FinanceCore.discount(yc::T, from, to) where {T <: AbstractYieldModel}
    # A log-native interval from 0 is `discount(yc, to)` exactly (see `__log_native_interval`).
    d = __log_native(yc) ? exp(-__log_native_interval(yc, from, to)) :
        FinanceCore.discount(yc, to) / FinanceCore.discount(yc, from)
    # The empty interval is the identity, also where L is infinite at both ends (from = to = Inf).
    # Under ForwardDiff 1.x `==` also compares partials, so this fires only where the
    # derivative is zero anyway.
    return from == to ? one(d) : d
end

"""
    forward(yc, from, to)

The forward `Rate` implied by the yield curve `yc` between times `from` and `to`.
"""
function FinanceCore.forward(yc::T, from, to = from + 1) where {T <: AbstractYieldModel}
    # forward = log(DF(from)/DF(to)) / (to-from) = (L(to) − L(from))/(to−from), with L the
    # cumulative log-discount (z(t)·t for a zero-native curve), which is 0 at t = 0.
    return Continuous((__log_discount(yc, to) - __log_discount(yc, from)) / (to - from))
end

"""
    par(curve,time;frequency=2)

Calculate the par yield for maturity `time` for the given `curve` and `frequency`. Returns a `Rate` object with periodicity corresponding to the `frequency`.

If `time` is shorter than one regular coupon period (e.g. `time=0.5` with `frequency=1`), the single stub payment implies a compounding frequency of `1/time`: the result is quoted as `Periodic(1/time)` when `1/time` is a (near-)integer, and otherwise an `ArgumentError` is thrown because the implied frequency cannot be represented as a `Periodic` rate.

If `time` is longer than one coupon period but not a whole number of periods, the schedule has a short first stub which accrues its actual length (see `Bond.coupon_times`); the result is the internal rate of return of the par-priced true-accrual schedule, quoted as `Periodic(frequency)` — so `par` of a flat curve recovers the curve's rate at any maturity. On such stub schedules this yield quote differs from the annualized par *coupon* `c` solving `c·Σᵢ δᵢ·DF(tᵢ) + DF(T) = 1`, which is what `InterestRateSwap` uses for its fixed leg.

# Examples

```julia-repl
julia> c = Yield.Constant(0.04);

julia> par(c,4)
Periodic(0.03960780543711406, 2)

julia> par(c,4;frequency=1)
Periodic(0.040000000000000036, 1)

julia> par(c,0.6;frequency=4)
Periodic(0.039413626195875295, 4)

julia> par(c,0.2;frequency=4)
Periodic(0.039374942589460726, 5)

julia> par(c,2.5)
Periodic(0.03960780543711406, 2)
```
"""
function par(curve, time; frequency = 2)
    coup_times = coupon_times(time, frequency)
    mat_disc = discount(curve, time)
    coupon_pv = sum(discount(curve, t) for t in coup_times)
    Δt = step(coup_times)
    r = (1 - mat_disc) / coupon_pv

    # A maturity that is not a whole number of periods leaves a short first stub
    # which accrues its actual length: such schedules pay `c·δᵢ` with the
    # annualized coupon `c` solving the true-accrual par condition, so the IRR
    # below is that of the correctly-accrued par-priced bond. Whole-period
    # schedules pay `r` every period (the historical construction, unchanged).
    stub = length(coup_times) > 1 && !__regular_schedule(time, frequency)
    c = stub ? __par_coupon(curve, time, frequency) : r

    # Build cash flows: initial outflow of -1, then coupons, final coupon+principal 1+coupon
    # Pre-allocate arrays for better performance
    n = length(coup_times)
    cfs = Vector{typeof(r)}(undef, n + 1)
    times = Vector{typeof(Δt)}(undef, n + 1)

    cfs[1] = -one(r)
    times[1] = zero(Δt)

    @inbounds for i in 1:n
        coup = stub ? c * (i == 1 ? first(coup_times) : 1 / frequency) : r
        cfs[i + 1] = i == n ? 1 + coup : coup
        times[i + 1] = coup_times[i]
    end

    r = FinanceCore.internal_rate_of_return(cfs, times)
    frequency_inner = 1 / Δt  # Simplified from min(1 / Δt, max(1 / Δt, frequency))
    if !isinteger(round(frequency_inner, digits = 8))
        throw(
            ArgumentError(
                "par(curve, $time; frequency=$frequency) implies a coupon period of $Δt and a compounding frequency of 1/Δt = $frequency_inner, which is not an integer and cannot be represented as a `Periodic` rate. Choose a maturity commensurate with the coupon frequency."
            )
        )
    end
    r = convert(Periodic(frequency_inner), r)
    return r
end

"""
    zero(curve,time)

Return the zero rate for the curve at the given time.
"""
function Base.zero(c::YC, time) where {YC <: AbstractYieldModel}
    # L/t; for a curve that defines only `discount`, L is -log(discount(c, time))
    return Continuous(__log_discount(c, time) / time)
end

"""
    accumulation(yc, from, to)

The accumulation factor for the yield curve `yc` for times `from` through `to`.
"""
function FinanceCore.accumulation(yc::AbstractYieldModel, time)
    return 1 ./ discount(yc, time)
end

# the reversed interval, so each curve's own interval method applies
FinanceCore.accumulation(yc::AbstractYieldModel, from, to) = FinanceCore.discount(yc, to, from)

## Curve Manipulations
"""
    CompositeYield(curve1, curve2, op)

Combines two yield curves by adding (`op = +`) or subtracting (`op = -`) their continuous zero
rates. Created via `+` and `-` on `AbstractYieldModel` objects; for scalar multiplication or
division, see [`ScaledYield`](@ref).

Given discount factors `DF₁(t)` and `DF₂(t)` with continuous zero rates `z₁` and `z₂`, the
composite discount factor is `exp(-(z₁ ± z₂) t)`:

- `+` gives `DF(t) = DF₁(t) × DF₂(t)`, the product of the two discount factors (for example a
  base curve and a spread);
- `-` gives `DF(t) = DF₁(t) / DF₂(t)`, their quotient.

These are the operations that compose the curves' factors: every interval factor of the result
is the product (or quotient) of the components' interval factors, so, for example,
`ForwardStarting(a + b, τ)` prices like `ForwardStarting(a, τ) + ForwardStarting(b, τ)`. Other
operations on the zero rates (`max`, `*`, …) are not accepted: they would build a new curve from
zero rates measured from time 0, not a composition of the two curves. For a pointwise
transformation of a curve's zero rates, use a [`TenorShift`](@ref), `curve + ((z, t) -> ...)`,
or define a curve type with its own `zero`.

Composition is performed in continuous-zero-rate space: a `+`/`-` composite reads each
component's zero rate, combines them, and applies a single `exp` to form the discount
factor — it no longer pays the `log`/`exp` round-trip that earlier versions did. Composing
many curves in a hot loop is still marginally slower than pre-fitting a single combined
curve, but the gap is small.

Curves can be added or subtracted together, but note that this is not always the same thing
as adding or subtracting spreads with rates. If spreads and base rates are expressed as zero
rates, then the curve addition/subtraction has the same effect as re-fitting the yield model
with the rate+spread inputs added together first. Non-zero rates (e.g. par rates) do not have
this same property.

## Examples

```julia
rates = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
spreads = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
mats = [1 / 12, 2 / 12, 3 / 12, 6 / 12, 1, 2, 3, 5, 7, 10, 20, 30]


### Zero coupon rates/spreads

q_rf_z = ZCBYield.(rates,mats)
q_s_z = ZCBYield.(spreads,mats)
q_y_z = ZCBYield.(rates + spreads,mats)

c_rf_z = fit(Spline.Linear(),q_rf_z,Fit.Bootstrap())
c_s_z = fit(Spline.Linear(),q_s_z,Fit.Bootstrap())
c_y_z = fit(Spline.Linear(),q_y_z,Fit.Bootstrap())

# adding curves when the spreads were zero spreads works
@test discount(c_rf_z+c_s_z,20) ≈ discount(c_y_z,20)


### Par coupon rates/spreads

q_rf = CMTYield.(rates,mats)
q_s = CMTYield.(spreads,mats)
q_y = CMTYield.(rates + spreads,mats)

c_rf = fit(Spline.Linear(),q_rf,Fit.Bootstrap())
c_s = fit(Spline.Linear(),q_s,Fit.Bootstrap())
c_y = fit(Spline.Linear(),q_y,Fit.Bootstrap())

# adding curves when the spreads were par spreads does not work
@test !(discount(c_rf+c_s,20) ≈ discount(c_y,20))
```
"""
struct CompositeYield{T, U, V <: Union{typeof(+), typeof(-)}} <: AbstractYieldModel
    r1::T
    r2::U
    op::V
end


# Composition happens in log-discount space: `+` adds the components' cumulative log-discounts
# (multiplying their discount factors) and `-` subtracts them (dividing), then a single `exp`
# forms the discount factor. The zero rate is the same combination of the components' zero rates.
function Base.zero(rc::CompositeYield, time)
    z1 = FinanceCore.rate(Base.zero(rc.r1, time))
    z2 = FinanceCore.rate(Base.zero(rc.r2, time))
    return Continuous(rc.op(z1, z2))
end
# At t = Inf the components' tails are combined before the limit is taken. The components' L is
# inlined at these calls: a composite inlined into a loop (`pv`'s `map`) otherwise copies each
# component onto the stack for every call. Inlining a Spline's L into every caller instead makes
# the generic interval too large to inline into bond pricing.
function __log_discount(rc::CompositeYield, time)
    __at_infinity(time) && return __log_discount_at_infinity(rc)
    return rc.op((@inline __log_discount(rc.r1, time)), (@inline __log_discount(rc.r2, time)))
end
function FinanceCore.discount(rc::CompositeYield, time)
    __at_infinity(time) && return __discount_at_infinity(rc)
    return exp(-__log_discount(rc, time))
end
__log_tail(rc::CompositeYield) = __combine_tails(rc.op, __log_tail(rc.r1), __log_tail(rc.r2))

# A wrapper's interval combines its components' intervals, so each keeps its own form: a log-native
# curve's difference of log-discounts (L(to) alone from 0), Smith–Wilson's signed ratio, a rebased
# curve's base interval, the ratio D(to)/D(from) of a curve that defines only `discount` (whose D(0)
# need not be 1). Those stay finite where both of a wrapper's discount factors underflow. At an
# infinite endpoint the components' limits alone lose information (flat forwards of 4% and −2% have
# L = Inf and −Inf), so the wrapper's own L there comes from its combined tail.
__at_infinity_either(from, to) = isinf(from) || isinf(to)
function __log_interval(rc::CompositeYield, from, to)
    __at_infinity_either(from, to) && return __log_discount(rc, to) - __log_discount(rc, from)
    return rc.op(__log_interval(rc.r1, from, to), __log_interval(rc.r2, from, to))
end

# ─── Yield shifts (TenorShift, ProjectedShift) ────────────────────────────

"""
    AbstractYieldShift <: AbstractYieldModel

Supertype for lazy zero-rate shift models: a curve produced by transforming
a base yield curve's zero rate via a user-supplied rule.

Two concrete subtypes:

- [`TenorShift`](@ref) — shift depends only on the tenor `t`. Use for parallel
  bumps, twists, butterflies — static curve transformations.
- [`ProjectedShift`](@ref) — shift depends on the tenor `t` *and* on a second
  time axis `τ` (projection / as-of / valuation-date time). Use for phase-in
  profiles (BMA SBA, IFRS17 macro scenarios) and any shift whose shape evolves
  across a projection horizon.

Both subtypes implement the standard `AbstractYieldModel` interface (`zero`,
`discount`, `forward`, `pv`).
"""
abstract type AbstractYieldShift <: AbstractYieldModel end

"""
    TenorShift(base, rule)

Lazy zero-rate transformation depending only on the tenor:
`z_new(t) = rule(z_base(t), t)`.

The `rule` function receives the base curve's `Continuous` zero rate and
the tenor, and returns a new rate. It is evaluated on demand — no discretization
or refitting. The base curve's analytic structure is fully preserved.

The rule function must have the signature `(z::Rate, t) -> Rate`. The return
value is type-asserted as `Rate` so that compounding convention is always
carried explicitly; rules returning a plain `Real` will raise a `TypeError`
at call time. Use `z + Continuous(0.01)`, `Periodic(0.04, 2)`, etc., and let
`Rate` arithmetic handle conversion to the curve's continuous representation.

Use this for static shifts — parallel bumps, twists, butterflies — that don't
depend on where you are in projection time. For shifts whose shape evolves
across a projection horizon, see [`ProjectedShift`](@ref).

# Constructing

The most ergonomic way to create a `TenorShift` is via the `+` operator
with an `AbstractYieldModel` and a two-argument function:

```julia
base = Yield.Constant(0.05)

# Parallel shift (+100 bp)
base + (z, t) -> z + Periodic(0.01, 1)

# Tenor-dependent twist (steepener that fades at 30y)
base + (z, t) -> z + Continuous(0.02 * max(0.0, 1.0 - t/30.0))
```

You can also construct directly:

```julia
TenorShift(base, (z, t) -> z + Continuous(0.01))
```

Note: The `+` operator dispatches on `Function`. For callable objects that are
not `Function` subtypes (e.g. custom structs with call syntax), use the direct
constructor: `TenorShift(base, my_callable)`.

`TenorShift` is a post-processing wrapper — it is not a fitting target.
ForwardDiff propagates correctly through the transform for sensitivity analysis,
but the rule function itself should be differentiable if used in an AD context.

`TransformedYield` is retained as a deprecated alias for `TenorShift`.

See also: [`ProjectedShift`](@ref), [`AbstractYieldShift`](@ref),
[`CompositeYield`](@ref), [`ScaledYield`](@ref).
"""
struct TenorShift{C <: AbstractYieldModel, F} <: AbstractYieldShift
    base::C
    rule::F
end

function Base.zero(s::TenorShift, t)
    z = Base.zero(s.base, t)
    return convert(Continuous(), s.rule(z, t)::FinanceCore.Rate)
end

"""
    ProjectedShift(base, rule, time)

Lazy zero-rate transformation that depends on the tenor `t` *and* a second
time axis `τ`:

    z_new(t) = rule(τ, z_base(t), t)

`τ` (stored as the `.time` field) is the **projection time** — the as-of or
valuation-date offset at which this curve is being evaluated. It is distinct
from the tenor `t`, which is time-to-maturity from `τ`.

The rule function must have the signature `(τ, z::Rate, t) -> Rate`. The
return value is type-asserted as `Rate` so that compounding convention is
always carried explicitly; rules returning a plain `Real` will raise a
`TypeError` at call time.

Use this for shifts whose shape evolves across a projection horizon — phase-in
profiles (BMA SBA, IFRS17 macro scenarios), runoff schedules, calendar-rolling
shocks. For static, tenor-only shifts, see [`TenorShift`](@ref).

# Constructing

There is no `+` operator sugar — `ProjectedShift` needs an explicit `τ`, which
fixing at composition time would defeat the purpose of storing the rule as a
year-independent first-class value. Always use the direct constructor:

```julia
base = Yield.Constant(0.05)

# −150 bp parallel shift, phased in linearly over 10 projection years.
phase_in = (τ, z, _) -> z + Continuous(-0.015 * min(τ, 10) / 10)

# Curve as seen at projection year 3 (30% phased in → -45 bp).
c3 = ProjectedShift(base, phase_in, 3.0)

# Curve as seen at projection year 10 (fully phased in → -150 bp).
c10 = ProjectedShift(base, phase_in, 10.0)
```

The intended pattern: store `rule` once as a first-class value, then call
`ProjectedShift(base, rule, τ)` at each `τ` in a projection loop.

See also: [`TenorShift`](@ref), [`AbstractYieldShift`](@ref).
"""
struct ProjectedShift{C <: AbstractYieldModel, F, T} <: AbstractYieldShift
    base::C
    rule::F
    time::T
end

function Base.zero(s::ProjectedShift, t)
    z = Base.zero(s.base, t)
    return convert(Continuous(), s.rule(s.time, z, t)::FinanceCore.Rate)
end
FinanceCore.discount(s::AbstractYieldShift, t) = _discount_from_zero(s, t)
__log_discount(s::AbstractYieldShift, t) = __zero_log_discount(s, t)
__log_native(::AbstractYieldShift) = true

# Deprecated alias for the previous name. Slated for removal one minor release after introduction.
Base.@deprecate_binding TransformedYield TenorShift

"""
    ScaledYield(curve, factor)

A yield model that scales the continuous zero rates of `curve` by a `Real` scalar `factor`.

Created via `curve * scalar` or `curve / scalar`. For example, `curve * 0.79` scales
all continuous zero rates by 0.79, which is useful for after-tax yield calculations.
"""
struct ScaledYield{T <: AbstractYieldModel, S <: Real} <: AbstractYieldModel
    curve::T
    factor::S
end

# Scaling multiplies the zero rate, and so the cumulative log-discount, by `factor` (it raises
# the discount factor to that power); the discount factor is a single `exp` of it.
function Base.zero(sy::ScaledYield, time)
    z = FinanceCore.rate(Base.zero(sy.curve, time))
    return Continuous(z * sy.factor)
end
# The curve's L is inlined at this call, as for `CompositeYield`.
function __log_discount(sy::ScaledYield, time)
    __at_infinity(time) && return __log_discount_at_infinity(sy)
    return sy.factor * @inline(__log_discount(sy.curve, time))
end
function FinanceCore.discount(sy::ScaledYield, time)
    __at_infinity(time) && return __discount_at_infinity(sy)
    return exp(-__log_discount(sy, time))
end
__log_tail(sy::ScaledYield) = __scale_tail(sy.factor, __log_tail(sy.curve))
function __log_interval(sy::ScaledYield, from, to)
    __at_infinity_either(from, to) && return __log_discount(sy, to) - __log_discount(sy, from)
    return sy.factor * __log_interval(sy.curve, from, to)
end
function FinanceCore.discount(w::Union{CompositeYield, ScaledYield}, from, to)
    d = exp(-__log_interval(w, from, to))
    return from == to ? one(d) : d
end

"""
    ForwardStarting(curve,forwardstart)

Rebase a `curve` so that `discount`/`accumulation`/etc. are re-based so that time zero from the new curves perspective is the given `forwardstart` time.

# Examples

```julia-repl
julia> zero = [5.0, 5.8, 6.4, 6.8] ./ 100
julia> maturity = [0.5, 1.0, 1.5, 2.0]
julia> curve = ZeroRateCurve(zero, maturity)
julia> fwd = Yield.ForwardStarting(curve, 1.0)

julia> discount(curve,1,2)
0.9275624570410582

julia> discount(fwd,1) # `curve` has effectively been reindexed to `1.0`
0.9275624570410582
```

# Extended Help

While `ForwardStarting` could be nested so that, e.g. the third period's curve is the one-period forward of the second period's curve, it will be more efficient to reuse the initial curve from a runtime and compiler perspective.

`ForwardStarting` is not used to construct a curve based on forward rates. 
"""
struct ForwardStarting{T, U} <: AbstractYieldModel
    curve::U
    forwardstart::T
end

# The rebased curve's factors are the base curve's interval factors from `forwardstart`, so the
# base curve's own interval method applies (finite in the far tail for the built-in curves).
FinanceCore.discount(c::ForwardStarting, to) = FinanceCore.discount(c.curve, c.forwardstart, to + c.forwardstart)
FinanceCore.discount(c::ForwardStarting, from, to) =
    FinanceCore.discount(c.curve, from + c.forwardstart, to + c.forwardstart)
# L of the rebased curve is the base curve's log-discount over [forwardstart, forwardstart + t], so it
# is real wherever the rebased discount factor is positive, even where the base curve's own discount
# factors are negative.
function __log_discount(c::ForwardStarting, t)
    __at_infinity(t) && return __log_discount_at_infinity(c)
    return __log_interval(c.curve, c.forwardstart, t + c.forwardstart)
end
__log_interval(c::ForwardStarting, from, to) = __log_interval(c.curve, from + c.forwardstart, to + c.forwardstart)
__log_tail(c::ForwardStarting) = __shift_tail(__log_tail(c.curve), c.forwardstart, __log_discount(c.curve, c.forwardstart))

"""
    Yield.AbstractYieldModel + Yield.AbstractYieldModel

The addition of two yields will create a `CompositeYield`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will be added together.
"""
function Base.:+(a::AbstractYieldModel, b::AbstractYieldModel)
    return CompositeYield(a, b, +)
end

function Base.:+(a::Constant, b::Constant)
    z_a = FinanceCore.rate(convert(Continuous(), a.rate))
    z_b = FinanceCore.rate(convert(Continuous(), b.rate))
    return Constant(Continuous(z_a + z_b))
end

function Base.:+(a::T, b::Union{Real, Rate}) where {T <: AbstractYieldModel}
    return a + Constant(b)
end

function Base.:+(a::Union{Real, Rate}, b::T) where {T <: AbstractYieldModel}
    return Constant(a) + b
end

function Base.:+(a::AbstractYieldModel, f::Function)
    return TenorShift(a, f)
end

function Base.:+(f::Function, a::AbstractYieldModel)
    return TenorShift(a, f)
end

"""
    curve * scalar
    scalar * curve

Scale the continuous zero rates of `curve` by a `Real` scalar. Returns a [`ScaledYield`](@ref).

This is useful for after-tax yield calculations. For example, `curve * 0.79` produces a
curve whose continuous zero rate at every point is 79% of the original.

# Examples

```julia-repl
julia> m = Yield.Constant(Continuous(0.05)) * 0.79;

julia> discount(m, 1) ≈ exp(-0.05 * 0.79)
true
```
"""
function Base.:*(a::AbstractYieldModel, b::Real)
    return ScaledYield(a, b)
end

function Base.:*(a::Real, b::AbstractYieldModel)
    return ScaledYield(b, a)
end

"""
    Yield.AbstractYieldModel - Yield.AbstractYieldModel

The subtraction of two yields will create a `CompositeYield`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the second curves will be subtracted from the first.
"""
function Base.:-(a::AbstractYieldModel, b::AbstractYieldModel)
    return CompositeYield(a, b, -)
end

function Base.:-(a::Constant, b::Constant)
    z_a = FinanceCore.rate(convert(Continuous(), a.rate))
    z_b = FinanceCore.rate(convert(Continuous(), b.rate))
    return Constant(Continuous(z_a - z_b))
end

function Base.:-(a::T, b::Union{Real, Rate}) where {T <: AbstractYieldModel}
    return a - Constant(b)
end

function Base.:-(a::Union{Real, Rate}, b::T) where {T <: AbstractYieldModel}
    return Constant(a) - b
end

"""
    curve / scalar

Scale the continuous zero rates of `curve` by `1/scalar`. Returns a [`ScaledYield`](@ref).

This is useful for grossing-up a yield to a pre-tax equivalent.

# Examples

```julia-repl
julia> m = Yield.Constant(Continuous(0.05)) / 0.79;

julia> discount(m, 1) ≈ exp(-0.05 / 0.79)
true
```
"""
function Base.:/(a::AbstractYieldModel, b::Real)
    return ScaledYield(a, inv(b))
end


end
