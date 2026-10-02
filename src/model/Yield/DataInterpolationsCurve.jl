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
struct Spline{S <: Sp.SplineCurve, R, T, E, I, X} <: AbstractInterpolatedZeroCurve
    spline::S                  # the requested method (not the order reduced for short grids)
    rates::ReadOnlyVector{R}   # continuously-compounded zero rates at the knots
    tenors::ReadOnlyVector{T}  # finite, ≥ 0, strictly increasing
    extrapolation::E           # validated long-end policy
    _interp::I                 # t -> continuous zero rate through the last knot (and the flat short end)
    _tail::X                   # the long-end policy's `CurveTail` (unused with `:extension`)
    _extend::Bool              # `:extension`: the interpolant continues past the last knot
    # Internal: the builders below pass an owned grid and a validated policy. Every policy has the
    # same tail type, so a runtime policy Symbol does not change the curve's type.
    Spline(::Unchecked, spline::S, g::KnotGrid{R, T}, extrapolation::E, interp::I, tail::X, extend::Bool) where {S, R, T, E, I, X} =
        new{S, R, T, E, I, X}(spline, ReadOnlyVector(g.rates), ReadOnlyVector(g.tenors), extrapolation, interp, tail, extend)
end

# The knot-curve skeleton (InterpolatedCurves.jl) over the interpolant. The zero rate is a direct read
# of it, without a `-log(discount)/t` round trip (and its 0/0 at t = 0). L is inlined, so that the
# interpolant is evaluated in the body of `discount` rather than behind a second call (which made
# bond pricing on a Linear curve 40% slower). The forward is L's right-hand derivative: ForwardDiff's
# dual time at a knot evaluates the piece that starts there (DataInterpolations' own `derivative`
# takes the piece that ends there).
__tail(c::Spline) = c._tail
__extends(c::Spline) = c._extend
@inline __interior_zero(c::Spline, t) = c._interp(t)
@inline __interior_log_discount(c::Spline, t) = c._interp(t) * t
__interior_forward(c::Spline, t) = ForwardDiff.derivative(s -> __interior_log_discount(c, s), t)

# Public, validating form: every direct construction copies and checks its inputs.
Spline(spline::Union{Sp.PolynomialSpline, Sp.BSpline, Sp.PCHIP, Sp.Akima}, tenors, rates; extrapolation = :flat_forward) =
    __build_public(spline, KnotGrid(rates, tenors, spline; who = "Yield.Spline"); extrapolation)

# The one builder per interpolation method, over an owned grid: validated, or an internal
# optimizer trial grid. `ZeroRateCurve`, `Yield.Spline`, `reconstruct`, and `fit` all end here,
# so the descriptor → curve-type mapping lives in one place.
function __build(s::Sp.SplineCurve, g::KnotGrid; extrapolation = :flat_forward)
    interpolant = __interpolant(s, __interpolation_grid(g), __interpolation_extrapolation(extrapolation))
    return Spline(Unchecked(), s, g, extrapolation, interpolant, __spline_tail(interpolant, g, extrapolation), extrapolation === :extension)
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

# DataInterpolations handles the short end (a flat zero rate before the first knot) and, only
# for `:extension`, the long end; every other long-end policy is a `CurveTail` (Extrapolation.jl).
function __interpolation_extrapolation(method)
    E = DataInterpolations.ExtrapolationType
    return method === :extension ? (; extrapolation_left = E.Constant, extrapolation_right = E.Extension) :
        (; extrapolation_left = E.Constant)
end

# The average (discrete) continuously compounded forward over the last knot interval,
# `(zₙtₙ - zₙ₋₁tₙ₋₁) / (tₙ - tₙ₋₁)`: the `:flat_forward` anchor for DataInterpolations-backed
# curves. With a single knot, the only forward the data imply is the zero rate itself.
function __last_discrete_forward(rates, tenors)
    n = length(tenors)
    n == 1 && return only(rates)
    zₙ, tₙ, zₘ, tₘ = rates[n], tenors[n], rates[n - 1], tenors[n - 1]
    return (zₙ * tₙ - zₘ * tₘ) / (tₙ - tₘ)
end

# The long-end tail of a DataInterpolations-backed curve. `:extension` continues the interpolant
# instead, and gets a tail of the same type as `:flat_zero` that is never evaluated.
function __spline_tail(interpolant, g::KnotGrid, method)
    t, z = last(g.tenors), last(g.rates)
    method === :extension && return __curve_tail(t, z, zero(z), zero(z), false)
    forward = () -> __last_discrete_forward(g.rates, g.tenors)
    slope = () -> DataInterpolations.derivative(interpolant, t)
    return __build_tail(method, t, z, forward, slope)
end
