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
and promoted to one concrete float element type each. `rates` takes numbers (continuously
compounded), `Rate`s, or a mixture, as [`ZeroRateCurve`](@ref) does. An `ArgumentError` is
raised for a length mismatch, empty or non-finite inputs, negative/unsorted/duplicate tenors, or
fewer knots than the interpolant needs. `Yield.Spline` is a [`Yield.AbstractInterpolatedZeroCurve`](@ref): read
its knots with [`knot_rates`](@ref)/[`knot_tenors`](@ref) and change it with
[`reconstruct`](@ref).
"""
struct Spline{S <: Sp.SplineCurve, R, T, E, I, X} <: AbstractInterpolatedZeroCurve
    spline::S                     # the requested method (not the order reduced for short grids)
    _rates::ReadOnlyVector{R, R}  # continuously-compounded zero rates at the knots; property `rates` is `knot_rates`
    tenors::ReadOnlyVector{T, T}  # finite, ≥ 0, strictly increasing
    extrapolation::E              # validated long-end policy
    _interp::I                    # t -> continuous zero rate through the last knot (and the flat short end)
    _tail::X                      # the long-end policy's `CurveTail` (unused with `:extension`)
    _extend::Bool                 # `:extension`: the interpolant continues past the last knot
    # Internal: the builders below pass an owned grid and a validated policy. Every policy has the
    # same tail type, so a runtime policy Symbol does not change the curve's type.
    Spline(::Unchecked, spline::S, g::KnotGrid{R, T}, extrapolation::E, interp::I, tail::X, extend::Bool) where {S, R, T, E, I, X} =
        new{S, R, T, E, I, X}(spline, ReadOnlyVector(g.rates), ReadOnlyVector(g.tenors), extrapolation, interp, tail, extend)
end

# The knot-curve skeleton (InterpolatedCurves.jl) over the interpolant. The zero rate is a direct read
# of it, with no L(t)/t round trip. L is inlined, so that the
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
# so the descriptor → curve-type mapping lives in one place. `interpolation` is
# `__interpolation(s, length(g.tenors))`.
function __build(s::Sp.SplineCurve, interpolation, g::KnotGrid; extrapolation = :flat_forward)
    interpolant = __interpolant(interpolation, __interpolation_grid(g), __interpolation_extrapolation(extrapolation))
    return Spline(Unchecked(), s, g, extrapolation, interpolant, __spline_tail(interpolant, g, extrapolation), extrapolation === :extension)
end
__build(s::Sp.SplineCurve, g::KnotGrid; extrapolation = :flat_forward) = __build(s, __interpolation(s, length(g.tenors)), g; extrapolation)
# Branching on the order keeps a public polynomial-spline curve's type inferable as a union of the
# three curve types.
function __build(b::Sp.PolynomialSpline, g::KnotGrid; extrapolation = :flat_forward)
    order = __polynomial_order(b, length(g.tenors))
    order == 1 && return __build(b, Val(1), g; extrapolation)
    order == 2 && return __build(b, Val(2), g; extrapolation)
    return __build(b, Val(3), g; extrapolation)
end

# The interpolation a descriptor uses on `n` knots. Code that builds many curves on one grid (a
# fit's trial curves and Jacobian, a bootstrap's steps) resolves it once and builds through it, so
# that each of those curves has one concrete type. A polynomial spline's order is a runtime field,
# so its interpolation is a `Val` of the order: branching on the order inside an AD-traced loss
# made every trial curve a union of the linear, quadratic, and cubic curve types, and all three
# were compiled for each width of dual number.
__interpolation(s::Sp.SplineCurve, n) = s
__interpolation(b::Sp.PolynomialSpline, n) = Val(__polynomial_order(b, n))
# A built curve's interpolation, read from its type so that it is inferable.
__interpolation(c::AbstractInterpolatedZeroCurve) = __interpolation(c.spline, length(knot_tenors(c)))
__interpolation(c::Spline{<:Sp.PolynomialSpline}) = __polynomial_interpolation(c._interp)
__polynomial_interpolation(::DataInterpolations.LinearInterpolation) = Val(1)
__polynomial_interpolation(::DataInterpolations.QuadraticSpline) = Val(2)
__polynomial_interpolation(::DataInterpolations.CubicSpline) = Val(3)

# documented short-grid rule (`Spline.PolynomialSpline`): the order reduces to `k - 1` on `k`
# knots, and a single knot is interpolated over two points (`__interpolation_grid`)
__polynomial_order(b::Sp.PolynomialSpline, n) = min(max(n, 2) - 1, b.order)

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

# A polynomial spline of the resolved order (`__interpolation`). `cache_parameters = true` precomputes
# per-segment parameters at construction so that evaluation is read-only and therefore thread-safe —
# notably for `QuadraticSpline`, which is B-spline-based and otherwise overwrites a shared internal
# coefficient buffer on every call. It also avoids recomputing parameters on each evaluation. (Curves
# here are built immutably, so the "do not mutate u/t" caveat of cached parameters does not apply.)
__interpolant(::Val{1}, g::KnotGrid, extrapolation_kwargs) =
    DataInterpolations.LinearInterpolation(g.rates, g.tenors; extrapolation_kwargs..., cache_parameters = true)
__interpolant(::Val{2}, g::KnotGrid, extrapolation_kwargs) =
    DataInterpolations.QuadraticSpline(g.rates, g.tenors; extrapolation_kwargs..., cache_parameters = true)
__interpolant(::Val{3}, g::KnotGrid, extrapolation_kwargs) =
    DataInterpolations.CubicSpline(g.rates, g.tenors; extrapolation_kwargs..., cache_parameters = true)

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

# The interpolant's left-hand slope at the last knot `t`: `DataInterpolations.derivative(interpolant, t)`
# at a time inside the grid, which is `_derivative(interpolant, t, interpolant.iguesser)` (DataInterpolations
# 9.2 through 10.1). `derivative` branches on a runtime `order`, so calling it compiled its second-order
# branch, a ForwardDiff derivative of the interpolant at a dual time, for every interpolant and element type
# a curve is built with: each policy's curve compiles this slope, `:linear` or not. That made the first
# gradient through `ZeroRateCurve(rates, tenors, Spline.Linear())` take 6.8 s instead of 1.9 s.
__boundary_slope(interpolant, t) = DataInterpolations._derivative(interpolant, t, interpolant.iguesser)

# The long-end tail of a DataInterpolations-backed curve. `:extension` continues the interpolant
# instead, and gets a tail of the same type as `:flat_zero` that is never evaluated.
function __spline_tail(interpolant, g::KnotGrid, method)
    t, z = last(g.tenors), last(g.rates)
    method === :extension && return __curve_tail(t, z, zero(z), zero(z), false)
    forward = () -> __last_discrete_forward(g.rates, g.tenors)
    slope = () -> __boundary_slope(interpolant, t)
    return __build_tail(method, t, z, forward, slope)
end
