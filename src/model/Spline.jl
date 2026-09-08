"""
Spline is a module which offers various degree splines used for fitting or bootstraping curves via the [`fit`](@ref FinanceModels.fit) function.

Available methods:

- `Spline.PolynomialSpline(n)` where n is the nth order. A piecewise polynomial interpolating spline (order 1/2/3 → linear / quadratic / natural cubic). These are fast and thread-safe to evaluate. Knot-rate sensitivities depend on the order: linear is local within the knot grid, while natural cubic coefficients depend on the whole grid.
- `Spline.BSpline(d)` where d is the polynomial degree. A degree-d B-spline produces (d-1)th-order-continuous piecewise polynomials. That is, degree 2/3 is very similar to a quadratic/cubic spline respectively. BSplines are global in that a change in one point affects the entire spline (though the spline still passes through the other given points still). Useful as a basis for least-squares fitting, but **not** thread-safe for concurrent evaluation — see [`Spline.BSpline`](@ref).

This object is not a fitted spline itself, rather it is a placeholder object which will be a spline representing the data only after using within [`fit`](@ref FinanceModels.fit).

Convenience methods which create a `Spline.PolynomialSpline` of the appropriate order (recommended for interpolating curves — fast and safe to evaluate concurrently):

- `Spline.Linear()` equals `PolynomialSpline(1)` (numerically identical to `BSpline(1)`)
- `Spline.Quadratic()` equals `PolynomialSpline(2)`
- `Spline.Cubic()` equals `PolynomialSpline(3)`

For a *global* B-spline (e.g. as a basis for smooth least-squares fitting) use `Spline.BSpline(d)` explicitly, noting its thread-safety caveat.

Knot-based yield curves built from these descriptors extrapolate flat-forward beyond their
last knot by default. Pass `extrapolation=:flat_zero`, `:linear`, `:extension`, or
`Yield.FlatForwardAt(forward)` to
[`Yield.Spline`](@ref FinanceModels.Yield.Spline), `ZeroRateCurve`, or spline `fit` methods to
select a different long-end policy. `:extension` is unavailable for MonotoneConvex.

Notes on Fitting:
- `fit(spline,quotes)` will fit entire curve at once, with knots equal to the maturity points of the `Quote`s
- `fit(spline, quotes, Fit.Bootstrap())` will curve one knot at a time, with knots equal to the maturity points of the `Quote`s

Generally, the former will be preferred for performance reasons.

## Examples

```julia
using FinanceModels
using BenchmarkTools
rates = [0.07, 0.16, 0.35, 0.92, 1.4, 1.74, 2.31, 2.41] ./ 100
mats = [1, 2, 3, 5, 7, 10, 20, 30]

qs = CMTYield.(rates, mats)
c = fit(Spline.Linear(), qs) # will fit entire curve at once, with knots equal to the maturity points of the `Quote`s
c = fit(Spline.Linear(), qs, Fit.Bootstrap()) # will curve one knot at a time, with knots equal to the maturity points of the `Quote`s

```
"""
module Spline
import ..FinanceCore
import ..AbstractModel

abstract type SplineCurve end

"""
    Spline.PolynomialSpline(order)

A piecewise polynomial interpolating spline of the given `order`, backed by DataInterpolations
(`order` 1 → `LinearInterpolation`, 2 → `QuadraticSpline`, 3 → natural `CubicSpline`).
These interpolants are thread-safe to evaluate concurrently. Linear interpolation has local
knot-rate sensitivities within the knot grid. Natural cubic coefficients depend on the whole
grid, so bumping one knot can affect other intervals and the default long-end forward.

The convenience constructors [`Spline.Linear`](@ref), [`Spline.Quadratic`](@ref), and [`Spline.Cubic`](@ref)
return `PolynomialSpline(1/2/3)`.
"""
struct PolynomialSpline <: SplineCurve
    order::Int
    function PolynomialSpline(order::Integer)
        order >= 1 || throw(ArgumentError("Spline.PolynomialSpline: order must be ≥ 1 (got $order)."))
        return new(order)
    end
end

"""
    Spline.BSpline(d)

A degree-`d` *global* B-spline, used primarily as a basis for least-squares curve *fitting*. A degree-`d`
B-spline produces `(d-1)`th-order-continuous piecewise polynomials, so degree 2/3 resembles a
quadratic/cubic spline. B-splines are *global*: changing one input point perturbs the entire curve (though
it still passes through the other given points).

!!! warning "Not thread-safe for concurrent evaluation"
    A `BSpline`-backed curve is **not safe to evaluate from multiple threads at once**. The underlying
    `DataInterpolations.BSplineInterpolation` reuses a single internal coefficient buffer that it
    overwrites on every evaluation, so concurrent `discount`/`zero`/`forward` calls on one shared curve can
    silently return wrong values. For multithreaded valuation, use a thread-safe interpolant (`Spline.Linear()`,
    `Spline.Quadratic()`, `Spline.Cubic()`, `Spline.PCHIP()`, or `Spline.MonotoneConvex()`), or give each
    thread its own copy of the curve.

For interpolating an already-known curve, prefer the piecewise polynomial convenience constructors
(`Spline.Cubic()` etc.): they build faster and are thread-safe. Knot-rate locality depends on the order.
"""
struct BSpline <: SplineCurve
    order::Int
    function BSpline(order::Integer)
        order >= 1 || throw(ArgumentError("Spline.BSpline: order must be ≥ 1 (got $order)."))
        return new(order)
    end
end

"""
    Spline.PCHIP()

Piecewise Cubic Hermite Interpolating Polynomial (PCHIP). Local and monotonicity-preserving:
each segment depends only on its immediate neighbors, so bumping one rate has bounded effect.
Produces C1-continuous curves (continuous first derivative), giving smooth forward rates
without the non-local coupling of cubic splines.

The default interpolation for `ZeroRateCurve` is `Spline.MonotoneConvex`; PCHIP is a good local,
monotonicity-preserving alternative.
"""
struct PCHIP <: SplineCurve end

"""
    Spline.Akima()

Akima (1970) interpolation. Local and resistant to outlier-induced oscillation:
each segment depends on a few neighboring points. Produces C1-continuous curves.

Compared to PCHIP, Akima can produce slightly different shapes near inflection points.
Both are local; PCHIP additionally preserves monotonicity.
"""
struct Akima <: SplineCurve end

"""
    Spline.MonotoneConvex()

Hagan-West (2006) monotone convex interpolation. With the default `:flat_forward`
extrapolation, it guarantees positive continuous forward rates when input rates imply
positive discrete forwards, and matches
discrete forward rates at knot points. Produces the best KRD locality among smooth methods.

Other extrapolation policies change only the tail beyond the last knot. They can introduce
a forward jump at that knot or negative forwards in the tail; the interior guarantees remain.

Unlike other `SplineCurve` types that wrap DataInterpolations, this dispatches to
`Yield.MonotoneConvex` which implements the Hagan-West sector-based polynomial construction.
Fitting this descriptor returns that native curve for every supported extrapolation policy.

# References
- Hagan & West, "Interpolation Methods for Curve Construction", Applied Mathematical Finance (2006)
"""
struct MonotoneConvex <: SplineCurve end

"""
    Spline.Linear()

Create a local linear spline (returns `PolynomialSpline(1)`, backed by `DataInterpolations.LinearInterpolation`).
This object is not a fitted spline itself, rather it is a placeholder which becomes a spline only after use
within [`fit`](@ref FinanceModels.fit),
or when passed to `ZeroRateCurve`.

Numerically **identical** to `BSpline(1)`, with local knot-rate sensitivities within the knot grid
and thread-safe evaluation (`Spline.BSpline` carries a thread-safety caveat for concurrent evaluation).

# Returns
- A `PolynomialSpline` object representing a linear spline.

# Examples
```julia
julia> Spline.Linear()
PolynomialSpline(1)
```
"""
Linear() = PolynomialSpline(1)

"""
    Spline.Quadratic()

Create a piecewise quadratic spline (returns `PolynomialSpline(2)`, backed by `DataInterpolations.QuadraticSpline`).
This object is not a fitted spline itself, rather it is a placeholder which becomes a spline only after use
within [`fit`](@ref FinanceModels.fit),
or when passed to `ZeroRateCurve`.

Differs numerically from `BSpline(2)` (a global quadratic B-spline); use `Spline.BSpline(2)` to recover the
previous behavior. This piecewise polynomial form is thread-safe.

# Returns
- A `PolynomialSpline` object representing a quadratic spline.

# Examples
```julia
julia> Spline.Quadratic()
PolynomialSpline(2)
```
"""
Quadratic() = PolynomialSpline(2)

"""
    Spline.Cubic()

Create a piecewise natural cubic spline (returns `PolynomialSpline(3)`, backed by `DataInterpolations.CubicSpline`).
This object is not a fitted spline itself, rather it is a placeholder which becomes a spline only after use
within [`fit`](@ref FinanceModels.fit),
or when passed to `ZeroRateCurve`.

Differs numerically from `BSpline(3)` (a global cubic B-spline); use `Spline.BSpline(3)` to recover the
previous behavior. This form builds faster and is thread-safe. Its coefficients depend on the whole
knot grid, so bumping one knot can affect other intervals and the default long-end forward.

# Returns
- A `PolynomialSpline` object representing a cubic spline.

# Examples
```julia
julia> Spline.Cubic()
PolynomialSpline(3)
```
"""
Cubic() = PolynomialSpline(3)


# used as the object which gets optmized before finally returning a completed spline

end
