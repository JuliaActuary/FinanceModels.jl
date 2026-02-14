"""
Spline is a module which offers various degree splines used for fitting or bootstraping curves via the [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss) function.

Available methods:

- `Spline.BSpline(n)` where n is the nth order. A nth-order B-Spline is analogous to an (n-1)th order polynomial spline. That is, a 3rd/4th order BSpline is very similar to a quadratic/cubic spline respectively. BSplines are global in that a change in one point affects the entire spline (though the spline still passes through the other given points still).
- `Spline.PolynomialSpline(n)` where n is the nth order.

This object is not a fitted spline itself, rather it is a placeholder object which will be a spline representing the data only after using within [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss).

Convenience methods which create a `Spline.BSpline` object of the appropriate order:

- `Spline.Linear()` equals `BSpline(2)`
- `Spline.Quadratic()` equals `BSpline(3)`
- `Spline.Cubic()` equals `BSpline(4)`

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

struct PolynomialSpline <: SplineCurve
    order::Int
end

struct BSpline <: SplineCurve
    order::Int
end

"""
    Spline.PCHIP()

Piecewise Cubic Hermite Interpolating Polynomial (PCHIP). Local and monotonicity-preserving:
each segment depends only on its immediate neighbors, so bumping one rate has bounded effect.
Produces C1-continuous curves (continuous first derivative), giving smooth forward rates
without the non-local coupling of cubic splines.

PCHIP is the default interpolation for `ZeroRateCurve`.
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

Hagan-West (2006) monotone convex interpolation. Finance-aware: guarantees positive
continuous forward rates (when input rates imply positive forwards) and matches
discrete forward rates at knot points. Produces the best KRD locality among smooth methods.

Unlike other `SplineCurve` types that wrap DataInterpolations, this dispatches to
`Yield.MonotoneConvex` which implements the Hagan-West sector-based polynomial construction.

# References
- Hagan & West, "Interpolation Methods for Curve Construction", Applied Mathematical Finance (2006)
"""
struct MonotoneConvex <: SplineCurve end

"""
    Spline.Linear()

Create a linear B-spline. This object is not a fitted spline itself, rather it is a placeholder object which will be a spline representing the data only after using within [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss).

# Returns
- A `BSpline` object representing a linear B-spline.



# Examples
```julia
julia> Spline.Linear()
BSpline(2)
```
"""
Linear() = BSpline(1)

"""
    Spline.Quadratic()

Create a quadratic B-spline. This object is not a fitted spline itself, rather it is a placeholder object which will be a spline representing the data only after using within [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss).

# Returns
- A `BSpline` object representing a quadratic B-spline.

# Examples
```julia
julia> Spline.Quadratic()
BSpline(3)
```
"""
Quadratic() = BSpline(2)

"""
    Spline.Cubic()

Create a cubic B-spline. This object is not a fitted spline itself, rather it is a placeholder object which will be a spline representing the data only after using within [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss).

# Returns
- A `BSpline` object representing a cubic B-spline.

# Examples
```julia
julia> Spline.Cubic()
BSpline(4)
```
"""
Cubic() = BSpline(3)


# used as the object which gets optmized before finally returning a completed spline

end
