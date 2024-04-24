"""
Spline is a module which offers various degree splines used for fitting or bootstraping curves via the [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss) function.

Available methods:

- `Spline.BSpline(n)` where n is the nth order. A nth-order B-Spline is analagous to an (n-1)th order polynomial spline. That is, a 3rd/4th order BSpline is very similar to a quadratic/cubic spline respectively. BSplines are global in that a change in one point affects the entire spline (though the spline still passes through the other given points still).
- `Spline.PolynomialSpline(n)` where n is the nth order.

This object is not a fitted spline itself, rather it is a placeholder object which will be a spline representing the data only after using within [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss).

and convienience methods which create a `Spline.BSpline` object of the appropriate order.

- `Spline.Linear()` equals `BSpline(2)`
- `Spline.Quadratic()` equals `BSpline(3)`
- `Spline.Cubic()` equals `BSpline(4)`
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