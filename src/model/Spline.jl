"""
Spline is a module which offers various degree splines used for fitting or bootstraping curves via the [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss) function.

Available methods:

- `Spline.BSpline(n)` where n is the nth order. A spline function of order n is a piecewise polynomial function of degree n − 1. This means that, e.g., cubic polynomial is a fourth degree B-Spline.

This object is not a fitted spline itself, rather it is a placeholder object which will be a spline representing the data only after using within [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss).

and convienience methods which create a `Spline.BSpline` object of the appropriate order.

- `Spline.Linear()`
- `Spline.Quadratic()`
- `Spline.Cubic()`
"""
module Spline
import ..FinanceCore
import ..BSplineKit
import ..AbstractModel


struct BSpline
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
Linear() = BSpline(2)

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
Quadratic() = BSpline(3)

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
Cubic() = BSpline(4)


# used as the object which gets optmized before finally returning a completed spline

end