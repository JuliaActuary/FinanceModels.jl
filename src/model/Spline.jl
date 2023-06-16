module Spline
import ..BSplineKit
abstract type SplineType end

struct BSpline <: SplineType
    order::Int
end

Linear() = BSpline(2)
Quadratic() = BSpline(3)
Cubic() = BSpline(4)


# struct Curve{F,U,V}
#     fn::F
#     xs::Vector{U}
#     ys::Vector{V}
# end

function Curve(b::BSpline, xs, ys)
    @show order = min(length(xs), b.order) # in case the length of xs is less than the spline order
    @show xs, ys
    int = BSplineKit.interpolate(xs, ys, BSplineKit.BSplineOrder(order))
    return BSplineKit.extrapolate(int, BSplineKit.Smooth())
end

end