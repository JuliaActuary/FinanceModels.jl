module Spline
import ..FinanceCore
import ..BSplineKit
import ..AbstractModel


struct BSpline
    order::Int
end

Linear() = BSpline(2)
Quadratic() = BSpline(3)
Cubic() = BSpline(4)


# used as the object which gets optmized before finally returning a completed spline

end