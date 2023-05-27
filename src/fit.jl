
# using Optimization, OptimizationOptimJL
# rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
# cons = (res, x, p) -> res .= [x[1]^2 + x[2]^2]
# x0 = zeros(2)
# p = [1.0, 100.0]
# prob = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff(); cons = cons)
# prob = Optimization.OptimizationProblem(prob, x0, p, lcons = [-5.0], ucons = [10.0])
# sol = solve(prob, IPNewton())

module Fit 
    abstract type FitMethod end

    struct Loss{T} <: FitMethod
        fn::T
    end

    struct Bootstrap <: FitMethod
        # spline method
    end


end


function fit(::Type{Yield.Constant},method::Fit.Loss,quotes)
    # find the rate that minimizes the loss function w.r.t. the calculated price vs the quotes
    function outer(x,p)
        m = Yield.Constant(x[1])
        mapreduce(+,quotes) do q
            method.fn(pv(m,q.instrument) - q.price)
        end
    end

    x0 = [0.01]

    prob = OptimizationFunction(outer, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(prob, x0, [0.])
    sol = solve(prob, Newton())
    return Rate(only(sol))



end
