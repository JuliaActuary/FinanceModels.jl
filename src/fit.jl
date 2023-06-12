module Fit
abstract type FitMethod end

struct Loss{T} <: FitMethod
    fn::T
end

struct Bootstrap <: FitMethod
    # spline method
end


end


function fit(::Type{Yield.Constant}, method::Fit.Loss, quotes)
    # find the rate that minimizes the loss function w.r.t. the calculated price vs the quotes
    function outer(x, p)
        m = Yield.Constant(x[1])
        mapreduce(+, quotes) do q
            method.fn(value(m, q.instrument) - q.price)
        end
    end

    x0 = [0.01]

    prob = OptimizationFunction(outer, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(prob, x0, [0.0])
    sol = solve(prob, Newton())
    return Yield.Constant(only(sol))

end
