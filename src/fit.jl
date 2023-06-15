module Fit
abstract type FitMethod end

struct Loss{T} <: FitMethod
    fn::T
end

struct Bootstrap <: FitMethod
    # spline method
end


end


__default_optic(m::Yield.Constant) = @optic(_.rate.value) => -1.0 .. 1.0
__default_optic(m::Equity.BlackScholesMerton) = __default_optic(m.σ)
__default_optic(m::Volatility.Constant) = @optic(_.σ) => -0.0 .. 10.0

function fit(mod0, quotes, method::F=Fit.Loss(x -> x^2);
    variables=OptArgs(__default_optic(mod0))
) where
{F<:Fit.Loss}
    # find the rate that minimizes the loss function w.r.t. the calculated price vs the quotes
    loss(m, quotes) =
        mapreduce(+, quotes) do q
            method.fn(pv(m, q.instrument) - q.price)
        end

    f = Base.Fix2(loss, quotes)
    ops = OptProblemSpec(f, SVector, mod0, variables)
    sol = solve(ops, ECA(), maxiters=300)

    return sol.uobj

end