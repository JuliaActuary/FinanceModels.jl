module Fit

abstract type FitMethod end

struct Loss{T} <: FitMethod
    fn::T
end

struct Bootstrap <: FitMethod
    # spline method
end


end


__default_optic(m::Yield.Constant) = OptArgs(@optic(_.rate.value) => -1.0 .. 1.0)
__default_optic(m::Yield.IntermediateYieldCurve) = OptArgs(@optic(_.ys[end]) => 0.0 .. 1.0)
__default_optic(m::Yield.NelsonSiegel) = OptArgs([
    @optic(_.τ₁) => 0.0 .. 100.0
    @optic(_.β₀) => -10.0 .. 10.0
    @optic(_.β₁) => -10.0 .. 10.0
    @optic(_.β₂) => -10.0 .. 10.0
]...)
__default_optic(m::Yield.NelsonSiegelSvensson) = OptArgs([
    @optic(_.τ₁) => 0.0 .. 100.0,
    @optic(_.τ₂) => 0.0 .. 100.0,
    @optic(_.β₀) => -10.0 .. 10.0,
    @optic(_.β₁) => -10.0 .. 10.0,
    @optic(_.β₂) => -10.0 .. 10.0,
    @optic(_.β₃) => -10.0 .. 10.0,
]...)
__default_optic(m::Equity.BlackScholesMerton) = __default_optic(m.σ)
__default_optic(m::Volatility.Constant) = OptArgs(@optic(_.σ) => -0.0 .. 10.0)


__default_optim(m) = ECA()

function fit(mod0, quotes, method::F=Fit.Loss(x -> x^2);
    variables=__default_optic(mod0),
    optimizer=__default_optim(mod0)
) where
{F<:Fit.Loss}
    # find the rate that minimizes the loss function w.r.t. the calculated price vs the quotes
    f = __loss_single_function(method, quotes)
    ops = OptProblemSpec(f, Vector, mod0, variables)
    sol = solve(ops, optimizer)
    return sol.uobj

end

function fit(mod0::Spline.BSpline, quotes, method::Fit.Bootstrap)
    discount_vector = [0.0]
    times = [maturity(quotes[1])]

    discount_vector[1] = let
        m = fit(Yield.Constant(), [quotes[1]], Fit.Loss(x -> x^2))
        discount(m, times[1])
    end

    for i in eachindex(quotes)[2:end]
        q = quotes[i]
        push!(times, maturity(q))
        push!(discount_vector, 0.0)
        m = Yield.IntermediateYieldCurve(mod0, times, discount_vector)
        discount_vector[i] = let
            m = fit(m, [q], Fit.Loss(x -> x^2))
            discount(m, times[i])
        end

    end
    zero_vec = -log.(clamp.(discount_vector, 0.00001, 1)) ./ times
    return Yield.Spline(mod0, [zero(eltype(times)); times], [first(zero_vec); zero_vec])
    # return Yield.Spline(mod0, times, zero_vec)

end

function fit(mod0::Yield.SmithWilson, quotes)
    cm, ts = cashflows_timepoints(quotes)
    prices = [q.price for q in quotes]

    return Yield.SmithWilson(ts, cm, prices; ufr=mod0.ufr, α=mod0.α)

end

function __loss_single_function(loss_method, quotes)
    function loss(m, quotes)
        return mapreduce(+, quotes) do q
            loss_method.fn(present_value(m, q.instrument) - q.price)
        end
    end
    return Base.Fix2(loss, quotes) # a function that takes a model and returns the loss
end