# Least-squares fits of any contract and of any model. Continues least_squares.jl.
using ForwardDiff
using Accessors, IntervalSets

include("helpers/least_squares_setup.jl")   # `lsq_index`, `LsqOpaque`, `zr`, `errors`, `certified`

@testset "least-squares fits: any contract, any model" begin
    qs = lsq_index(200)

    @testset "any contract" begin
        sub = qs[1:60]
        c5 = ZeroRateCurve(fill(0.04, 5), [1.0, 3.0, 5.0, 10.0, 20.0], Spline.Cubic())
        ref = zr(fit(c5, sub, Fit.LeastSquares()))
        # contracts that `fit` cannot look inside
        @test zr(fit(c5, [Quote(q.price, LsqOpaque(q.instrument)) for q in sub], Fit.LeastSquares())) == ref
        # a quote vector of abstract element type, mixing contract types
        mixed = Any[sub..., ZCBPrice(0.97, 0.5), ParYield(0.05, 12.0), Quote(1.9, Composite(Cashflow(1.0, 1.0), Cashflow(1.0, 2.5)))]
        @test certified(errors(Spline.Cubic(), mixed, [1.0, 3.0, 5.0, 10.0, 20.0]), zr(fit(c5, mixed, Fit.LeastSquares())), mixed)
        # option quotes: one volatility for five strikes that imply different ones
        vol(k) = 0.2 + 0.02 * (k - 1.0)
        calls = [Option.EuroCall(CommonEquity(), k, 1.0) for k in 0.8:0.1:1.2]
        oq = [Quote(pv(Equity.BlackScholesMerton(0.03, 0.01, Volatility.Constant(vol(c.strike))), c), c) for c in calls]
        σ = fit(Equity.BlackScholesMerton(0.03, 0.01, Volatility.Constant()), oq, Fit.LeastSquares()).σ.σ
        sse(s) = sum(abs2, pv(Equity.BlackScholesMerton(0.03, 0.01, Volatility.Constant(s)), q.instrument) - q.price for q in oq) / 2
        @test abs(ForwardDiff.derivative(sse, σ)) <= 1.0e-12 * ForwardDiff.derivative(s -> ForwardDiff.derivative(sse, s), σ)
        @test σ ≈ fit(Equity.BlackScholesMerton(0.03, 0.01, Volatility.Constant()), oq).σ.σ atol = 1.0e-6
    end

    @testset "any model" begin
        zq(m) = [Quote(pv(m, Cashflow(1.0, t)), Cashflow(1.0, t)) for t in 0.5:0.5:30]
        # one parameter, four prices
        q4 = ZCBPrice([0.9, 0.8, 0.7, 0.6])
        m = fit(Yield.Constant(), q4, Fit.LeastSquares())
        half(r) = sum(abs2, exp(-r * t) - q.price for (t, q) in zip(1:4, q4)) / 2
        r = m.rate.continuous_value
        @test abs(ForwardDiff.derivative(half, r)) <= 1.0e-12
        @test m isa Yield.Constant
        @test r ≈ fit(Yield.Constant(), q4).rate.continuous_value atol = 1.0e-8
        # Nelson–Siegel recovers the curve that priced the quotes, from the package's default start
        truth = Yield.NelsonSiegel(2.0, 0.058, -0.02, 0.03)
        ns = fit(Yield.NelsonSiegel(), zq(truth), Fit.LeastSquares())
        @test [ns.τ₁, ns.β₀, ns.β₁, ns.β₂] ≈ [2.0, 0.058, -0.02, 0.03] atol = 1.0e-10
        # chosen variables: the decay time stays where the model has it
        betas = (@optic(_.β₀) => -10.0 .. 10.0, @optic(_.β₁) => -10.0 .. 10.0, @optic(_.β₂) => -10.0 .. 10.0)
        nb = fit(Yield.NelsonSiegel(2.0, 0.0, 0.0, 0.0), zq(truth), Fit.LeastSquares(); variables = betas)
        @test nb.τ₁ == 2.0
        @test [nb.β₀, nb.β₁, nb.β₂] ≈ [0.058, -0.02, 0.03] atol = 1.0e-10
        # a short-rate model: at least as good a fit as the loss fit's
        vq = zq(ShortRate.Vasicek(0.3, 0.045, 0.012, Continuous(0.02)))
        v0 = ShortRate.Vasicek(0.1, 0.03, 0.01, Continuous(0.03))
        vsse(v) = sum(abs2, pv(v, q.instrument) - q.price for q in vq)
        @test vsse(fit(v0, vq, Fit.LeastSquares())) <= vsse(fit(v0, vq)) * (1 + 1.0e-8)
        # a minimum beyond a variable's bound throws
        tight = (@optic(_.rate.continuous_value) => 0.0 .. 0.05,)
        @test_throws "upper end of its interval" fit(Yield.Constant(), q4, Fit.LeastSquares(); variables = tight)
        # So does one on the bound, which rounding would otherwise put just inside, with a two-sided
        # derivative: a one-year zero at a price of 1 is fitted by a rate of exactly 0.
        bounded(p) = fit(Yield.Constant(Continuous(0.02)), [ZCBPrice(p, 1.0)], Fit.LeastSquares(); variables = tight).rate.continuous_value
        @test_throws "lower end of its interval" bounded(1.0)
        @test_throws "lower end of its interval" ForwardDiff.derivative(bounded, 1.0)
        near = exp(-1.0e-5)
        @test bounded(near) ≈ 1.0e-5 atol = 1.0e-14
        @test ForwardDiff.derivative(bounded, near) ≈ -1 / near rtol = 1.0e-10
        # This forward requires a spot of exactly 100.001, the upper bound. A quote just
        # inside that bound has a two-sided spot sensitivity of one.
        eurusd = FX.Pair(:EUR, :USD)
        norates = Yield.Constant(Continuous(0.0))
        capped(p) = fit(FX.Forwards(eurusd, 100.0, norates, norates), [Quote(p, FX.Forward(eurusd, 100.0, 1.0))], Fit.LeastSquares(); variables = (@optic(_.spot) => 99.0 .. 100.001,)).spot
        @test_throws "upper end of its interval" capped(100.001 - 100)
        @test_throws "upper end of its interval" ForwardDiff.derivative(capped, 100.001 - 100)
        @test capped(0.00099) ≈ 100.00099 atol = 1.0e-10
        @test ForwardDiff.derivative(capped, 0.00099) ≈ 1 rtol = 1.0e-10
        # Widening the rate interval must preserve the fitted price of the 30-year zero.
        # A full Newton step from a zero rate would overshoot to 54.6.
        wide = (@optic(_.rate.continuous_value) => -1.0e9 .. 1.0e9,)
        @test fit(Yield.Constant(Continuous(0.0)), [ZCBPrice(1.8, 30.0)], Fit.LeastSquares(); variables = wide).rate.continuous_value ≈ -log(1.8) / 30 atol = 1.0e-14
        # Other number types: the solve works in Float64 at least, and keeps a wider type.
        function optimum(prices)
            r = big"0.1"
            for _ in 1:60
                d = [exp(-r * t) for t in 1:3]
                r -= sum(-t * d[t] * (d[t] - prices[t]) for t in 1:3) / sum(t^2 * d[t] * (2 * d[t] - prices[t]) for t in 1:3)
            end
            return r
        end
        for T in (Float32, BigFloat)
            prices = T.([0.91, 0.82, 0.75])
            vars = (@optic(_.rate.continuous_value) => T(-1) .. T(1),)
            fitted = fit(Yield.Constant(Continuous(T(0.05))), ZCBPrice.(prices, T.(1:3)), Fit.LeastSquares(); variables = vars).rate.continuous_value
            @test fitted isa (T == BigFloat ? BigFloat : Float64)
            @test abs(fitted - optimum(prices)) <= 1.0e-12
        end
        for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            ts = [1.0, 2.0, 4.0, 6.0]
            curve = ZeroRateCurve(fill(big"0.04", 4), ts, spline)
            prices = [exp(-t * (big"0.03" + big"0.001" * t)) for t in ts]
            first_rate(p) = first(zr(fit(curve, ZCBPrice.([p; prices[2:end]], ts), Fit.LeastSquares())))
            @test first_rate(prices[1]) isa BigFloat
            sensitivity = ForwardDiff.derivative(first_rate, prices[1])
            @test sensitivity isa BigFloat
            @test sensitivity ≈ -1 / prices[1] rtol = 1.0e-10
        end
        # More parameters than one ForwardDiff chunk, including the last partial chunk.
        ts17 = collect(1.0:17.0)
        prices17 = exp.(-0.03 .* ts17)
        long_curve = ZeroRateCurve(fill(0.04, 17), ts17, Spline.Linear())
        last_rate(p) = last(zr(fit(long_curve, ZCBPrice.([prices17[1:(end - 1)]; p], ts17), Fit.LeastSquares())))
        @test ForwardDiff.derivative(last_rate, last(prices17)) ≈ -1 / (17 * last(prices17)) rtol = 1.0e-10
        # A discrete set excludes the continuous optimum of 0.03 between its first two values.
        three = [ZCBPrice(exp(-0.03), 1.0)]
        @test_throws "variable 1 is not given an interval" fit(Yield.Constant(Continuous(0.05)), three, Fit.LeastSquares(); variables = (@optic(_.rate.continuous_value) => [0.0, 0.05, 0.1],))
        @test fit(Yield.Constant(Continuous(0.05)), three, Fit.LeastSquares(); variables = (@optic(_.rate.continuous_value) => 0.0 .. 0.1,)).rate.continuous_value ≈ 0.03 atol = 1.0e-12
        # The certificate measures a correction by the prices it moves, not in the variable's units: the
        # same contract in years and with cashflows 1e11 years out, where the rate that fits is 7e-12.
        # The two cashflows are worth d² - 3d for d = discount(t), and -1.23 at d = 0.49005.
        for t in (1.0, 1.0e11)
            far = [Quote(-1.23, Composite(Cashflow(-3.0, t), Cashflow(1.0, 2t)))]
            d(p) = discount(fit(Yield.Constant(Continuous(0.0)), [Quote(p, far[1].instrument)], Fit.LeastSquares()), t)
            @test d(-1.23) ≈ 0.4900495061637922 rtol = 1.0e-10
            @test ForwardDiff.derivative(d, -1.23) ≈ 1 / (2 * 0.4900495061637922 - 3) rtol = 1.0e-8
        end
        # What the method needs of a model: `variables`, and prices that ForwardDiff differentiates twice.
        @test_throws MethodError fit(Yield.SmithWilson(ufr = 0.03, α = 0.1), ZCBPrice([0.97, 0.94]), Fit.LeastSquares())
        hw = ShortRate.HullWhite(0.1, 0.015, Yield.Constant(Continuous(0.03)))
        swaptions = [Option.Swaption(e, e + 5.0, 0.03, 2) for e in 1.0:3.0]
        @test_throws "first-order" fit(ShortRate.HullWhite(0.05, 0.01, Yield.Constant(Continuous(0.03))), [Quote(pv(hw, s), s) for s in swaptions], Fit.LeastSquares())
        # more variables than the quotes determine: these prices came from a Nelson–Siegel curve
        @test_throws FitConvergenceError fit(Yield.NelsonSiegelSvensson(), zq(truth), Fit.LeastSquares())
    end
end
