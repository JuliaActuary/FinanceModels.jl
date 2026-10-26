@testset "exact spline fits" begin
    # the largest pricing residual, relative to each quote's price
    reprice(c, qs) = maximum(abs(present_value(c, q.instrument) - q.price) / abs(q.price) for q in qs)
    cmt = CMTYield.(
        [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.4, 1.74, 2.31, 2.41] ./ 100,
        [1 / 12, 2 / 12, 3 / 12, 6 / 12, 1, 2, 3, 5, 7, 10, 20, 30],
    )
    zcb = ZCBPrice.([0.99, 0.97, 0.94, 0.9], [1.0, 2.0, 3.0, 5.0])
    mixed = [ZCBYield(0.03, 1.0); CMTYield.([0.032, 0.035], [2, 5]); OISYield(0.036, 7.0); ParSwapYield(0.037, 10; frequency = 2)]
    splines = (
        Spline.Linear(), Spline.Quadratic(), Spline.Cubic(), Spline.BSpline(2), Spline.BSpline(3),
        Spline.PCHIP(), Spline.Akima(), Spline.MonotoneConvex(),
    )

    @testset "$spline reprices its quotes" for spline in splines
        for qs in (cmt, zcb, mixed)
            c = fit(spline, qs, Fit.Exact())
            @test reprice(c, qs) < 1.0e-12
            @test knot_tenors(c) == sort(float.(maturity.(qs)))
        end
        # a loss fit converges to the same curve
        c, cl = fit(spline, zcb, Fit.Exact()), fit(spline, zcb, Fit.Loss(abs2); solve_kwargs = (; g_tol = 1.0e-12))
        @test rate.(knot_rates(c)) ≈ rate.(knot_rates(cl)) atol = 1.0e-7
        for extrapolation in (:flat_zero, :linear)
            @test reprice(fit(spline, mixed, Fit.Exact(); extrapolation), mixed) < 1.0e-12
        end
    end

    @testset "a linear exact fit is the bootstrap" begin
        for qs in (cmt, zcb, mixed)
            @test rate.(knot_rates(fit(Spline.Linear(), qs, Fit.Exact()))) ≈ rate.(knot_rates(fit(Spline.Linear(), qs, Fit.Bootstrap()))) atol = 1.0e-13
        end
    end

    @testset "notionals do not matter" begin
        unit = [Quote(p, Cashflow(1.0, t)) for (p, t) in zip([0.99, 0.97, 0.94, 0.9], [1.0, 2.0, 3.0, 5.0])]
        large = [Quote(1.0e9 * q.price, Cashflow(1.0e9, q.instrument.time)) for q in unit]
        for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            @test rate.(knot_rates(fit(spline, large, Fit.Exact()))) ≈ rate.(knot_rates(fit(spline, unit, Fit.Exact()))) rtol = 1.0e-13
        end
    end

    @testset "a refit on the quote maturities is the spline fit" begin
        tenors = float.(maturity.(cmt))
        for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            c0 = ZeroRateCurve(fill(0.02, length(tenors)), tenors, spline)
            @test isequal(fit(c0, cmt), fit(spline, cmt))
            @test isequal(fit(c0, cmt), fit(c0, cmt, Fit.Loss()))
            @test isequal(fit(spline, cmt), fit(spline, cmt, Fit.Loss()))
            @test isequal(fit(c0, cmt, Fit.Exact()), fit(spline, cmt, Fit.Exact()))
        end
        # other knot grids minimize the squared loss, as before
        c0 = ZeroRateCurve(fill(0.02, 4), [0.5, 2.0, 10.0, 30.0], Spline.Linear())
        @test isequal(fit(c0, cmt), fit(c0, cmt, FinanceModels.__default_loss(c0)))
    end

    @testset "an optimizer or its settings fit the squared loss" begin
        c0 = ZeroRateCurve(fill(0.02, 4), [1.0, 2.0, 3.0, 5.0], Spline.Cubic())
        # without a method, they select Fit.Loss(), bit for bit, and reach the optimizer
        @test isequal(
            fit(Spline.Cubic(), zcb; solve_kwargs = (; g_tol = 1.0e-12)),
            fit(Spline.Cubic(), zcb, Fit.Loss(); solve_kwargs = (; g_tol = 1.0e-12)),
        )
        @test isequal(
            fit(c0, zcb; solve_kwargs = (; g_tol = 1.0e-12)),
            fit(c0, zcb, Fit.Loss(); solve_kwargs = (; g_tol = 1.0e-12)),
        )
        @test_throws FitConvergenceError fit(Spline.Cubic(), zcb; solve_kwargs = (; maxiters = 1))
        @test_throws FitConvergenceError fit(Spline.Cubic(), zcb; optimizer = FailingFitOptimizer())
        @test_throws FitConvergenceError fit(c0, zcb; solve_kwargs = (; maxiters = 1))
        @test_throws FitConvergenceError fit(Spline.Cubic(), zcb, Fit.Loss(); solve_kwargs = (; maxiters = 1))
        # Fit.Loss() is the squared loss, bit for bit
        @test isequal(
            fit(Spline.Cubic(), zcb, Fit.Loss(); solve_kwargs = (; g_tol = 1.0e-12)),
            fit(Spline.Cubic(), zcb, Fit.Loss(x -> x^2); solve_kwargs = (; g_tol = 1.0e-12)),
        )
        @test isequal(fit(Yield.NelsonSiegel(), cmt), fit(Yield.NelsonSiegel(), cmt, Fit.Loss(x -> x^2)))
    end

    @testset "Fit.Exact" begin
        for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            @test reprice(fit(spline, cmt, Fit.Exact()), cmt) < 1.0e-12
        end
        @test_throws MethodError Fit.Exact(10)        # keyword only
        @test_throws FitConvergenceError fit(Spline.Cubic(), cmt, Fit.Exact(maxiters = 1))
        # a refit with Fit.Exact needs one quote per knot; its knots need not be the maturities
        @test_throws "one quote per knot" fit(ZeroRateCurve(fill(0.02, 3), [1.0, 5.0, 10.0], Spline.Linear()), zcb, Fit.Exact())
        c4 = fit(ZeroRateCurve(fill(0.02, 4), [1.5, 2.5, 3.5, 5.0], Spline.Linear()), zcb, Fit.Exact())
        @test reprice(c4, zcb) < 1.0e-12
    end

    @testset "residuals are measured against each quote's size on the fitted curve" begin
        # a long zero-coupon price far below its value at the seed reprices to rounding
        tiny = [ZCBPrice(1.0e-12, 50.0)]
        for f in (() -> fit(Spline.Linear(), tiny, Fit.Exact()), () -> fit(Spline.Linear(), tiny, Fit.Bootstrap()))
            @test reprice(f(), tiny) < 1.0e-12
        end
        # a negative price, which no discount factor reprices, throws
        @test_throws FitConvergenceError fit(Spline.Linear(), [ZCBPrice(-1.0e-20, 50.0)], Fit.Exact())
        @test_throws ArgumentError fit(Spline.Linear(), [ZCBPrice(-1.0e-20, 50.0)], Fit.Bootstrap())
        # a contract of fixed cashflows worth almost nothing is measured against its cashflows, not its price
        spread = [ZCBPrice(0.98, 1.0), Quote(1.0e-10, FinanceCore.Composite(Cashflow(1.0, 1.0), Cashflow(-1.01, 2.0)))]
        c = fit(Spline.Linear(), spread, Fit.Exact())
        @test abs(pv(c, spread[2].instrument) - spread[2].price) < 1.0e-14
    end

    @testset "quotes that no curve reprices throw" begin
        # a negative price needs a negative discount factor
        bad = ZCBPrice.([0.97, -0.1], [1.0, 2.0])
        for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            err = try
                fit(spline, bad, Fit.Exact())
                nothing
            catch e
                e
            end
            @test err isa FitConvergenceError
            @test occursin("Fit.Loss()", sprint(showerror, err))
        end
    end

    @testset "derivatives through an exact fit" begin
        tenors = [1.0, 2.0, 5.0, 10.0]
        rates = [0.03, 0.032, 0.035, 0.037]
        cfs = Cashflow.([3.0, 3.0, 103.0], [2.0, 5.0, 10.0])
        value(spline) = x -> pv(fit(spline, OISYield.(x, tenors), Fit.Exact()), cfs)
        @test ForwardDiff.gradient(value(Spline.Linear()), rates) ≈
            ForwardDiff.gradient(x -> pv(fit(Spline.Linear(), OISYield.(x, tenors), Fit.Bootstrap()), cfs), rates) rtol = 1.0e-12
        h = 1.0e-6
        for spline in (Spline.Cubic(), Spline.MonotoneConvex())
            g = ForwardDiff.gradient(value(spline), rates)
            fd = [(value(spline)(rates .+ h .* (1:4 .== i)) - value(spline)(rates .- h .* (1:4 .== i))) / 2h for i in 1:4]
            @test g ≈ fd rtol = 1.0e-6
        end
        # a dual number hidden in a contract the fit does not strip reaches the solve and throws
        @test_throws "reached the solve through data it does not differentiate" ForwardDiff.derivative(
            a -> discount(fit(Spline.Linear(), [Quote(0.97, WrappedContract(Cashflow(a, 1.0))), ZCBPrice(0.94, 2.0)], Fit.Exact()), 1.5), 1.0
        )
    end

    @testset "FX forwards fit their spline foreign curve exactly" begin
        eurusd = FX.Pair(:EUR, :USD)
        usd = fit(Spline.Linear(), ZCBYield.([0.04, 0.045, 0.05, 0.052], [1.0, 2.0, 5.0, 10.0]), Fit.Bootstrap())
        S = 1.1
        truth = FX.Forwards(eurusd, S, usd, Yield.Constant(Continuous(0.025)))
        quotes = [FX.Outright(eurusd, forward(truth, t), t) for t in [0.5, 1.0, 2.0, 5.0]]
        for spline in (Spline.Linear(), Spline.Cubic())
            m = fit(FX.Forwards(eurusd, S, usd, spline), quotes, Fit.Exact())
            @test isequal(fit(FX.Forwards(eurusd, S, usd, spline), quotes), fit(FX.Forwards(eurusd, S, usd, spline), quotes, Fit.Loss()))
            @test all(abs(pv(m, q.instrument)) < 1.0e-12 for q in quotes)
            @test isequal(m.foreign, fit(spline, map(q -> FinanceModels.FX.__implied_foreign_quote(truth, q), quotes), Fit.Exact()))
        end
    end
end
