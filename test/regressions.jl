# Test-only solver that deterministically exercises `fit`'s unsuccessful-retcode path.
struct FailingFitOptimizer end
function FinanceModels.Optimization.solve(
        prob::FinanceModels.Optimization.OptimizationProblem,
        ::FailingFitOptimizer
    )
    return (
        u = prob.u0,
        retcode = FinanceModels.Optimization.SciMLBase.ReturnCode.Failure,
    )
end

# Regression tests from the 2026-06 ecosystem audit
@testset "audit regressions" begin
    @testset "fit(spline, quotes, Fit.Loss) uses the supplied loss" begin
        qs = ZCBPrice([0.9, 0.8, 0.7])
        calls = Ref(0)
        counted = Fit.Loss(x -> (calls[] += 1; x^2))
        fit(Spline.Linear(), qs, counted)
        # the user-supplied loss was previously silently replaced with the default
        @test calls[] > 0
    end

    @testset "par with a non-representable stub maturity" begin
        c = Yield.Constant(0.04)
        @test rate(par(c, 4)) ≈ 0.03960780543711406
        @test par(c, 0.2; frequency = 4).compounding == Periodic(5)
        # 1/0.3 is not an integer → informative error instead of an InexactError
        @test_throws ArgumentError par(c, 0.3)
    end

    @testset "present_value past maturity is 0, not an error" begin
        b = Bond.Fixed(0.05, Periodic(1), 3)
        m = Yield.Constant(0.03)
        @test present_value(m, b, 4.0) == 0.0
        @test present_value(m, Projection(b, m, CashflowProjection()), 5.0) == 0.0
        # (n.b. `present_value(m, ::Cashflow, t)` hits FinanceCore's Cashflow method,
        # which discounts to time 0 regardless of the third argument — so the
        # empty-fold path is exercised via a single-cashflow Projection instead)
        @test present_value(m, Projection(Cashflow(10.0, 1.0), m, CashflowProjection()), 2.0) == 0.0
    end

    @testset "Forward contract shifts cashflow times" begin
        b = Bond.Fixed(0.05, Periodic(1), 2)
        fwd = Forward(1.0, b)
        cfs = collect(Projection(fwd, NullModel(), CashflowProjection()))
        @test [cf.time for cf in cfs] ≈ [2.0, 3.0]
        @test [cf.amount for cf in cfs] ≈ [0.05, 1.05]
        # for a flat curve, shifting all cashflows by Δ scales the PV by discount(Δ)
        m = Yield.Constant(0.03)
        @test present_value(m, fwd) ≈ present_value(m, b) * discount(m, 1.0)
    end

    @testset "ParSwapYield" begin
        # Swap fixed-leg conventions differ by market, so the frequency is required.
        @test_throws UndefKeywordError ParSwapYield(0.04, 5)
        q = ParSwapYield(0.04, 5; frequency = 4)
        @test q.price ≈ 1.0
        @test q.instrument.frequency == Periodic(4)
        @test ParSwapYield(0.04, 5; frequency = Periodic(1)).instrument.frequency == Periodic(1)
        # a Periodic Rate input must agree with an explicit frequency
        q2 = ParSwapYield(Periodic(0.04, 2), 5; frequency = 2)
        @test q2.instrument.frequency == Periodic(2)
        @test_throws ArgumentError ParSwapYield(Periodic(0.04, 2), 5; frequency = 4)
        # round-trip: a curve fit to par-swap quotes reprices them to root-finder
        # precision (the fitting-time curve and the returned curve are identical
        # for local interpolants, including the pinned t=0 knot)
        swap_rates = [0.02, 0.025, 0.03]
        for frequency in (1, 4)
            qs = ParSwapYield.(swap_rates, [1, 2, 3]; frequency)
            c = fit(Spline.Linear(), qs, Fit.Bootstrap())
            for (r, t) in zip(swap_rates, [1, 2, 3])
                @test rate(par(c, t; frequency)) ≈ r atol = 1.0e-10
            end
        end
        @test_throws UndefKeywordError ParSwapYield(swap_rates)
        @test ParSwapYield(swap_rates; frequency = 1) == ParSwapYield.(swap_rates, [1.0, 2.0, 3.0]; frequency = 1)
    end

    @testset "ParYield frequency" begin
        @test ParYield(0.04, 5).instrument.frequency == Periodic(2)
        @test ParYield(0.04, 5; frequency = 1).instrument.frequency == Periodic(1)
        # A Periodic rate sets its own frequency; a conflicting frequency is an error,
        # not a silently ignored keyword.
        @test ParYield(Periodic(0.04, 4), 5).instrument.frequency == Periodic(4)
        @test ParYield(Periodic(0.04, 4), 5; frequency = 4) == ParYield(Periodic(0.04, 4), 5)
        @test_throws ArgumentError ParYield(Periodic(0.04, 4), 5; frequency = 2)
        # Other rates convert to the requested frequency.
        @test ParYield(Continuous(0.04), 5; frequency = 1).instrument.coupon_rate ≈ exp(0.04) - 1
    end

    @testset "OISYield conventions" begin
        # One year or less: a single payment. Longer: annual payments on both legs.
        @test OISYield(0.04, 0.5).instrument == Bond.Fixed(0.0, Periodic(1), 0.5)
        @test OISYield(0.04, 0.5).price ≈ 1.04^-0.5
        # annual coupons at the quoted rate (a FinanceCore older than 2.6 converts the rate to
        # itself through continuous compounding, off by an ulp)
        b = OISYield(0.04, 5).instrument
        @test b isa Bond.Fixed && b.frequency == Periodic(1) && b.maturity == 5
        @test b.coupon_rate ≈ 0.04 rtol = 1.0e-15
        @test OISYield(0.04, 5).price == 1.0
        # A one-year single payment equals a one-coupon annual par swap.
        @test pv(Yield.Constant(0.04), OISYield(0.04, 1).instrument) ≈ OISYield(0.04, 1).price
        @test pv(Yield.Constant(Periodic(0.04, 1)), OISYield(0.04, 5).instrument) ≈ 1.0 atol = 1.0e-14
    end

    @testset "PCHIP and Akima ZeroRateCurve" begin
        zrates = [0.02, 0.025, 0.03, 0.035]
        tenors = [1.0, 2.0, 5.0, 10.0]
        @testset "$spline" for spline in (Spline.PCHIP(), Spline.Akima())
            zrc = ZeroRateCurve(zrates, tenors, spline)
            for (r, t) in zip(zrates, tenors)
                @test rate(zero(zrc, t)) ≈ r atol = 1.0e-10
            end
            @test 0 < discount(zrc, 3.0) < 1
        end
    end

    @testset "PCHIP and Akima loss fits do not stall at a flat seed" begin
        tenors = [1 / 12, 2 / 12, 3 / 12, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
        rates = [
            0.0375, 0.0372, 0.0372, 0.0372, 0.0364,
            0.0368, 0.0369, 0.038, 0.04, 0.0423,
        ]

        for quote_type in (CMTYield, ParYield), spline in (Spline.PCHIP(), Spline.Akima())
            qs = quote_type.(rates, tenors)
            curve = fit(spline, qs)
            max_price_error = maximum(
                abs,
                present_value(curve, q.instrument) - q.price for q in qs
            )
            @test max_price_error < 1.0e-6
        end
    end

    @testset "optimizer failures are never returned as fitted models" begin
        qs = ZCBPrice.([0.97, 0.93, 0.88], [1.0, 2.0, 3.0])
        fits = (
            () -> fit(Yield.Constant(), qs; optimizer = FailingFitOptimizer()),
            () -> fit(Spline.MonotoneConvex(), qs; optimizer = FailingFitOptimizer()),
            () -> fit(Spline.Linear(), qs, Fit.Loss(abs2); optimizer = FailingFitOptimizer()),
        )

        for run_fit in fits
            err = try
                run_fit()
                nothing
            catch e
                e
            end
            @test err isa FitConvergenceError
            @test err isa Exception
            if err isa FitConvergenceError
                @test err.retcode == FinanceModels.Optimization.SciMLBase.ReturnCode.Failure
                msg = sprint(showerror, err)
                @test startswith(msg, "FitConvergenceError: ")
                @test occursin("optimizer return code Failure", msg)
            end
        end
    end

    @testset "ZeroRateCurve negative time" begin
        zrc = ZeroRateCurve([0.02, 0.03], [1.0, 2.0])
        @test discount(zrc, 0.0) == 1.0
        @test_throws DomainError discount(zrc, -1.0)
    end

    @testset "bootstrap accepts unsorted quotes" begin
        sorted = fit(Spline.Linear(), ZCBPrice.([0.9, 0.8, 0.7], [1.0, 2.0, 3.0]), Fit.Bootstrap())
        unsorted = fit(Spline.Linear(), ZCBPrice.([0.8, 0.9, 0.7], [2.0, 1.0, 3.0]), Fit.Bootstrap())
        for t in 0.5:0.5:3.0
            @test discount(sorted, t) ≈ discount(unsorted, t)
        end
    end

    @testset "bootstrap rejects duplicate maturities loudly" begin
        # previously surfaced as a cryptic root-bracketing failure
        @test_throws ArgumentError fit(Spline.Linear(), ZCBPrice.([0.9, 0.89, 0.7], [1.0, 1.0, 3.0]), Fit.Bootstrap())
    end

    @testset "MonotoneConvex fit with a single quote" begin
        # range(0.01, 0.05, length=1) used to throw before any solving happened
        c = fit(Spline.MonotoneConvex(), [ZCBPrice(0.95, 2.0)])
        @test discount(c, 2.0) ≈ 0.95 atol = 1.0e-8
    end

    @testset "second-order optimizer fit does not warn about ADtype" begin
        # The generic and MonotoneConvex loss functions declared first-order
        # `AutoForwardDiff`, so a second-order optimizer made OptimizationBase warn and
        # auto-promote to `SecondOrder` once per fit. They now declare `SecondOrder` AD.
        # The bare `@test_logs` (default `min_level = Warn`) asserts no warnings; the
        # default `LBFGS` path only requests gradients, so it is unaffected.
        qs = ZCBPrice.([0.97, 0.93, 0.88, 0.8], [1.0, 2.0, 3.0, 5.0])
        reprice(c) = maximum(abs, present_value(c, q.instrument) - q.price for q in qs)

        # generic (bounded) path: IPNewton is the bounds-compatible second-order method
        # (plain Newton is rejected by Optim's Fminbox for box-constrained problems).
        ipn = FinanceModels.OptimizationOptimJL.IPNewton()
        cp = @test_logs fit(Yield.CairnsPritchard(), qs; optimizer = ipn)
        @test reprice(cp) < 1.0e-8
        @test_logs fit(Yield.NelsonSiegel(), qs; optimizer = ipn)  # smooth 4-param fit: assert only quiet

        # MonotoneConvex (unbounded) path: Newton is the canonical second-order method.
        mc = @test_logs fit(Spline.MonotoneConvex(), qs; optimizer = FinanceModels.OptimizationOptimJL.Newton())
        @test reprice(mc) < 1.0e-10
    end
end

@testset "batch knot-rate optic rebuilds once per candidate" begin
    # One `@optic(_.rates[i])` per knot rebuilt the whole curve once per knot for every
    # optimizer candidate, so allocation grew with the square of the knot count. The
    # batch optic rebuilds once, so allocation should grow roughly linearly.
    setall_bytes(c, vals) = @allocated Accessors.setall(c, FinanceModels.KnotRatesOptic(), vals)
    builders = (
        n -> Yield.MonotoneConvex(collect(range(0.02, 0.05; length = n)), collect(1.0:n)),
        n -> ZeroRateCurve(collect(range(0.02, 0.05; length = n)), collect(1.0:n), Spline.Linear()),
    )
    for build in builders
        bytes = map((5, 50, 500)) do n
            c = build(n)
            vals = collect(c.rates) .+ 0.001
            setall_bytes(c, vals)   # compile
            setall_bytes(c, vals)
        end
        # Linear growth gives ratios of at most about 10 per tenfold increase in knots;
        # per-knot rebuilds gave ratios of about 30 (5 → 50) and 80 (50 → 500).
        @test bytes[2] / bytes[1] < 20
        @test bytes[3] / bytes[2] < 20
    end
end
