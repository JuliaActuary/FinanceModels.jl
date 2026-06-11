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
        q = ParSwapYield(0.04, 5)
        @test q.price ≈ 1.0
        @test q.instrument.frequency == Periodic(4) # quarterly by default
        # a Periodic Rate input carries its own frequency (overrides the kwarg default)
        q2 = ParSwapYield(Periodic(0.04, 2), 5)
        @test q2.instrument.frequency == Periodic(2)
        # round-trip: a curve fit to par-swap quotes reprices them to root-finder
        # precision (the fitting-time curve and the returned curve are identical
        # for local interpolants, including the pinned t=0 knot)
        swap_rates = [0.02, 0.025, 0.03]
        qs = ParSwapYield.(swap_rates, [1, 2, 3])
        c = fit(Spline.Linear(), qs, Fit.Bootstrap())
        for (r, t) in zip(swap_rates, [1, 2, 3])
            @test rate(par(c, t; frequency = 4)) ≈ r atol = 1.0e-10
        end
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
        c = fit(Yield.MonotoneConvex(), [ZCBPrice(0.95, 2.0)])
        @test discount(c, 2.0) ≈ 0.95 atol = 1.0e-8
    end
end
