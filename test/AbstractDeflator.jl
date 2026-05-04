@testset "AbstractDeflator integration" begin
    @testset "AbstractYieldModel <: AbstractDeflator" begin
        @test FinanceModels.Yield.AbstractYieldModel <: FinanceCore.AbstractDeflator
        c = FinanceModels.Yield.Constant(0.04)
        @test c isa FinanceCore.AbstractDeflator
    end

    @testset "factor agrees with discount on yield curves" begin
        # Use a continuous-force constant yield so we can compare to exp(-r*t)
        c = FinanceModels.Yield.Constant(FinanceCore.Continuous(0.04))
        @test FinanceCore.factor(c, 5) ≈ FinanceCore.discount(c, 5)
        @test FinanceCore.factor(c, 1, 6) ≈ FinanceCore.discount(c, 1, 6)
        @test FinanceCore.factor(c, 1, 6) ≈ exp(-0.04 * 5)

        # Spline yield curve, a non-trivial term structure
        zc = FinanceModels.fit(
            FinanceModels.Spline.Linear(),
            FinanceModels.ZCBYield.([0.03, 0.035, 0.04, 0.045], [1, 2, 5, 10]),
            FinanceModels.Fit.Bootstrap()
        )
        @test FinanceCore.factor(zc, 5) ≈ FinanceCore.discount(zc, 5)
        @test FinanceCore.factor(zc, 1, 6) ≈ FinanceCore.discount(zc, 1, 6)
    end

    @testset "compose yield × decrement" begin
        # The headline workflow: yield × default-survival composition.
        # Use Continuous(0.03) so the composite force is exactly 0.042 = 0.03 + 0.012.
        yield = FinanceModels.Yield.Constant(FinanceCore.Continuous(0.03))
        λ_d   = FinanceCore.Continuous(0.012)         # 1.2% default hazard
        deflator = FinanceCore.compose(yield, λ_d)

        @test deflator isa FinanceCore.CompositeDeflator
        @test FinanceCore.factor(deflator, 0, 5) ≈
            FinanceCore.factor(yield, 0, 5) * FinanceCore.factor(λ_d, 0, 5)
        @test FinanceCore.factor(deflator, 0, 5) ≈ exp(-0.042 * 5)

        # pv integration with a Cashflow vector
        cashflows = [FinanceCore.Cashflow(100.0, t) for t in 1.0:5.0]
        expected = sum(100.0 * exp(-0.042 * t) for t in 1.0:5.0)
        @test FinanceCore.pv(deflator, cashflows) ≈ expected
    end

    @testset "term-structure yield × constant decrement" begin
        # Real term-structure yield curve composed with mortality force
        zc = FinanceModels.fit(
            FinanceModels.Spline.Cubic(),
            FinanceModels.ZCBYield.([0.03, 0.035, 0.04, 0.045, 0.05], [1, 2, 3, 5, 10]),
            FinanceModels.Fit.Bootstrap()
        )
        μ = FinanceCore.Continuous(0.012)
        deflator = FinanceCore.compose(zc, μ)
        # Each cashflow's factor combines the curve's spot discount with
        # the mortality survival from time 0
        for t in (1.0, 5.0, 10.0)
            @test FinanceCore.factor(deflator, t) ≈
                FinanceCore.factor(zc, t) * FinanceCore.factor(μ, t)
        end
    end

    @testset "Quote dispatch works for yield models (AbstractDeflator)" begin
        # Quote(model, contract) should still work for yield models even though
        # they no longer subtype AbstractModel — we add an AbstractDeflator method.
        yc = FinanceModels.Yield.Constant(FinanceCore.Continuous(0.04))
        # Use the Bond.Fixed contract from a ZCBYield quote as the contract argument
        zcb_quote = FinanceModels.ZCBYield(0.04, 5)
        bond = zcb_quote.instrument
        q = FinanceCore.Quote(yc, bond)
        @test q isa FinanceCore.Quote
        @test q.price ≈ FinanceCore.pv(yc, bond)
    end
end
