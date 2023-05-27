@testset "Rate Combinations" begin
    riskfree_maturities = [0.5, 1.0, 1.5, 2.0]
    riskfree = [5.0, 5.8, 6.4, 6.8] ./ 100 # spot rates
    rf_curve = FinanceModels.Zero(riskfree, riskfree_maturities)

    @testset "base + spread" begin

        spread_maturities = [0.5, 1.0, 1.5, 3.0] # different maturities
        spread = [1.0, 1.8, 1.4, 1.8] ./ 100 # spot spreads

        spread_curve = FinanceModels.Zero(spread, spread_maturities)

        yield = rf_curve + spread_curve

        @test rate(zero(yield, 0.5)) ≈ first(riskfree) + first(spread)

        @test discount(yield, 1.0) ≈ 1 / (1 + riskfree[2] + spread[2])^1
        @test discount(yield, 1.5) ≈ 1 / (1 + riskfree[3] + spread[3])^1.5
    end

    @testset "multiplicaiton and division" begin
        @testset "multiplication" begin
            factor = .79
            c = rf_curve * factor
            (discount(c,10)-1 * factor) ≈ discount(rf_curve,10)
            (accumulation(c,10)-1 * factor) ≈ accumulation(rf_curve,10)
            forward(c,5,10) * factor ≈ forward(rf_curve,5,10)
            FinanceModels.par(c,10) * factor ≈ FinanceModels.par(rf_curve,10)

            c = factor * rf_curve
            (discount(c,10)-1 * factor) ≈ discount(rf_curve,10)
            (accumulation(c,10)-1 * factor) ≈ accumulation(rf_curve,10)
            forward(c,5,10) * factor ≈ forward(rf_curve,5,10)
            FinanceModels.par(c,10) * factor ≈ FinanceModels.par(rf_curve,10)

            @test discount(FinanceModels.Constant(0.1) * FinanceModels.Constant(0.1),10) ≈ discount(FinanceModels.Constant(0.01),10)
        end

        @testset "division" begin
            factor = .79
            c = rf_curve / (factor^-1)
            (discount(c,10)-1 * factor) ≈ discount(rf_curve,10)
            (accumulation(c,10)-1 * factor) ≈ accumulation(rf_curve,10)
            forward(c,5,10) * factor ≈ forward(rf_curve,5,10)
            FinanceModels.par(c,10) * factor ≈ FinanceModels.par(rf_curve,10)
            @test discount(FinanceModels.Constant(0.1) / FinanceModels.Constant(0.5),10) ≈ discount(FinanceModels.Constant(0.2),10)
            @test discount(0.1 / FinanceModels.Constant(0.5),10) ≈ discount(FinanceModels.Constant(0.2),10)
        end
    end
end