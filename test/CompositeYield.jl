@testset "Rate Combinations" begin
    riskfree_maturities = [0.5, 1.0, 1.5, 2.0]
    riskfree = [5.0, 5.8, 6.4, 6.8] ./ 100 # spot rates
    rf_curve = fit(Spline.Cubic(), ZCBYield.(riskfree, riskfree_maturities), Fit.Bootstrap())

    @testset "base + spread" begin

        spread_maturities = [0.5, 1.0, 1.5, 3.0] # different maturities
        spread = [1.0, 1.8, 1.4, 1.8] ./ 100 # spot spreads


        spread_curve = rf_curve = fit(Spline.Cubic(), ZCBYield.(spread, spread_maturities), Fit.Bootstrap())

        yield = rf_curve + spread_curve

        @test zero(yield, 0.5) ≈ Periodic(first(riskfree) + first(spread), 1)

        @test discount(yield, 1.0) ≈ 1 / (1 + riskfree[2] + spread[2])^1
        @test discount(yield, 1.5) ≈ 1 / (1 + riskfree[3] + spread[3])^1.5
    end

    @testset "multiplicaiton and division" begin
        @testset "multiplication" begin
            factor = 0.79
            c = rf_curve * factor
            (discount(c, 10) - 1 * factor) ≈ discount(rf_curve, 10)
            (accumulation(c, 10) - 1 * factor) ≈ accumulation(rf_curve, 10)
            forward(c, 5, 10) * factor ≈ forward(rf_curve, 5, 10)
            par(c, 10) * factor ≈ par(rf_curve, 10)

            c = factor * rf_curve
            (discount(c, 10) - 1 * factor) ≈ discount(rf_curve, 10)
            (accumulation(c, 10) - 1 * factor) ≈ accumulation(rf_curve, 10)
            forward(c, 5, 10) * factor ≈ forward(rf_curve, 5, 10)
            par(c, 10) * factor ≈ par(rf_curve, 10)

            @test discount(Yield.Constant(0.1) * Yield.Constant(0.1), 10) ≈ discount(Yield.Constant(0.01), 10)
        end

        @testset "division" begin
            factor = 0.79
            c = rf_curve / (factor^-1)
            (discount(c, 10) - 1 * factor) ≈ discount(rf_curve, 10)
            (accumulation(c, 10) - 1 * factor) ≈ accumulation(rf_curve, 10)
            forward(c, 5, 10) * factor ≈ forward(rf_curve, 5, 10)
            par(c, 10) * factor ≈ par(rf_curve, 10)
            @test discount(Yield.Constant(0.1) / Yield.Constant(0.5), 10) ≈ discount(Yield.Constant(0.2), 10)
            @test discount(0.1 / Yield.Constant(0.5), 10) ≈ discount(Yield.Constant(0.2), 10)
        end
    end
end