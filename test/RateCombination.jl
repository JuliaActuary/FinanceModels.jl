@testset "Rate Combinations" begin
    @testset "base + spread" begin
        riskfree_maturities = [0.5, 1.0, 1.5, 2.0]
        riskfree = [5.0, 5.8, 6.4, 6.8] ./ 100 # spot rates

        spread_maturities = [0.5, 1.0, 1.5, 3.0] # different maturities
        spread = [1.0, 1.8, 1.4, 1.8] ./ 100 # spot spreads

        rf_curve = Yields.Zero(riskfree, riskfree_maturities)
        spread_curve = Yields.Zero(spread, spread_maturities)

        yield = rf_curve + spread_curve

        @test rate(zero(yield, 0.5)) ≈ first(riskfree) + first(spread)

        @test discount(yield, 1.0) ≈ 1 / (1 + riskfree[2] + spread[2])^1
        @test discount(yield, 1.5) ≈ 1 / (1 + riskfree[3] + spread[3])^1.5
    end
end