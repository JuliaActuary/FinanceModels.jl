@testset "Rate Combinations" begin
    riskfree_maturities = [0.5, 1.0, 1.5, 2.0]
    riskfree = [5.0, 5.8, 6.4, 6.8] ./ 100 # spot rates
    rf_curve = fit(Spline.Linear(), ZCBYield.(riskfree, riskfree_maturities), Fit.Bootstrap())

    @testset "base + spread" begin

        spread_maturities = [0.5, 1.0, 1.5, 3.0] # different maturities
        spread = [1.0, 1.8, 1.4, 1.8] ./ 100 # spot spreads


        spread_curve = fit(Spline.Linear(), ZCBYield.(spread, spread_maturities), Fit.Bootstrap())

        yield = rf_curve + spread_curve

        @test zero(yield, 0.5) ≈ Periodic(first(riskfree) + first(spread), 1)

        @test discount(yield, 1.0) ≈ 1 / (1 + riskfree[2] + spread[2])^1
        @test discount(yield, 1.5) ≈ 1 / (1 + riskfree[3] + spread[3])^1.5

        rates = ZCBYield([0.01, 0.02, 0.03])
        spreads = ZCBYield([0.02, 0.03, 0.04])
        yields = ZCBYield([0.03, 0.05, 0.07])

        r = fit(Spline.Linear(), rates, Fit.Bootstrap())
        s = fit(Spline.Linear(), spreads, Fit.Bootstrap())
        y = fit(Spline.Linear(), yields, Fit.Bootstrap())

        @test discount(r + s, 1) ≈ discount(y, 1)
        @test discount(r + s, 2) ≈ discount(y, 2)
        @test discount(r + s, 3) ≈ discount(y, 3)

        rates = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.4, 1.74, 2.31, 2.41] ./ 100
        spreads = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.4, 1.74, 2.31, 2.41] ./ 100
        mats = [1 / 12, 2 / 12, 3 / 12, 6 / 12, 1, 2, 3, 5, 7, 10, 20, 30]


        ### Zero coupon rates/spreads

        q_rf_z = ZCBYield.(rates, mats)
        q_s_z = ZCBYield.(spreads, mats)
        q_y_z = ZCBYield.(rates + spreads, mats)

        c_rf_z = fit(Spline.Linear(), q_rf_z, Fit.Bootstrap())
        c_s_z = fit(Spline.Linear(), q_s_z, Fit.Bootstrap())
        c_y_z = fit(Spline.Linear(), q_y_z, Fit.Bootstrap())

        # adding curves when the spreads were zero spreads works
        @test discount(c_rf_z + c_s_z, 20) ≈ discount(c_y_z, 20)


        ### Par coupon rates/spreads


        q_rf = CMTYield.(rates, mats)
        q_s = CMTYield.(spreads, mats)
        q_y = CMTYield.(rates + spreads, mats)

        c_rf = fit(Spline.Linear(), q_rf, Fit.Bootstrap())
        c_s = fit(Spline.Linear(), q_s, Fit.Bootstrap())
        c_y = fit(Spline.Linear(), q_y, Fit.Bootstrap())

        # adding curves when the spreads were par spreads does not work
        @test !(discount(c_rf + c_s, 20) ≈ discount(c_y, 20))


    end

    @testset "multiplicaiton and division" begin
        @testset "multiplication" begin
            factor = 0.79
            c = rf_curve * factor
            target_curve = fit(Spline.Linear(), ZCBYield.(riskfree .* factor, riskfree_maturities), Fit.Bootstrap())
            @test discount(c, 2) ≈ discount(target_curve, 2)
            @test accumulation(c, 2) ≈ accumulation(target_curve, 2)
            @test forward(c, 1, 2) ≈ forward(target_curve, 1, 2)
            @test par(c, 2) ≈ par(target_curve, 2)

            c = factor * rf_curve
            @test discount(c, 2) ≈ discount(target_curve, 2)
            @test accumulation(c, 2) ≈ accumulation(target_curve, 2)
            @test forward(c, 1, 2) ≈ forward(target_curve, 1, 2)
            @test par(c, 2) ≈ par(target_curve, 2)

            @test discount(Yield.Constant(0.1) * Yield.Constant(0.1), 10) ≈ discount(Yield.Constant(0.01), 10)
        end

        @testset "division" begin
            factor = 0.79
            c = rf_curve / factor
            target_curve = fit(Spline.Linear(), ZCBYield.(riskfree ./ factor, riskfree_maturities), Fit.Bootstrap())
            @test discount(c, 2) ≈ discount(target_curve, 2)
            @test accumulation(c, 2) ≈ accumulation(target_curve, 2)
            @test forward(c, 1, 2) ≈ forward(target_curve, 1, 2)
            @test par(c, 2) ≈ par(target_curve, 2)

            @test discount(Yield.Constant(0.1) / Yield.Constant(0.5), 10) ≈ discount(Yield.Constant(0.2), 10)
            @test discount(0.1 / Yield.Constant(0.5), 10) ≈ discount(Yield.Constant(0.2), 10)
        end
    end
end
