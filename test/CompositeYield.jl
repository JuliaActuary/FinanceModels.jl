@testset "Rate Combinations" begin
    riskfree_maturities = [0.5, 1.0, 1.5, 2.0]
    riskfree = [5.0, 5.8, 6.4, 6.8] ./ 100 # spot rates
    rf_curve = fit(Spline.Linear(), ZCBYield.(riskfree, riskfree_maturities), Fit.Bootstrap())

    @testset "base + spread" begin

        spread_maturities = [0.5, 1.0, 1.5, 3.0] # different maturities
        spread = [1.0, 1.8, 1.4, 1.8] ./ 100 # spot spreads


        spread_curve = fit(Spline.Linear(), ZCBYield.(spread, spread_maturities), Fit.Bootstrap())

        yield = rf_curve + spread_curve

        # + operates on continuous zero rates, equivalent to multiplicative discount factors
        @test discount(yield, 0.5) ≈ discount(rf_curve, 0.5) * discount(spread_curve, 0.5)
        @test discount(yield, 1.0) ≈ discount(rf_curve, 1.0) * discount(spread_curve, 1.0)
        @test discount(yield, 1.5) ≈ discount(rf_curve, 1.5) * discount(spread_curve, 1.5)

        rates = ZCBYield([0.01, 0.02, 0.03])
        spreads = ZCBYield([0.02, 0.03, 0.04])

        r = fit(Spline.Linear(), rates, Fit.Bootstrap())
        s = fit(Spline.Linear(), spreads, Fit.Bootstrap())

        @test discount(r + s, 1) ≈ discount(r, 1) * discount(s, 1)
        @test discount(r + s, 2) ≈ discount(r, 2) * discount(s, 2)
        @test discount(r + s, 3) ≈ discount(r, 3) * discount(s, 3)

        rates = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.4, 1.74, 2.31, 2.41] ./ 100
        spreads = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.4, 1.74, 2.31, 2.41] ./ 100
        mats = [1 / 12, 2 / 12, 3 / 12, 6 / 12, 1, 2, 3, 5, 7, 10, 20, 30]


        ### Zero coupon rates/spreads

        q_rf_z = ZCBYield.(rates, mats)
        q_s_z = ZCBYield.(spreads, mats)

        c_rf_z = fit(Spline.Linear(), q_rf_z, Fit.Bootstrap())
        c_s_z = fit(Spline.Linear(), q_s_z, Fit.Bootstrap())

        # adding curves produces multiplicative discount factors
        @test discount(c_rf_z + c_s_z, 20) ≈ discount(c_rf_z, 20) * discount(c_s_z, 20)


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
            # In continuous zero-rate space: z_total(t) = z_rf(t) * z_factor
            # where z_factor = log(1 + factor) for Periodic(1) factor
            z_factor = log(1 + factor)
            for t in riskfree_maturities
                z_rf = -log(discount(rf_curve, t)) / t
                @test discount(c, t) ≈ exp(-(z_rf * z_factor) * t)
            end
            @test accumulation(c, 2) ≈ 1 / discount(c, 2)

            c = factor * rf_curve
            for t in riskfree_maturities
                z_rf = -log(discount(rf_curve, t)) / t
                @test discount(c, t) ≈ exp(-(z_rf * z_factor) * t)
            end

            # Constant * Constant also operates in continuous zero-rate space
            z1 = log(1.1)  # continuous zero rate of Constant(0.1)
            @test discount(Yield.Constant(0.1) * Yield.Constant(0.1), 10) ≈ exp(-(z1 * z1) * 10)
        end

        @testset "division" begin
            factor = 0.79
            c = rf_curve / factor
            z_factor = log(1 + factor)
            for t in riskfree_maturities
                z_rf = -log(discount(rf_curve, t)) / t
                @test discount(c, t) ≈ exp(-(z_rf / z_factor) * t)
            end

            # Constant / Constant also operates in continuous zero-rate space
            z_a = log(1.1)  # continuous zero rate of Constant(0.1)
            z_b = log(1.5)  # continuous zero rate of Constant(0.5)
            @test discount(Yield.Constant(0.1) / Yield.Constant(0.5), 10) ≈ exp(-(z_a / z_b) * 10)
            @test discount(0.1 / Yield.Constant(0.5), 10) ≈ exp(-(z_a / z_b) * 10)
        end
    end
end
