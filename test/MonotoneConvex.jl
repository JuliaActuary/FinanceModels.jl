@testset "Monotone Convex" begin

    # Interpolation of the yield curve, Gustaf Dehlbom
    # http://uu.diva-portal.org/smash/get/diva2:1477828/FULLTEXT01.pdf

    prices = [0.98, 0.955, 0.92, 0.88, 0.83]
    times = [1, 2, 3, 4, 5]
    quotes = ZCBPrice.(prices, times)
    rates = @. -log(prices) / times

    curves = [
        ("direct constructor", Yield.MonotoneConvex(rates, times)),
        ("fit", fit(Yield.MonotoneConvex(), quotes)),
    ]

    f, fᵈ = Yield.__monotone_convex_fs(rates, times)

    @testset "$name" for (name, c) in curves
        @test all(isapprox.(f, c.f; atol=1e-8))
        @test all(isapprox.(fᵈ, c.fᵈ; atol=1e-8))

        @test fᵈ[1] ≈ 0.0202 atol = 0.0001
        @test fᵈ[2] ≈ 0.0258 atol = 0.0001
        @test fᵈ[3] ≈ 0.0373 atol = 0.0001
        @test fᵈ[4] ≈ 0.0445 atol = 0.0001
        @test fᵈ[5] ≈ 0.0585 atol = 0.0001

        @test f[1] ≈ 0.0188 atol = 0.0001
        @test f[2] ≈ 0.023 atol = 0.0001
        @test f[3] ≈ 0.0316 atol = 0.0001
        @test f[4] ≈ 0.0409 atol = 0.0001
        @test f[5] ≈ 0.0515 atol = 0.0001
        @test f[6] ≈ 0.062 atol = 0.0001

        @test FinanceModels.Yield.g(0.5, 0.018793076350927487, 0.023021969250703423, 0.020202707317519466) ≈ 0.0042 * (0.5)^2 - 0.0014 atol = 0.0001
        @test FinanceModels.Yield.g(0, 0.023021969250703423, 0.03158945081076577, 0.02584123118388738) ≈ -0.0028 atol = 0.0001
        @test FinanceModels.Yield.g(0.5, 0.023021969250703423, 0.03158945081076577, 0.02584123118388738) ≈ 0.0087 * (0.5)^2 − 0.0002 * 0.5 - 0.0028 atol = 0.0001
        @test FinanceModels.Yield.g(0.5, 0.03158945081076577, 0.04089471650423902, 0.03733767043764417) ≈ -0.0063 * (0.5)^2 + 0.0156 * 0.5 - 0.0057 atol = 0.0001
        @test FinanceModels.Yield.g(0.5, 0.04089471650423902, 0.051473984626221235, 0.04445176257083387) ≈ 0.0102 * (0.5)^2 + 0.0004 * 0.5 - 0.0036 atol = 0.0001
    end

    # Reference function for zero rates (tolerance needed due to rounded coefficients)
    function r(t)
        if 0 <= t <= 1
            return 0.0014t^2 + 0.0188
        elseif 1 <= t <= 1.0233
            return -0.0028 / t + 0.023
        elseif 1.0233 <= t <= 2
            return 0.0029t^2 - 0.0088t - 0.0058 / t + 0.0319
        elseif 2 <= t <= 3
            return -0.0022t^2 + 0.0212t + 0.0324 / t - 0.0268
        elseif 3 <= t <= 4
            return 0.0031t^2 - 0.0274t - 0.1188 / t + 0.1217
        elseif 4 <= t <= 5
            return -0.0035t^2 + 0.0525t + 0.314 / t - 0.2005
        else
            error("t is out of the defined range")
        end
    end

    @testset "zero rates, knots, discounts - $name" for (name, c) in curves
        # Test zero rates against reference function
        @testset "zero rate at t=$t" for t in range(0, 5, 30)
            @test rate(zero(c, t)) ≈ r(t) atol = 0.0001
        end

        # Test that zero rates match exactly at the knot points
        @testset "knot points" begin
            for (i, t) in enumerate(times)
                @test rate(zero(c, t)) ≈ rates[i]
            end
        end

        # Test that discount factors match original prices at knot points
        @testset "discount factors" begin
            for (i, t) in enumerate(times)
                @test discount(c, t) ≈ prices[i]
            end
        end
    end


    # Tests from Google TF Quant Finance
    # https://github.com/google/tf-quant-finance/blob/master/tf_quant_finance/rates/hagan_west/monotone_convex_test.py

    @testset "TF Quant Finance - interpolated yields" begin
        # test_interpolated_yields_with_yields Example1
        # Reference times and yields (zero rates in percent)
        ref_times = [1.0, 2.0, 3.0, 4.0]
        yields = [5.0, 4.75, 4.533333, 4.775] ./ 100  # convert to decimal

        c = Yield.MonotoneConvex(yields, ref_times)

        # Interpolation times and expected yields
        interp_times = [0.25, 0.5, 1.0, 2.0, 3.0, 1.1, 2.5, 2.9, 3.6, 4.0]
        expected = [5.1171875, 5.09375, 5.0, 4.75, 4.533333, 4.9746, 4.624082, 4.535422, 4.661777, 4.775] ./ 100

        @testset "t=$t" for (i, t) in enumerate(interp_times)
            @test rate(zero(c, t)) ≈ expected[i] atol = 0.0001
        end
    end

    @testset "TF Quant Finance - yields at knot points" begin
        # test_interpolated_yields_with_yields Example2
        # Verifies exact interpolation at knot points
        ref_times = [0.1, 0.2, 0.21]
        yields = [0.1, 0.2, 0.3]

        c = Yield.MonotoneConvex(yields, ref_times)

        @testset "t=$t" for (i, t) in enumerate(ref_times)
            @test rate(zero(c, t)) ≈ yields[i] atol = 1.0e-10
        end
    end

    @testset "fit with ZCBPrice" begin
        prices = [0.98, 0.955, 0.92, 0.88, 0.83]
        times = [1, 2, 3, 4, 5]
        quotes = ZCBPrice.(prices, times)

        c = fit(Yield.MonotoneConvex(), quotes)

        # Verify discount factors match prices
        @testset "discount at t=$t" for (i, t) in enumerate(times)
            @test discount(c, t) ≈ prices[i] atol = 1.0e-6
        end
    end

    @testset "fit with ParYield" begin
        par_rates = [0.02, 0.025, 0.03, 0.035, 0.04]
        times = [1, 2, 3, 4, 5]
        quotes = ParYield.(par_rates, times)

        c = fit(Yield.MonotoneConvex(), quotes)

        # Verify par rates match
        @testset "par at t=$t" for (i, t) in enumerate(times)
            @test rate(par(c, t)) ≈ par_rates[i] atol = 0.0001
        end
    end

    @testset "forward rates" begin
        prices = [0.98, 0.955, 0.92, 0.88, 0.83]
        times = [1, 2, 3, 4, 5]
        rates = @. -log(prices) / times
        quotes = ZCBPrice.(prices, times)

        curves = [
            ("direct constructor", Yield.MonotoneConvex(rates, times)),
            ("fit", fit(Yield.MonotoneConvex(), quotes)),
        ]

        @testset "$name" for (name, c) in curves
            # Forward at t=0 should equal instantaneous forward f[1]
            @test Yield.forward(c, 0.0) ≈ c.f[1] atol = 1.0e-10

            # Forward at internal knot points equals the instantaneous forward f[i+1]
            # (at x=1 of the interval ending at that knot)
            @testset "forward at knot t=$t" for (i, t) in enumerate(times[1:(end - 1)])
                @test Yield.forward(c, t) ≈ c.f[i + 1] atol = 1.0e-10
            end

            # Forward at/beyond last time equals the last discrete forward (extrapolation)
            @test Yield.forward(c, times[end]) ≈ c.fᵈ[end] atol = 1.0e-10

            # Forward should be continuous and positive across intervals
            @testset "forward positive and continuous" begin
                ts = range(0.01, 4.99, 100)  # Stay within last interval to avoid extrapolation jump
                fwds = [Yield.forward(c, t) for t in ts]
                @test all(fwds .> 0)

                # Check continuity: adjacent forward values shouldn't jump excessively
                for i in 2:length(fwds)
                    @test abs(fwds[i] - fwds[i - 1]) < 0.01
                end
            end

            # Extrapolation: forward beyond last time equals last discrete forward
            @test Yield.forward(c, 6.0) ≈ c.fᵈ[end] atol = 1.0e-10
            @test Yield.forward(c, 10.0) ≈ c.fᵈ[end] atol = 1.0e-10

            # Mid-interval forwards should be bounded by fᵈ ± max deviation
            @testset "mid-interval forward bounds" begin
                # First interval: t ∈ (0, 1)
                f_mid = Yield.forward(c, 0.5)
                @test f_mid > 0  # Must be positive
                @test abs(f_mid - c.fᵈ[1]) < max(abs(c.f[1] - c.fᵈ[1]), abs(c.f[2] - c.fᵈ[1])) + 0.001

                # Second interval: t ∈ (1, 2)
                f_mid = Yield.forward(c, 1.5)
                @test f_mid > 0  # Must be positive
                @test abs(f_mid - c.fᵈ[2]) < max(abs(c.f[2] - c.fᵈ[2]), abs(c.f[3] - c.fᵈ[2])) + 0.001
            end
        end
    end
end
