@testset "SmithWilson" begin

    @testset "SmithWilson" begin

        ufr = 0.03
        α = 0.1
        w = SmithWilson(ufr, α)
        u = [5.0, 7.0]
        qb = [2.3, -1.2]

        # Basic behaviour
        sw = Yields.SmithWilsonCurve(w,u, qb)
        @test sw.ufr == ufr
        @test sw.α == α
        @test sw.u == u
        @test sw.qb == qb
        @test_throws DomainError Yields.SmithWilsonCurve(w,u, [2.4, -3.4, 8.9])

        # Empty u and Qb should result in a flat yield curve
        # Use this to test methods expected from <:AbstractYieldCurve
        # Only discount and zero are explicitly implemented, so the others should follow automatically
        sw_flat = Yields.SmithWilsonCurve(w,Float64[], Float64[])
        @test discount(sw_flat, 10.0) == exp(-ufr * 10.0)
        @test accumulation(sw_flat, 10.0) ≈ exp(ufr * 10.0)
        @test rate(convert(Yields.Continuous(), zero(sw_flat, 8.0))) ≈ ufr
        @test discount.(sw_flat, [5.0, 10.0]) ≈ exp.(-ufr .* [5.0, 10.0])
        @test rate(convert(Yields.Continuous(), forward(sw_flat, 5.0, 8.0))) ≈ ufr

        # A trivial Qb vector (=0) should result in a flat yield curve
        ufr_curve = Yields.SmithWilsonCurve(w,u, [0.0, 0.0])
        @test discount(ufr_curve, 10.0) == exp(-ufr * 10.0)

        # A single payment at time 4, zero interest
        obs = [ZCBPrice(1.0,4)]
        curve_with_zero_yield = curve(w,obs)
        @test discount(curve_with_zero_yield, 4.0) == 1.0

        # In the long end it's still just UFR
        @test rate(convert(Yields.Continuous(), forward(curve_with_zero_yield, 1000.0, 2000.0))) ≈ ufr

        # Three maturities have known discount factors
        times = [1.0, 2.5, 5.6]
        prices = [0.9, 0.7, 0.5]
        cfs = [1 0 0
            0 1 0
            0 0 1]

        curve_three = Yields.SmithWilson(times, cfs, prices, ufr = ufr, α = α)
        @test transpose(cfs) * discount.(curve_three, times) ≈ prices

        # Two cash flows with payments at three times
        prices = [1.0, 0.9]
        cfs = [0.1 0.1
            1.0 0.1
            0.0 1.0]
        curve_nondiag = Yields.SmithWilson(times, cfs, prices, ufr = ufr, α = α)
        @test transpose(cfs) * discount.(curve_nondiag, times) ≈ prices

        # Round-trip zero coupon quotes
        zcq_times = [1.2, 4.5, 5.6]
        zcq_prices = [1.0, 0.9, 1.2]
        zcq = Yields.ZeroCouponQuote.(zcq_prices, zcq_times)
        sw_zcq = Yields.SmithWilson(zcq, ufr = ufr, α = α)
        @testset "ZeroCouponQuotes round-trip" for idx = 1:length(zcq_times)
            @test discount(sw_zcq, zcq_times[idx]) ≈ zcq_prices[idx]
        end

        # uneven frequencies
        swq_maturities = [1.2, 2.5, 3.6]
        swq_interests = [-0.02, 0.3, 0.04]
        frequency = [2, 1, 2]
        swq = Yields.SwapQuote.(swq_interests, swq_maturities, frequency)
        swq_times = 0.5:0.5:3.5   # Maturities are rounded down to multiples of 1/frequency, [1.0, 2.5, 3.5]
        swq_payments = [
            -0.01 0.3 0.02
            0.99 0 0.02
            0 0.3 0.02
            0 0 0.02
            0 1.3 0.02
            0 0 0.02
            0 0 1.02
        ]

        # Round-trip swap quotes
        swq_maturities = [1.2, 2.5, 3.6]
        swq_interests = [-0.02, 0.3, 0.04]
        frequency = 2
        swq = Yields.SwapQuote.(swq_interests, swq_maturities, frequency)
        swq_times = 0.5:0.5:3.5   # Maturities are rounded down to multiples of 1/frequency, [1.0, 2.5, 3.5]
        swq_payments = [-0.01 0.15 0.02
            0.99 0.15 0.02
            0.0 0.15 0.02
            0.0 0.15 0.02
            0.0 1.15 0.02
            0.0 0.0 0.02
            0.0 0.0 1.02]
        sw_swq = Yields.SmithWilson(swq, ufr = ufr, α = α)
        @testset "SwapQuotes round-trip" for swapIdx = 1:length(swq_interests)
            @test sum(discount.(sw_swq, swq_times) .* swq_payments[:, swapIdx]) ≈ 1.0
        end

        @test Yields.__ratetype(sw_swq) == Yields.Rate{Float64,typeof(Yields.DEFAULT_COMPOUNDING)}

        # Round-trip bullet bond quotes (reuse data from swap quotes)
        bbq_prices = [1.3, 0.1, 4.5]
        bbq = Yields.BulletBondQuote.(swq_interests, bbq_prices, swq_maturities, frequency)
        sw_bbq = Yields.SmithWilson(bbq, ufr = ufr, α = α)
        @testset "BulletBondQuotes round-trip" for bondIdx = 1:length(swq_interests)
            @test sum(discount.(sw_bbq, swq_times) .* swq_payments[:, bondIdx]) ≈ bbq_prices[bondIdx]
        end

        @testset "SW ForwardStarting" begin
            fwd_time = 1.0
            fwd = Yields.ForwardStarting(sw_swq, fwd_time)

            @test discount(fwd, 3.7) ≈ discount(sw_swq, fwd_time, fwd_time + 3.7)
        end

        # EIOPA risk free rate (no VA), 31 August 2021.
        # https://www.eiopa.europa.eu/sites/default/files/risk_free_interest_rate/eiopa_rfr_20210831.zip
        eiopa_output_qb = [-0.59556534586390800
            -0.07442224713453920
            -0.34193181987682400
            1.54054875814153000
            -2.15552046042343000
            0.73559290752221900
            1.89365225129089000
            -2.75927773116240000
            2.24893737130629000
            -1.51625404117395000
            0.19284859623817400
            1.13410725406271000
            0.00153268224642171
            0.00147942301778158
            -1.85022125156483000
            0.00336230229850928
            0.00324546553910162
            0.00313268874430658
            0.00302383083427276
            1.36047951448615000]
        eiopa_output_u = 1:20
        eiopa_ufr = log(1.036)
        eiopa_α = 0.133394
        sw_eiopa_expected = Yields.SmithWilson(eiopa_output_u, eiopa_output_qb; ufr = eiopa_ufr, α = eiopa_α)

        eiopa_eurswap_maturities = [1:12; 15; 20]
        eiopa_eurswap_rates = [-0.00615, -0.00575, -0.00535, -0.00485, -0.00425, -0.00375, -0.003145,
            -0.00245, -0.00185, -0.00125, -0.000711, -0.00019, 0.00111, 0.00215]   # Reverse engineered from output curve. This is the full precision of market quotes.
        eiopa_eurswap_quotes = Yields.SwapQuote.(eiopa_eurswap_rates, eiopa_eurswap_maturities, 1)
        sw_eiopa_actual = Yields.SmithWilson(eiopa_eurswap_quotes, ufr = eiopa_ufr, α = eiopa_α)

        @testset "Match EIOPA calculation" begin
            @test sw_eiopa_expected.u ≈ sw_eiopa_actual.u
            @test sw_eiopa_expected.qb ≈ sw_eiopa_actual.qb
        end
    end
end