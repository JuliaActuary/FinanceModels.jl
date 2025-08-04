@testset "SmithWilson" begin

    ufr = 0.03
    α = 0.1
    u = [5.0, 7.0]
    qb = [2.3, -1.2]

    # Basic behaviour
    sw = Yield.SmithWilson(u, qb; ufr = ufr, α = α)
    @test sw.ufr == ufr
    @test sw.α == α
    @test sw.u == u
    @test sw.qb == qb
    @test_throws DomainError Yield.SmithWilson(u, [2.4, -3.4, 8.9], ufr = ufr, α = α)

    # Empty u and Qb should result in a flat yield curve
    # Use this to test methods expected from <:AbstractYieldCurve
    # Only discount and zero are explicitly implemented, so the others should follow automatically
    sw_flat = Yield.SmithWilson(Float64[], Float64[]; ufr = ufr, α = α)
    @test discount(sw_flat, 10.0) == exp(-ufr * 10.0)
    @test accumulation(sw_flat, 10.0) ≈ exp(ufr * 10.0)
    @test zero(sw_flat, 8.0) ≈ Continuous(ufr)
    @test discount.(sw_flat, [5.0, 10.0]) ≈ exp.(-ufr .* [5.0, 10.0])
    @test forward(sw_flat, 5.0, 8.0) ≈ Continuous(ufr)

    # A trivial Qb vector (=0) should result in a flat yield curve
    ufr_curve = Yield.SmithWilson(u, [0.0, 0.0]; ufr = ufr, α = α)
    @test discount(ufr_curve, 10.0) == exp(-ufr * 10.0)

    # A single payment at time 4, zero interest
    curve_with_zero_yield = Yield.SmithWilson([4.0], reshape([1.0], 1, 1), [1.0]; ufr = ufr, α = α)
    @test discount(curve_with_zero_yield, 4.0) == 1.0

    # In the long end it's still just UFR
    @test forward(curve_with_zero_yield, 1000.0, 2000.0) ≈ Continuous(ufr)

    # Three maturities have known discount factors
    times = [1.0, 2.5, 5.6]
    prices = [0.9, 0.7, 0.5]
    qs = ZCBPrice.(prices, times)


    curve_three = fit(Yield.SmithWilson(ufr = ufr, α = α), qs)
    @test [pv(curve_three, q.instrument) for q in qs] ≈ prices

    # Two cash flows with payments at three times
    prices = [1.0, 0.9]
    times = [1.0, 2.5, 5.6]
    cfs = [
        0.1 0.1
        1.1 0.1
        0.0 1.1
    ]
    qs = [Quote(q[1], Cashflow.(q[2], times)) for q in zip(prices, eachcol(cfs))]
    curve_nondiag = Yield.SmithWilson(times, cfs, prices; ufr = ufr, α = α)
    @test transpose(cfs) * discount.(curve_nondiag, times) ≈ prices
    curve_nondiag = fit(Yield.SmithWilson(ufr = ufr, α = α), qs)
    @test transpose(cfs) * discount.(curve_nondiag, times) ≈ prices

    # Round-trip zero coupon quotes
    zcq_times = [1.2, 4.5, 5.6]
    zcq_prices = [1.0, 0.9, 1.2]
    qs = ZCBPrice.(zcq_prices, zcq_times)
    sw_zcq = fit(Yield.SmithWilson(ufr = ufr, α = α), qs)
    @testset "ZeroCouponQuotes round-trip" for idx in 1:length(zcq_times)
        @test discount(sw_zcq, zcq_times[idx]) ≈ zcq_prices[idx]
    end

    # Round-trip swap quotes
    maturities = [1.2, 2.5, 3.6]
    coupon = [-0.02, 0.3, 0.04]
    frequency = Periodic.(2)
    qs = Quote.(
        ones(length(maturities)),
        Bond.Fixed.(coupon, frequency, maturities)
    )

    sw_swq = fit(Yield.SmithWilson(ufr = ufr, α = α), qs)
    swq_payments, swq_times = FinanceModels.cashflows_timepoints(qs)
    @testset "SwapQuotes round-trip" for swapIdx in 1:length(coupon)
        @test sum(discount.(sw_swq, swq_times) .* swq_payments[:, swapIdx]) ≈ 1.0
    end
    @testset "SW ForwardStarting" begin
        fwd_time = 1.0
        fwd = Yield.ForwardStarting(sw_swq, fwd_time)

        @test discount(fwd, 3.7) ≈ discount(sw_swq, fwd_time, fwd_time + 3.7)
    end

    # Round-trip bullet bond quotes (reuse data from swap quotes)
    bbq_prices = [1.3, 0.1, 4.5]
    qs = Quote.(
        bbq_prices,
        Bond.Fixed.(coupon, frequency, maturities)
    )
    sw_bbq = fit(Yield.SmithWilson(ufr = ufr, α = α), qs)
    @testset "BulletBondQuotes round-trip" for bondIdx in 1:length(bbq_prices)
        @test sum(discount.(sw_bbq, swq_times) .* swq_payments[:, bondIdx]) ≈ bbq_prices[bondIdx]
    end


    # EIOPA risk free rate (no VA), 31 August 2021.
    # https://www.eiopa.europa.eu/sites/default/files/risk_free_interest_rate/eiopa_rfr_20210831.zip
    eiopa_output_qb = [
        -0.595565345863908
        -0.0744222471345392
        -0.341931819876824
        1.54054875814153
        -2.15552046042343
        0.735592907522219
        1.89365225129089
        -2.7592777311624
        2.24893737130629
        -1.51625404117395
        0.192848596238174
        1.13410725406271
        0.00153268224642171
        0.00147942301778158
        -1.85022125156483
        0.00336230229850928
        0.00324546553910162
        0.00313268874430658
        0.00302383083427276
        1.36047951448615
    ]
    eiopa_output_u = 1:20
    eiopa_ufr = log(1.036)
    eiopa_α = 0.133394
    sw_eiopa_expected = Yield.SmithWilson(eiopa_output_u, eiopa_output_qb; ufr = eiopa_ufr, α = eiopa_α)

    eiopa_eurswap_maturities = [1:12; 15; 20]
    # Reverse engineered from output curve. This is the full precision of market quotes.
    eiopa_eurswap_rates = [
        -0.00615, -0.00575, -0.00535, -0.00485, -0.00425, -0.00375, -0.003145,
        -0.00245, -0.00185, -0.00125, -0.000711, -0.00019, 0.00111, 0.00215,
    ]
    eiopa_eurswap_quotes = Quote.(1.0, Bond.Fixed.(eiopa_eurswap_rates, Periodic(1), eiopa_eurswap_maturities))
    sw_eiopa_actual = fit(Yield.SmithWilson(ufr = eiopa_ufr, α = eiopa_α), eiopa_eurswap_quotes)

    @testset "Match EIOPA calculation" begin
        @test sw_eiopa_expected.u ≈ sw_eiopa_actual.u
        @test sw_eiopa_expected.qb ≈ sw_eiopa_actual.qb
    end
end
