@testset "FX" begin
    eurusd = FX.Pair(:EUR, :USD)

    @testset "Pair" begin
        @test eurusd === FX.Pair(:EUR, :USD)
        @test inv(eurusd) === FX.Pair(:USD, :EUR)
        @test inv(inv(eurusd)) === eurusd
    end

    usd = Yield.Constant(Continuous(0.05))
    eur = Yield.Constant(Continuous(0.03))
    S = 1.10
    m = FX.Forwards(eurusd, S, usd, eur)

    @testset "covered interest parity forwards" begin
        @test forward(m, 0.0) ≈ S
        @test forward(m, 1.0) ≈ S * exp(0.05 - 0.03)
        @test forward(m, 2.5) ≈ S * exp(0.02 * 2.5)
        @test m(1.0) == forward(m, 1.0)

        mi = inv(m)
        @test mi.pair === FX.Pair(:USD, :EUR)
        @test forward(mi, 0.0) ≈ 1 / S
        @test forward(mi, 3.0) ≈ 1 / forward(m, 3.0)
    end

    @testset "FX.Forward contract" begin
        F1 = forward(m, 1.0)
        atm = FX.Forward(eurusd, F1, 1.0)
        @test maturity(atm) == 1.0
        @test pv(m, atm) ≈ 0.0 atol = 1e-14
        # struck below the forward: worth the discounted difference in quote currency
        @test pv(m, FX.Forward(eurusd, F1 - 0.01, 1.0)) ≈ 0.01 * exp(-0.05)
        # pair mismatch errors loudly instead of pricing a crossed/inverted rate
        @test_throws ArgumentError pv(m, FX.Forward(FX.Pair(:GBP, :USD), 1.0, 1.0))
        @test_throws ArgumentError pv(m, FX.Forward(inv(eurusd), 1.0, 1.0))
    end

    @testset "quote conventions" begin
        q = FX.Outright(eurusd, 1.1225, 1.0)
        @test q.price == 0.0
        @test q.instrument.strike == 1.1225
        @test q.instrument.time == 1.0

        qp = FX.ForwardPoints(eurusd, 25.0, 0.5; spot = 1.10)
        @test qp.instrument.strike ≈ 1.10 + 25 / 10_000

        usdjpy = FX.Pair(:USD, :JPY)
        qj = FX.ForwardPoints(usdjpy, -30.0, 1.0; spot = 150.0, scale = 100)
        @test qj.instrument.strike ≈ 149.70

        # pair broadcasts as a scalar
        qs = FX.Outright.(eurusd, [1.11, 1.12], [1.0, 2.0])
        @test length(qs) == 2
        @test qs[2].instrument.time == 2.0
        qps = FX.ForwardPoints.(eurusd, [10.0, 20.0], [0.5, 1.0]; spot = 1.10)
        @test qps[1].instrument.strike ≈ 1.101
    end

    @testset "implied zero-coupon quotes" begin
        # a non-flat domestic curve proves DF_d enters and cancels correctly
        usd_boot = fit(Spline.Linear(), ZCBYield.([0.04, 0.045, 0.05, 0.052], [1.0, 2.0, 5.0, 10.0]), Fit.Bootstrap())
        m_true = FX.Forwards(eurusd, S, usd_boot, eur)
        ts = [0.5, 1.0, 2.0, 5.0, 10.0]
        quotes = [FX.Outright(eurusd, forward(m_true, t), t) for t in ts]

        implied = FX.implied_zcb_quotes(quotes, S, usd_boot)
        @test all(implied[i].price ≈ discount(eur, ts[i]) for i in eachindex(ts))
        implied_m = FX.implied_zcb_quotes(m_true, quotes)
        @test all(implied_m[i].price ≈ implied[i].price for i in eachindex(ts))

        # an off-market (nonzero-price) quote implies the same discount factor
        K = forward(m_true, 2.0) - 0.01
        off = FX.Forward(eurusd, K, 2.0)
        q_off = Quote(pv(m_true, off), off)
        @test FX.implied_zcb_quotes([q_off], S, usd_boot)[1].price ≈ discount(eur, 2.0)

        # mixed pairs are refused rather than blended into a meaningless curve
        gbp_quote = FX.Outright(FX.Pair(:GBP, :USD), 1.25, 1.0)
        @test_throws ArgumentError FX.implied_zcb_quotes([quotes[1], gbp_quote], S, usd_boot)
        @test_throws ArgumentError FX.implied_zcb_quotes(m_true, [gbp_quote])

        @testset "bootstrap fit round-trip" begin
            m_fit = fit(FX.Forwards(eurusd, S, usd_boot, Spline.Cubic()), quotes, Fit.Bootstrap())
            @test all(abs(pv(m_fit, q.instrument)) < 1e-10 for q in quotes)
            @test all(isapprox(forward(m_fit, t), forward(m_true, t)) for t in ts)
            # flat foreign truth → the interpolated curve is exact between knots too
            @test forward(m_fit, 3.3) ≈ forward(m_true, 3.3)
        end

        @testset "all-at-once spline fit (Fit.Loss) via implied quotes" begin
            m_fit = fit(FX.Forwards(eurusd, S, usd_boot, Spline.Cubic()), quotes)
            @test all(abs(pv(m_fit, q.instrument)) < 1e-6 for q in quotes)
        end
    end

    @testset "explicit cross-currency basis composition" begin
        basis = Yield.Constant(Continuous(-0.002)) # −20bp basis on the EUR leg
        mb = FX.Forwards(eurusd, S, usd, eur + basis)
        @test forward(mb, 2.0) ≈ S * exp(-(0.03 - 0.002) * 2.0) / exp(-0.05 * 2.0)
        # a negative basis cheapens base-currency discounting → higher forwards
        @test forward(mb, 2.0) > forward(m, 2.0)
    end

    @testset "parametric foreign curve via generic fit (composed optics)" begin
        ts = 0.5:0.5:5.0
        quotes = [FX.Outright(eurusd, forward(m, t), t) for t in ts]
        m_fit = fit(FX.Forwards(eurusd, S, usd, Yield.Constant()), quotes)
        @test discount(m_fit.foreign, 1.0) ≈ exp(-0.03) atol = 1e-6
        @test all(isapprox(forward(m_fit, t), forward(m, t); atol = 1e-6) for t in ts)
    end
end
