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

    @testset "foreign-curve optic composition (unbounded forms)" begin
        # the docstring-documented unbounded 1-tuple form and a bare optic both compose
        # through the `foreign` field (the bounded `optic => interval` form is covered
        # via the parametric `fit` test above)
        m0 = FX.Forwards(eurusd, 1.0, usd, eur)
        o_tuple = FinanceModels.__fx_foreign_optic((@optic(_.rate),))
        @test only(o_tuple)(m0) === m0.foreign.rate
        o_bare = FinanceModels.__fx_foreign_optic(@optic(_.rate))
        @test o_bare(m0) === m0.foreign.rate
    end

    @testset "textbook: covered interest parity (Hull §5.10)" begin
        # Hull, "Options, Futures, and Other Derivatives" (9th ed.), §5.10: with the
        # spot AUD/USD rate at 0.6200 and 2-year continuously-compounded rates of 5%
        # (AUD) and 7% (USD), the 2-year forward is 0.62·e^{(0.07−0.05)·2} = 0.6453.
        audusd = FX.Pair(:AUD, :USD)
        hull = FX.Forwards(audusd, 0.6200, Yield.Constant(Continuous(0.07)), Yield.Constant(Continuous(0.05)))
        @test forward(hull, 2.0) ≈ 0.6453 atol = 5e-5
        # covered interest arbitrage: a forward struck below/above the CIP rate has
        # positive/negative value to the buyer of AUD
        @test pv(hull, FX.Forward(audusd, 0.6300, 2.0)) > 0
        @test pv(hull, FX.Forward(audusd, 0.6600, 2.0)) < 0
    end

    @testset "interest rate parity with periodic compounding" begin
        # the discrete-compounding IRP form: F(n) = S·(1+r_d)ⁿ/(1+r_f)ⁿ with
        # annually-compounded rates
        m1 = FX.Forwards(eurusd, 1.10, Yield.Constant(Periodic(0.02, 1)), Yield.Constant(Periodic(0.01, 1)))
        @test forward(m1, 1.0) ≈ 1.10 * 1.02 / 1.01
        @test forward(m1, 2.0) ≈ 1.10 * 1.02^2 / 1.01^2
    end

    @testset "FX.Converted multi-currency projection" begin
        usd_boot = fit(Spline.Linear(), ZCBYield.([0.04, 0.045, 0.05], [1.0, 3.0, 5.0]), Fit.Bootstrap())
        fx = FX.Forwards(eurusd, 1.08, usd_boot, eur)
        store = Dict("EURUSD" => fx, "ESTR" => eur)
        bond = Bond.Fixed(0.04, Periodic(2), 5.0)

        # each converted cashflow is the original amount times the CIP forward at its time
        cfs = collect(Projection(FX.Converted(bond, "EURUSD"), store, CashflowProjection()))
        raw = collect(bond)
        @test length(cfs) == length(raw)
        @test all(cfs[i].time == raw[i].time for i in eachindex(raw))
        @test all(cfs[i].amount ≈ raw[i].amount * forward(fx, raw[i].time) for i in eachindex(raw))

        # fundamental identity: converting at CIP forwards then discounting on the
        # domestic curve equals discounting on the foreign curve then converting at spot
        @test pv(usd_boot, Projection(FX.Converted(bond, "EURUSD"), store, CashflowProjection())) ≈ 1.08 * pv(eur, bond)

        # a converted floating leg resolves its own reference curve from the same store
        frn = Bond.Floating(0.001, Periodic(4), 2.0, "ESTR")
        lhs = pv(usd_boot, Projection(FX.Converted(frn, "EURUSD"), store, CashflowProjection()))
        rhs = 1.08 * pv(eur, Projection(frn, store, CashflowProjection()))
        @test lhs ≈ rhs

        # transducer-modified (Eduction) contracts convert too
        scaled = pv(usd_boot, Projection(FX.Converted(bond |> Map(cf -> cf * 100.0), "EURUSD"), store, CashflowProjection()))
        @test scaled ≈ 100.0 * 1.08 * pv(eur, bond)

        @test maturity(FX.Converted(bond, "EURUSD")) == 5.0
        # a missing model key errors loudly
        @test_throws KeyError collect(Projection(FX.Converted(bond, "GBPUSD"), store, CashflowProjection()))
    end

    @testset "textbook: fixed-for-fixed currency swap (Hull §7.9)" begin
        # Hull, "Options, Futures, and Other Derivatives" (9th ed.), §7.9: flat term
        # structures of 9% (USD) and 4% (JPY), continuously compounded; spot ¥110 = $1.
        # An institution receives 5% on ¥1,200M and pays 8% on $10M annually for 3 more
        # years, with principals exchanged at maturity. Hull values the swap at $1.543M
        # two ways: as two bonds converted at spot, and as a portfolio of FX forwards.
        jpyusd = FX.Pair(:JPY, :USD)
        usd9 = Yield.Constant(Continuous(0.09))
        jpy4 = Yield.Constant(Continuous(0.04))
        spot = 1 / 110 # USD per JPY
        fx = FX.Forwards(jpyusd, spot, usd9, jpy4)

        # the forward exchange rates of Hull's Table 7.9
        @test forward(fx, 1.0) ≈ 0.009557 atol = 5e-7
        @test forward(fx, 2.0) ≈ 0.010047 atol = 5e-7
        @test forward(fx, 3.0) ≈ 0.010562 atol = 5e-7

        yen_leg = FX.Converted(Bond.Fixed(0.05, Periodic(1), 3) |> Map(cf -> cf * 1200.0), "JPYUSD")
        usd_leg = Bond.Fixed(0.08, Periodic(1), 3) |> Map(cf -> cf * -10.0)
        swap = Composite(yen_leg, usd_leg)
        p = Projection(swap, Dict("JPYUSD" => fx), CashflowProjection())

        # (a) the portfolio-of-forward-contracts valuation — our projection route
        @test pv(usd9, p) ≈ 1.5430 atol = 1e-4
        # (b) the two-bond valuation: B_F/110 − B_D = 1,230.55/110 − 9.6439
        B_D = 10.0 * pv(usd9, Bond.Fixed(0.08, Periodic(1), 3))
        B_F = 1200.0 * pv(jpy4, Bond.Fixed(0.05, Periodic(1), 3))
        @test B_D ≈ 9.6439 atol = 1e-4
        @test B_F ≈ 1230.55 atol = 1e-2
        # the two routes agree to numerical precision
        @test pv(usd9, p) ≈ B_F * spot - B_D
    end
end
