@testset "FX" begin
    eurusd = FX.Pair(:EUR, :USD)

    @testset "Pair" begin
        @test eurusd == FX.Pair(:EUR, :USD)
        @test eurusd === FX.Pair(:EUR, :USD) # immutable value: egal too
        @test inv(eurusd) == FX.Pair(:USD, :EUR)
        @test inv(inv(eurusd)) == eurusd
        @test repr(eurusd) == "FX.Pair(:EUR, :USD)"

        # currencies are not restricted to Symbols: any values work, e.g. ISO 4217
        # numeric codes, and the whole pipeline carries them along
        iso = FX.Pair(978, 840) # EUR/USD by ISO numeric code
        @test inv(iso) == FX.Pair(840, 978)
        @test repr(iso) == "FX.Pair(978, 840)"
        m_iso = FX.Forwards(iso, 1.10, Yield.Constant(Continuous(0.05)), Yield.Constant(Continuous(0.03)))
        @test forward(m_iso, 1.0) ≈ 1.10 * exp(0.02)
        @test pv(m_iso, FX.Forward(iso, forward(m_iso, 1.0), 1.0)) ≈ 0.0 atol = 1e-14
        # a Symbol pair and an integer-code pair are distinct pairs, per the
        # direction-guard design
        @test_throws ArgumentError pv(m_iso, FX.Forward(FX.Pair(:EUR, :USD), 1.0, 1.0))

        # plain strings work too (currencies are values, not type parameters)
        s = FX.Pair("EUR", "USD")
        @test inv(s) == FX.Pair("USD", "EUR")
        @test repr(s) == "FX.Pair(\"EUR\", \"USD\")"
        # equality and hashing are content-based across string storage types (the
        # fixed widths CSV readers assign per column, or substrings as here), so
        # contracts and models match whenever the currencies read the same
        sub = FX.Pair(SubString("EURO", 1, 3), SubString("USDX", 1, 3))
        @test sub == s
        @test hash(sub) == hash(s)
        m_str = FX.Forwards(s, 1.10, Yield.Constant(Continuous(0.05)), Yield.Constant(Continuous(0.03)))
        @test pv(m_str, FX.Forward(sub, forward(m_str, 1.0), 1.0)) ≈ 0.0 atol = 1e-14
        # ...but a Symbol pair and a string pair are distinct conventions
        @test s != eurusd
        @test_throws ArgumentError pv(m_str, FX.Forward(eurusd, 1.0, 1.0))
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
        @test mi.pair == FX.Pair(:USD, :EUR)
        @test forward(mi, 0.0) ≈ 1 / S
        @test forward(mi, 3.0) ≈ 1 / forward(m, 3.0)
        # pricing under the inverted model: the role swap is complete, so an
        # off-market USDEUR forward discounts on the (now-domestic) EUR curve
        K = forward(mi, 1.0) - 0.005
        @test pv(mi, FX.Forward(FX.Pair(:USD, :EUR), K, 1.0)) ≈ 0.005 * exp(-0.03)
    end

    @testset "FX.Forward contract" begin
        F1 = forward(m, 1.0)
        atm = FX.Forward(eurusd, F1, 1.0)
        @test maturity(atm) == 1.0
        @test pv(m, atm) ≈ 0.0 atol = 1e-14
        # struck below the forward: worth the discounted difference in quote currency
        off = FX.Forward(eurusd, F1 - 0.01, 1.0)
        @test pv(m, off) ≈ 0.01 * exp(-0.05)
        # valuation at cur_time > 0: same payoff, discounted over the remaining period
        @test pv(m, off, 0.5) ≈ 0.01 * exp(-0.05 * 0.5)
        # inclusive at cur_time == settlement: the undiscounted payoff (matches the
        # `cf.time >= cur_time` projection filter convention)
        @test pv(m, off, 1.0) ≈ 0.01
        # a forward that settled before the valuation time is worth zero
        @test pv(m, off, 2.0) == 0.0
        # pair mismatch errors loudly instead of pricing a crossed/inverted rate
        @test_throws ArgumentError pv(m, FX.Forward(FX.Pair(:GBP, :USD), 1.0, 1.0))
        @test_throws ArgumentError pv(m, FX.Forward(inv(eurusd), 1.0, 1.0))
        # ...including through the valuation-time form
        @test_throws ArgumentError pv(m, FX.Forward(FX.Pair(:GBP, :USD), 1.0, 1.0), 0.5)
    end

    @testset "quote conventions" begin
        q = FX.Outright(eurusd, 1.1225, 1.0)
        @test q.price == 0.0
        @test q.instrument.strike == 1.1225
        @test q.instrument.time == 1.0

        # pair broadcasts as a scalar
        qs = FX.Outright.(eurusd, [1.11, 1.12], [1.0, 2.0])
        @test length(qs) == 2
        @test qs[2].instrument.time == 2.0

        # points quotes have no constructor (the pair-dependent pip scale must stay
        # visible at the call site); the documented idiom is explicit arithmetic
        points = [10.0, 20.0]
        qps = FX.Outright.(eurusd, 1.10 .+ points ./ 10_000, [0.5, 1.0])
        @test qps[1].instrument.strike ≈ 1.101
    end

    @testset "implied zero-coupon quotes" begin
        # a non-flat domestic curve proves DF_d enters and cancels correctly
        usd_boot = fit(Spline.Linear(), ZCBYield.([0.04, 0.045, 0.05, 0.052], [1.0, 2.0, 5.0, 10.0]), Fit.Bootstrap())
        m_true = FX.Forwards(eurusd, S, usd_boot, eur)
        ts = [0.5, 1.0, 2.0, 5.0, 10.0]
        quotes = [FX.Outright(eurusd, forward(m_true, t), t) for t in ts]

        # inv under non-flat curves: the role swap plus spot inversion is exact
        @test forward(inv(m_true), 2.0) ≈ 1 / forward(m_true, 2.0)

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

        # a corrupt quote implying a non-positive discount factor fails at the source
        # rather than as NaN inside a downstream spline fit
        bad = Quote(-2.0, FX.Forward(eurusd, forward(m_true, 1.0), 1.0))
        @test_throws ArgumentError FX.implied_zcb_quotes([bad], S, usd_boot)
        @test_throws ArgumentError FX.__implied_foreign_quote(m_true, bad)

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

    @testset "textbook: covered interest parity (Hull)" begin
        audusd = FX.Pair(:AUD, :USD)
        # Hull, "Options, Futures, and Other Derivatives", the parameterization used
        # through the 9th edition (§5.10): spot 0.6200 USD per AUD and 2-year
        # continuously-compounded rates of 5% (AUD) and 7% (USD), so the 2-year
        # forward is 0.62·e^{(0.07−0.05)·2} = 0.6453.
        hull9 = FX.Forwards(audusd, 0.6200, Yield.Constant(Continuous(0.07)), Yield.Constant(Continuous(0.05)))
        @test forward(hull9, 2.0) ≈ 0.6453 atol = 5e-5
        # covered interest arbitrage: a forward struck below/above the CIP rate has
        # positive/negative value to the buyer of AUD
        @test pv(hull9, FX.Forward(audusd, 0.6300, 2.0)) > 0
        @test pv(hull9, FX.Forward(audusd, 0.6600, 2.0)) < 0

        # later editions restate the same example (Example 5.6): spot 0.7500 USD per
        # AUD, 2-year rates 3% (AUD) and 1% (USD) → forward 0.75·e^{(0.01−0.03)·2}
        # = 0.7206, and the two covered-interest-arbitrage strategies lock in riskless
        # profits of 21.87 and 55.79 USD at t = 2
        hull11 = FX.Forwards(audusd, 0.7500, Yield.Constant(Continuous(0.01)), Yield.Constant(Continuous(0.03)))
        @test forward(hull11, 2.0) ≈ 0.7206 atol = 5e-5
        # borrow 1,000 AUD (owing 1,061.84 at t=2) and buy them forward at 0.7000
        @test 1061.84 * (forward(hull11, 2.0) - 0.7000) ≈ 21.87 atol = 0.01
        # borrow 1,000 USD (1,415.79 AUD once grown at 3%) and sell forward at 0.7600
        @test 1415.79 * (0.7600 - forward(hull11, 2.0)) ≈ 55.79 atol = 0.01
        # the contract's time-0 value is the discounted locked-in profit
        @test 1061.84 * pv(hull11, FX.Forward(audusd, 0.7000, 2.0)) ≈ 21.87 * exp(-0.01 * 2) atol = 0.01
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

        # rolling the valuation date through the projection path: only cashflows at
        # or after cur_time remain, discounted from cur_time
        cur = 2.25
        manual = sum(cf.amount * forward(fx, cf.time) * discount(usd_boot, cur, cf.time) for cf in raw if cf.time >= cur)
        @test pv(usd_boot, Projection(FX.Converted(bond, "EURUSD"), store, CashflowProjection()), cur) ≈ manual

        # the identity pair converts at forward ≡ 1, letting generic multi-currency
        # code route domestic contracts through the same machinery
        idfx = FX.Forwards(FX.Pair(:USD, :USD), 1.0, usd_boot, usd_boot)
        @test forward(idfx, 2.0) == 1.0
        @test pv(usd_boot, Projection(FX.Converted(bond, "USDUSD"), Dict("USDUSD" => idfx), CashflowProjection())) ≈ pv(usd_boot, bond)

        @test maturity(FX.Converted(bond, "EURUSD")) == 5.0
        # a missing model key errors loudly
        @test_throws KeyError collect(Projection(FX.Converted(bond, "GBPUSD"), store, CashflowProjection()))
        # ...and so does a mis-keyed store: a yield curve where the FX model belongs
        # is named at the source instead of failing inside the transducer pipeline
        @test_throws ArgumentError collect(Projection(FX.Converted(bond, "EURUSD"), Dict("EURUSD" => usd_boot), CashflowProjection()))
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

    @testset "textbook: currency swap (Hull Examples 7.2 & 7.3, later editions)" begin
        # Later editions restate the §7.9 swap with updated rates: flat 2.5% (USD) and
        # 1.5% (JPY), continuously compounded; ¥110 = $1; receive 3% on ¥1,200M and
        # pay 4% on $10M annually for 3 more years. Example 7.2 values it as a
        # portfolio of FX forwards (Table forwards 0.009182 / 0.009275 / 0.009368);
        # Example 7.3 as two bonds: B_D = $10.4191M, B_F = ¥1,252.01M, so the swap is
        # worth 1,252.01/110 − 10.4191 ≈ $0.9629M (0.96279 unrounded; the book's
        # printed total carries rounding).
        jpyusd = FX.Pair(:JPY, :USD)
        usd25 = Yield.Constant(Continuous(0.025))
        jpy15 = Yield.Constant(Continuous(0.015))
        spot = 1 / 110
        fx = FX.Forwards(jpyusd, spot, usd25, jpy15)

        @test forward(fx, 1.0) ≈ 0.009182 atol = 5e-7
        @test forward(fx, 2.0) ≈ 0.009275 atol = 5e-7
        @test forward(fx, 3.0) ≈ 0.009368 atol = 5e-7

        # Example 7.2's year-1 row: net cashflow 36·F(1) − 0.4 = −0.0694 discounts
        # to −0.0677
        @test (36 * forward(fx, 1.0) - 0.4) * discount(usd25, 1.0) ≈ -0.0677 atol = 1e-4

        yen_leg = FX.Converted(Bond.Fixed(0.03, Periodic(1), 3) |> Map(cf -> cf * 1200.0), "JPYUSD")
        usd_leg = Bond.Fixed(0.04, Periodic(1), 3) |> Map(cf -> cf * -10.0)
        swap = Composite(yen_leg, usd_leg)
        p = Projection(swap, Dict("JPYUSD" => fx), CashflowProjection())
        @test pv(usd25, p) ≈ 0.9629 atol = 2e-4

        B_D = 10.0 * pv(usd25, Bond.Fixed(0.04, Periodic(1), 3))
        B_F = 1200.0 * pv(jpy15, Bond.Fixed(0.03, Periodic(1), 3))
        @test B_D ≈ 10.4191 atol = 1e-4
        @test B_F ≈ 1252.01 atol = 1e-2
        @test pv(usd25, p) ≈ B_F * spot - B_D
    end

    @testset "FX.ParBasisSwap long-end basis calibration" begin
        sofr = fit(Spline.Linear(), ZCBYield.([0.045, 0.047, 0.048, 0.05], [1.0, 2.0, 5.0, 10.0]), Fit.Bootstrap())
        estr = fit(Spline.Linear(), ZCBYield.([0.028, 0.03, 0.031, 0.032], [1.0, 2.0, 5.0, 10.0]), Fit.Bootstrap())

        @testset "hand-derived: flat projection + flat spread" begin
            # with an annually-paying leg off a flat 3% Periodic(1) projection curve and
            # a −20bp spread, every coupon is 0.028, so the implied CSA curve is exactly
            # flat 2.8% Periodic(1): DF_f(t) = 1.028^(−t)
            ref = Yield.Constant(Periodic(0.03, 1))
            q1 = FX.ParBasisSwap(eurusd, -0.002, 1.0; reference = ref, frequency = Periodic(1))
            @test q1.price == 1.0
            @test maturity(q1) == 1.0
            @test only(q1.instrument.cashflows).amount ≈ 1.028

            tenors = [1.0, 2.0, 5.0]
            qs = [FX.ParBasisSwap(eurusd, -0.002, T; reference = ref, frequency = Periodic(1)) for T in tenors]
            m_fit = fit(FX.Forwards(eurusd, S, usd, Spline.Linear()), qs, Fit.Bootstrap())
            @test all(isapprox(discount(m_fit.foreign, T), 1.028^-T; atol = 1e-8) for T in tenors)

            # at-market legs price to par (in base-currency units) under the true model
            csa28 = Yield.Constant(Periodic(0.028, 1))
            m_true = FX.Forwards(eurusd, S, usd, csa28)
            @test all(pv(m_true, q.instrument) ≈ 1.0 for q in qs)
            # and the leg composes with FX.Converted: domestic value = spot × foreign value
            leg = qs[2].instrument
            store = Dict("EURUSD" => m_true)
            @test pv(usd, Projection(FX.Converted(leg, "EURUSD"), store, CashflowProjection())) ≈ S * pv(csa28, leg)
        end

        @testset "materialization matches Bond.Floating projection" begin
            b = -0.0015
            q = FX.ParBasisSwap(eurusd, b, 3.0; reference = estr)
            frn = collect(Projection(Bond.Floating(b, Periodic(4), 3.0, "R"), Dict("R" => estr), CashflowProjection()))
            legcfs = q.instrument.cashflows
            @test length(legcfs) == length(frn) == 12
            @test all(legcfs[i].time == frn[i].time for i in eachindex(frn))
            @test all(legcfs[i].amount ≈ frn[i].amount for i in eachindex(frn))
        end

        @testset "mixed short-end forwards + long-end basis swaps" begin
            # truth: the CSA curve is the projection curve plus a flat −18bp basis
            basis = Yield.Constant(Continuous(-0.0018))
            csa_true = estr + basis
            m_true = FX.Forwards(eurusd, S, sofr, csa_true)
            freq = Periodic(4)
            f = freq.frequency

            # the par spread each tenor must carry for the swap to be at-market on the
            # truth curve: b = (1 − DF(T) − Σ fᵢ/f·DF(tᵢ)) / (Σ DF(tᵢ)/f)
            tenors = [2.0, 5.0, 10.0]
            spreads = map(tenors) do T
                ts = Bond.coupon_times(T, f)
                dfs = [discount(csa_true, t) for t in ts]
                fwds = [rate(freq(forward(estr, t - 1 / f, t))) for t in ts]
                (1 - dfs[end] - sum(fwds[i] / f * dfs[i] for i in eachindex(ts))) / (sum(dfs) / f)
            end
            # sanity: implied par spreads sit near the flat continuous basis
            @test all(abs(b - -0.0018) < 2e-4 for b in spreads)

            fx_qs = [FX.Outright(eurusd, forward(m_true, t), t) for t in [0.25, 0.5, 1.0]]
            bs_qs = [FX.ParBasisSwap(eurusd, spreads[i], tenors[i]; reference = estr) for i in eachindex(tenors)]
            # local (linear) interpolation keeps the sequential bootstrap exact for
            # coupon-bearing quotes; a global spline (e.g. Cubic) reshapes earlier
            # segments as later knots are added, drifting solved quotes off par
            m_fit = fit(FX.Forwards(eurusd, S, sofr, Spline.Linear()), vcat(fx_qs, bs_qs), Fit.Bootstrap())

            # both quote families reprice on the fitted model
            @test all(abs(pv(m_fit, q.instrument)) < 1e-9 for q in fx_qs)
            @test all(abs(pv(m_fit, q.instrument) - 1.0) < 1e-9 for q in bs_qs)

            # independent cross-check through the projection machinery: the actual
            # constant-notional basis swap — a floating EUR leg converted at the fitted
            # forwards against a spot-scaled floating USD leg — prices to zero
            T5 = tenors[2]
            swap = Composite(
                FX.Converted(Bond.Floating(spreads[2], freq, T5, "ESTR"), "EURUSD"),
                Bond.Floating(0.0, freq, T5, "SOFR") |> Map(cf -> cf * -S),
            )
            p = Projection(swap, Dict("EURUSD" => m_fit, "ESTR" => estr, "SOFR" => sofr), CashflowProjection())
            @test abs(pv(sofr, p)) < 1e-6
        end

        @testset "parametric foreign curve via generic fit" begin
            ref = Yield.Constant(Periodic(0.03, 1))
            qs = [FX.ParBasisSwap(eurusd, -0.002, T; reference = ref, frequency = Periodic(1)) for T in [1.0, 2.0, 3.0]]
            m_fit = fit(FX.Forwards(eurusd, S, usd, Yield.Constant()), qs)
            @test discount(m_fit.foreign, 1.0) ≈ 1 / 1.028 atol = 1e-6
            @test discount(m_fit.foreign, 3.0) ≈ 1.028^-3 atol = 1e-6
        end

        @testset "guards" begin
            gbpusd = FX.Pair(:GBP, :USD)
            q_gbp = FX.ParBasisSwap(gbpusd, -0.001, 2.0; reference = estr)
            m0 = FX.Forwards(eurusd, S, sofr, Spline.Cubic())
            @test_throws ArgumentError fit(m0, [q_gbp], Fit.Bootstrap())
            m_e = FX.Forwards(eurusd, S, sofr, estr)
            @test_throws ArgumentError pv(m_e, q_gbp.instrument)
            @test_throws ArgumentError pv(m_e, q_gbp.instrument, 0.5)
            # an outright and a basis swap at the same maturity are refused by the
            # bootstrap's distinct-maturities check rather than silently blended
            clash = [FX.Outright(eurusd, 1.12, 2.0), FX.ParBasisSwap(eurusd, -0.001, 2.0; reference = estr)]
            @test_throws ArgumentError fit(m0, clash, Fit.Bootstrap())
        end
    end
end
