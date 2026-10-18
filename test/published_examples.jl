# Worked examples from the literature, reproduced at their published precision. Each testset cites
# its source; long tables are in test/reference/, transcribed as printed.

@testset "published examples" begin
    read_table(file) = [
        parse.(Float64, split(line, ',')) for line in eachline(joinpath(@__DIR__, "reference", file))
            if !startswith(line, '#') && !startswith(line, "month")
    ]

    @testset "mortgage pools and amortization" begin
        # the PSA benchmark's CPR at a loan age in months
        psa(speed, age) = speed / 100 * 0.06 * min(age, 30) / 30
        runoff(loan, cpr) = collect(Projection(Bond.Pool(loan, "speed"), Models(Yield.Constant(0.05), Dict("speed" => cpr)), AssetProjection()))
        smm(row) = row.prepayment / (row.opening - row.principal)

        @testset "Fabozzi (2016), Handbook of Mortgage-Backed Securities, Ch. 3" begin
            # Table 3.1 (p. 94): the SMM of each month's PSA CPR, to six decimals, at 100, 90, 165 and 300
            # PSA; the last value is months 31–360 (a pool's last month has nothing left to prepay)
            table = Dict(
                100 => [
                    0.000167, 0.000334, 0.000501, 0.000669, 0.000837, 0.001006, 0.001174, 0.001343, 0.001513, 0.001682,
                    0.001852, 0.002022, 0.002193, 0.002364, 0.002535, 0.002707, 0.002878, 0.003051, 0.003223, 0.003396,
                    0.003569, 0.003743, 0.003917, 0.004091, 0.004265, 0.00444, 0.004615, 0.004791, 0.004967, 0.005143, 0.005143,
                ],
                90 => [
                    0.00015, 0.0003, 0.000451, 0.000602, 0.000753, 0.000904, 0.001056, 0.001208, 0.00136, 0.001513,
                    0.001665, 0.001818, 0.001971, 0.002125, 0.002278, 0.002432, 0.002586, 0.002741, 0.002896, 0.003051,
                    0.003206, 0.003361, 0.003517, 0.003673, 0.00383, 0.003986, 0.004143, 0.0043, 0.004458, 0.004615, 0.004615,
                ],
                165 => [
                    0.000275, 0.000552, 0.000829, 0.001107, 0.001386, 0.001665, 0.001946, 0.002227, 0.002509, 0.002792,
                    0.003077, 0.003361, 0.003647, 0.003934, 0.004222, 0.00451, 0.0048, 0.00509, 0.005381, 0.005674,
                    0.005967, 0.006261, 0.006556, 0.006852, 0.00715, 0.007448, 0.007747, 0.008047, 0.008348, 0.00865, 0.00865,
                ],
                300 => [
                    0.000501, 0.001006, 0.001513, 0.002022, 0.002535, 0.003051, 0.003569, 0.004091, 0.004615, 0.005143,
                    0.005674, 0.006208, 0.006745, 0.007285, 0.007828, 0.008375, 0.008925, 0.009479, 0.010036, 0.010596,
                    0.01116, 0.011728, 0.012299, 0.012873, 0.013452, 0.014034, 0.01462, 0.01521, 0.015804, 0.016402, 0.016402,
                ],
            )
            for (speed, smms) in table
                rows = runoff(Bond.LevelPayment(0.08125, 12, 30), Prepayment.CPR([psa(speed, t) for t in 1:360]))
                @test round.(smm.(rows[[1:31; 359]]); digits = 6) == [smms; last(smms)]
            end

            # Tables 3.2 and 3.3 (pp. 98, 100): a $400 million passthrough, gross WAC 8.125%, net coupon
            # 7.5%, WAM 357 months, seasoned three months, so its month t is loan age t + 3. The pool
            # pays interest at the gross rate; the servicing spread is outside the asset, so net interest
            # is the opening balance at 7.5%.
            loan = Bond.LevelPayment(0.08125, 12, 357 / 12)
            passthrough(speed) = runoff(loan, Prepayment.CPR([psa(speed, t + 3) for t in 1:357]))
            net(row) = row.opening * 0.075 / 12
            rows = passthrough(165)
            @test length(rows) == 357
            for (month, balance, s, payment, interest, scheduled, prepaid, principal, cashflow) in read_table("fabozzi_hmbs_table3_2.csv")
                r = rows[Int(month)]
                @test round.(4.0e8 .* [r.opening, r.interest + r.principal, net(r), r.principal, r.prepayment, r.principal + r.prepayment, net(r) + r.principal + r.prepayment]) ==
                    [balance, payment, interest, scheduled, prepaid, principal, cashflow]
                month < 357 && @test round(smm(r); digits = 5) == s
            end
            @test iszero(last(rows).prepayment) && last(rows).closing == 0
            rows = passthrough(100)
            for (month, prepaid, cashflow) in read_table("fabozzi_hmbs_table3_3.csv")
                r = rows[Int(month)]
                @test round.(4.0e8 .* [r.prepayment, net(r) + r.principal + r.prepayment]) == [prepaid, cashflow]
            end
        end

        @testset "Fabozzi and Kothari (2008), Introduction to Securitization, Table 3.3" begin
            # a $660 million pool, gross WAC 6%, WAM 358 months, seasoned two months: total principal
            # at 100 and 250 PSA
            loan = Bond.LevelPayment(0.06, 12, 358 / 12)
            total(speed) = [6.6e8 * (r.principal + r.prepayment) for r in runoff(loan, Prepayment.CPR([psa(speed, t + 2) for t in 1:358]))]
            slow, fast = total(100), total(250)
            for (month, at100, at250) in read_table("fabozzi_kothari_table3_3.csv")
                @test round.([slow[Int(month)], fast[Int(month)]]) == [at100, at250]
            end
        end

        @testset "SIFMA Uniform Practices, Ch. SF, Standard Formulas (1999)" begin
            # SF-4: the first month of a 9.5% 360-month passthrough at 150 PSA (net coupon 9%), per unit
            # of par: scheduled amortization, gross interest, prepayment, principal and cash flow
            r = first(runoff(Bond.LevelPayment(0.095, 12, 30), Prepayment.CPR([psa(150, t) for t in 1:360])))
            @test round.([r.principal, r.interest, r.prepayment, r.principal + r.prepayment, r.principal + r.prepayment + 0.09 / 12]; digits = 8) ==
                [0.00049188, 0.00791667, 0.00025022, 0.0007421, 0.0082421]

            # SF-7 and SF-8: a 9.5% pool amortizing over 359 months reports factors 0.85150625 and
            # 0.84732282 for consecutive months, with 344 and then 343 months remaining
            s = Bond.__schedule(Bond.LevelPayment(0.095, 12, 359 / 12))
            @test round.([s[15].closing, s[16].closing]; digits = 8) == [0.992133, 0.99157471]
            F1, F2 = 0.85150625, 0.84732282
            scheduled = F1 * s[16].closing / s[15].closing
            @test round.([scheduled, F1 - scheduled, scheduled - F2]; digits = 8) == [0.85102709, 0.00047916, 0.00370427]
            implied = 1 - F2 / scheduled
            @test round(100 * implied; digits = 6) == 0.43527
            @test round(100 * (1 - (1 - implied)^12); digits = 4) == 5.1
            # 150 PSA in the pool's 17th month reproduces the second factor
            rows = runoff(Bond.LevelPayment(0.095, 12, 359 / 12), Prepayment.CPR([t == 16 ? psa(150, 17) : 0.0 for t in 1:359]))
            @test round(F1 * rows[16].closing / rows[16].opening; digits = 8) == F2
        end

        @testset "level-payment schedules" begin
            # Carpenter, J., "Debt Instruments and Markets: The Mortgage Market", NYU Stern lecture notes,
            # slides 5, 6 and 10: $100,000 at 7.25% over 360 months (opening balance, payment, interest,
            # principal), and 100 at 5.5% compounded semiannually over two years, as a pool in which half
            # the loans prepay at 0.5 and a quarter of the rest at 1.5
            rows = runoff(Bond.LevelPayment(0.0725, 12, 30), Prepayment.CPR(0.0))
            @test round.(1.0e5 .* [rows[1].opening, rows[1].interest + rows[1].principal]; digits = 2) == [100_000.0, 682.18]
            for (k, opening, interest, principal) in [(1, 100_000.0, 604.17, 78.01), (2, 99_921.99, 603.7, 78.48), (3, 99_843.51, 603.22, 78.96), (4, 99_764.55, 602.74, 79.43), (360, 678.08, 4.1, 678.08)]
                @test round.(1.0e5 .* [rows[k].opening, rows[k].interest, rows[k].principal]; digits = 2) == [opening, interest, principal]
            end
            semiannual = Bond.LevelPayment(0.055, 2, 2)
            rows = runoff(semiannual, Prepayment.CPR(0.0))
            @test [round.(100 .* [r.opening, r.interest + r.principal, r.interest, r.principal, r.closing]; digits = 2) for r in rows] == [
                [100.0, 26.74, 2.75, 23.99, 76.01],
                [76.01, 26.74, 2.09, 24.65, 51.36],
                [51.36, 26.74, 1.41, 25.33, 26.03],
                [26.03, 26.74, 0.72, 26.03, 0.0],
            ]
            # the share prepaying over a half-year is 1 - (1 - CPR)^(1/2)
            rows = runoff(semiannual, Prepayment.CPR([0.75, 0.0, 0.4375, 0.0]))
            @test [round.(100 .* [r.opening, r.interest + r.principal, r.interest, r.principal, r.prepayment, r.closing]; digits = 2) for r in rows] == [
                [100.0, 26.74, 2.75, 23.99, 38.0, 38.0],
                [38.0, 13.37, 1.05, 12.33, 0.0, 25.68],
                [25.68, 13.37, 0.71, 12.66, 3.25, 9.76],
                [9.76, 10.03, 0.27, 9.76, 0.0, 0.0],
            ]

            # Berliner, Quinones and Bhattacharya, Ch. 1 of Fabozzi (2016), p. 13: $200,000 at 4.5% over
            # 30 years pays $1,013 a month: $750 interest and $263 principal in month 1, $685 and $328 in
            # month 60, and "both roughly $507" in month 176
            rows = runoff(Bond.LevelPayment(0.045, 12, 30), Prepayment.CPR(0.0))
            @test round(2.0e5 * (rows[1].interest + rows[1].principal)) == 1013
            @test round.(2.0e5 .* [rows[1].interest, rows[1].principal, rows[60].interest, rows[60].principal]) == [750, 263, 685, 328]
            @test all(x -> abs(2.0e5 * x - 507) < 1, (rows[176].interest, rows[176].principal))
        end
    end

    @testset "credit" begin
        risky(asset) = Credit.Defaultable(asset, "issuer")
        on(curve, hazard, recovery) = Models(curve, Dict("issuer" => Credit.Assumptions(hazard; recovery)))

        @testset "Hull (2017), Options, Futures, and Other Derivatives, 10th ed., Tables 25.1–25.3" begin
            # CDS valuation: hazard rate 2%, recovery 40%, risk-free 5% continuously compounded, five
            # annual premiums in arrears
            curve, hazard = Yield.Constant(Continuous(0.05)), Credit.ConstantHazard(0.02)
            # Table 25.1: survival and unconditional default probabilities by year
            @test round.([Credit.survival(hazard, t) for t in 1:5]; digits = 4) == [0.9802, 0.9608, 0.9418, 0.9231, 0.9048]
            @test round.(diff([Credit.default_probability(hazard, t) for t in 0:5]); digits = 4) == [0.0198, 0.0194, 0.019, 0.0186, 0.0183]
            # Table 25.2: the premium leg per unit of spread is the survival-weighted annuity, a
            # zero-recovery bond paying 1 a year less its principal
            zero = Bond.Fixed(0.0, Periodic(1), 5)
            pv(asset, recovery) = present_value(on(curve, hazard, recovery), risky(asset))
            @test round(pv(Bond.Fixed(1.0, Periodic(1), 5), 0.0) - pv(zero, 0.0); digits = 4) == 4.0728
            # Table 25.3 and p. 576: the protection leg, (1 - R) of par at default, and the binary CDS's,
            # par at default, are the recovery legs of a zero at recoveries 0.6 and 1. Hull's defaults
            # happen at mid-year; continuous-time default agrees at the published four decimals.
            @test round(pv(zero, 0.6) - pv(zero, 0.0); digits = 4) == 0.0506
            @test round(pv(zero, 1.0) - pv(zero, 0.0); digits = 4) == 0.0844
        end

        @testset "Hull (2002), Options, Futures, and Other Derivatives, 5th ed., Ch. 26" begin
            # Zero-coupon corporate yields of 5.25%, 5.50%, 5.70%, 5.85% and 5.95% for 1 to 5 years
            # over a flat 5% risk-free yield (continuously compounded) imply cumulative expected losses
            # of 0.2497%, 0.9950%, 2.0781%, 3.3428% and 4.6390%. With no recovery the expected loss is
            # the default probability, so bootstrapping a hazard curve from the bonds' prices gives it.
            yields = [0.0525, 0.055, 0.057, 0.0585, 0.0595]
            bonds = [risky(Bond.Fixed(0.0, Periodic(1), t)) for t in 1:5]
            curve = Yield.Constant(Continuous(0.05))
            fitted = fit(Credit.Assumptions(Credit.PiecewiseHazard(); recovery = 0.0), Quote.(exp.(-yields .* (1:5)), bonds), Fit.Bootstrap(); discount = curve)
            @test round.([100 * Credit.default_probability(fitted.hazard, t) for t in 1:5]; digits = 4) == [0.2497, 0.995, 2.0781, 3.3428, 4.639]
        end
    end

    @testset "short-rate models and trees" begin
        @testset "Hull and White (1996), Using Hull–White interest rate trees" begin
            # Journal of Derivatives 3(3); also Hull, Options, Futures, and Other Derivatives, 11th ed.,
            # Example 32.1. Table 1: the DM zero curve of July 8, 1994, continuously compounded, with
            # maturities in days of 1/365 year and linear interpolation; a = 0.1, σ = 0.01.
            days = [3, 31, 62, 94, 185, 367, 731, 1096, 1461, 1826, 2194, 2558, 2922, 3287, 3653]
            rates = [5.01772, 4.98284, 4.97234, 4.96157, 4.99058, 5.09389, 5.79733, 6.30595, 6.73464, 6.94816, 7.08807, 7.27527, 7.30852, 7.3979, 7.49015]
            curve = ZeroRateCurve(rates ./ 100, days ./ 365, Spline.Linear())
            hw = ShortRate.HullWhite(0.1, 0.01, curve)
            # Table 3: zero rates and discount factors at the tree's times and at 9 years, and the
            # instantaneous forward rate and α(t) = f(0, t) + σ²(1 - e^{-at})²/2a² (equation 8)
            ts = [1, 2, 3, 4, 9]
            @test round.([100 * rate(zero(curve, t)) for t in ts]; digits = 6) == [5.092755, 5.795397, 6.304557, 6.733466, 7.39741]
            @test round.([discount(curve, t) for t in ts]; digits = 6) == [0.950348, 0.890557, 0.827673, 0.763885, 0.513879]
            @test round.([100 * rate(Yield.instantaneous_forward(curve, t)) for t in 0:3]; digits = 6) == [5.01772, 5.299942, 7.206143, 7.830417]
            @test round.([100 * FinanceModels.__hw_alpha(hw, t) for t in 0:3]; digits = 6) == [5.01772, 5.30447, 7.222572, 7.864004]
            # Tables 2 and 3: the three-step tree (Δt = 1 year): the spacing Δr, the branching
            # probabilities (down, middle, up) of the nodes j = -1, 0, 1, and the α by which each step's
            # rates are raised to fit the discount factors. The paper truncates the tree at |j| = 2 with
            # non-standard branching; this lattice branches to the nearest nodes everywhere.
            lattice = FinanceModels.__build_lattice(hw, [0.0, 1.0, 2.0, 3.0, 4.0], 0.0)
            @test round(lattice.dx[2]; digits = 9) == 0.016489508
            @test [round.(p; digits = 6) for p in lattice.p[2]] == [(0.123613, 0.657611, 0.218776), (0.166667, 0.666667, 0.166667), (0.218776, 0.657611, 0.123613)]
            @test round.(100 .* lattice.α; digits = 5) == [5.09275, 6.50257, 7.33932, 8.05381]

            # Table 8: a 3-year (3 × 365 days) European put, strike 63, on a zero paying 100 at 9 years,
            # analytic price 1.8093
            @test round(100 * present_value(hw, Option.ZCBPut(3.0, 9.0, 0.63)); digits = 4) == 1.8093
            # On the lattice, the put is a zero with one put date less the zero. The paper's tree (Table
            # 8: 1.8491, 1.8179, 1.8060, 1.8128, 1.8089 at 10, 30, 50, 100, 200 steps) prices the bond
            # analytically at expiry and truncates at |j| = 2; this lattice equals QuantLib's on the same
            # grid (test/reference/quantlib_callable.py) and oscillates toward the analytic price as the
            # paper's does
            putable = Bond.Putable(Bond.Amortizing(0.0, [3.0, 9.0], [0.0, 1.0]), [3.0], 0.63)
            tree(n) = 100 * (present_value(ShortRate.HullWhiteTree(hw; max_step = 3 / n), putable) - discount(curve, 9.0))
            values = Dict(n => tree(n) for n in (10, 30, 50, 100, 200))
            for (n, quantlib) in ((10, 1.851067584587), (30, 1.818511838558), (50, 1.806380533112), (100, 1.813054215094), (200, 1.809025970671))
                @test values[n] ≈ quantlib atol = 1.0e-10
            end
            # at 200 steps the error is below that of the paper's tree (1.8089 against 1.8093)
            @test abs(values[200] - 100 * present_value(hw, Option.ZCBPut(3.0, 9.0, 0.63))) < 4.0e-4
        end

        @testset "Hull, Technical Note 15, Valuing Options on Coupon-Bearing Bonds in a One-Factor Interest Rate Model" begin
            # Jamshidian's decomposition: under Vasicek with a = 0.1, b = 0.1, σ = 0.02 and r₀ = 10%, a
            # 3-year European put, strike 98, on a bond paying 5 every six months and 100 at 5 years
            m = ShortRate.Vasicek(0.1, 0.1, 0.02, Continuous(0.1))
            ts, cs = [3.5, 4.0, 4.5, 5.0], [5, 5, 5, 105]
            # at 3 years the bond is Σ cᵢ A(3, tᵢ) e^{-B(3, tᵢ) r}
            @test round.([discount(m, 3.0, t, 0.0) for t in ts]; digits = 4) == [0.9988, 0.9952, 0.9895, 0.9819]
            @test round.([log(discount(m, 3.0, t, 0.0) / discount(m, 3.0, t, 1.0)) for t in ts]; digits = 4) == [0.4877, 0.9516, 1.3929, 1.8127]
            @test round.([discount(m, 3.0), discount(m, 5.0)]; digits = 4) == [0.7419, 0.6101]
            # the rate at which the bond is worth the strike, and the zero strikes it implies
            excess(r) = sum(c * discount(m, 3.0, t, r) for (c, t) in zip(cs, ts)) - 98
            lo, hi = 0.0, 0.3
            while hi - lo > 1.0e-14
                mid = (lo + hi) / 2
                excess(mid) > 0 ? (lo = mid) : (hi = mid)
            end
            @test round(lo; digits = 5) == 0.10952
            strikes = [c * discount(m, 3.0, t, lo) for (c, t) in zip(cs, ts)]
            @test round.(strikes; digits = 3) == [4.734, 4.484, 4.248, 84.535]
            # the zero-bond puts on the 4- and 4.5-year coupons. The note's other two, 0.0125 and 0.8085
            # (total 0.8752), carry its four-figure σ_P and h; unrounded they are 0.012449 and 0.808418
            # (total 0.875126)
            puts = [c * present_value(m, Option.ZCBPut(3.0, t, k / c)) for (c, t, k) in zip(cs, ts, strikes)]
            @test round.(puts[2:3]; digits = 4) == [0.0228, 0.0314]
        end
    end
end
