using Random

@testset "Stochastic Models" begin

    @testset "Vasicek" begin
        # Vasicek parameters from Brigo & Mercurio example
        a, b, σ, r0 = 0.136, 0.0168, 0.0119, 0.01
        v = ShortRate.Vasicek(a, b, σ, Continuous(r0))

        @testset "discount (closed-form ZCB)" begin
            # Verify against hand-computed values
            # B(T) = (1 - exp(-a*T))/a
            # lnA = (B-T)*(a^2*b - 0.5*σ^2)/a^2 - σ^2*B^2/(4a)
            # P = exp(lnA - B*r0)
            T = 5.0
            B = (1 - exp(-a * T)) / a
            lnA = (B - T) * (a^2 * b - 0.5 * σ^2) / a^2 - σ^2 * B^2 / (4a)
            expected = exp(lnA - B * r0)
            @test discount(v, T) ≈ expected

            # discount at time 0 should be 1
            @test discount(v, 0.0) ≈ 1.0

            # discount should be positive and ≤ 1 for reasonable rates
            @test 0 < discount(v, 10.0) < 1
        end

        @testset "zero, forward, par derived from discount" begin
            z = zero(v, 5.0)
            @test z isa FinanceCore.Rate
            @test FinanceCore.rate(z) > 0

            f = forward(v, 2.0, 3.0)
            @test f isa FinanceCore.Rate

            p = par(v, 5.0)
            @test p isa FinanceCore.Rate
        end

        @testset "present_value with fixed bond" begin
            bond = Bond.Fixed(0.05, Periodic(2), 10)
            pv_val = present_value(v, bond)
            @test pv_val > 0

            # A bond paying 5% coupons discounted at ~1-2% should be worth more than par
            @test pv_val > 1.0
        end

        @testset "simulate" begin
            rng = Random.MersenneTwister(42)
            scenarios = simulate(v; n_scenarios = 50, timestep = 1 / 12, horizon = 10.0, rng)
            @test length(scenarios) == 50
            @test scenarios[1] isa RatePath

            # Each scenario should function as a yield model
            @test discount(scenarios[1], 5.0) > 0
            @test discount(scenarios[1], 0.0) ≈ 1.0

            # PV per scenario should work
            bond = Bond.Fixed(0.05, Periodic(2), 5)
            pv1 = present_value(scenarios[1], bond)
            @test pv1 > 0
        end

        @testset "pv_mc" begin
            rng = Random.MersenneTwister(123)
            bond = Bond.Fixed(0.05, Periodic(2), 5)
            mc_pv = pv_mc(v, bond; n_scenarios = 500, timestep = 1 / 12, rng)
            analytical_pv = present_value(v, bond)

            # MC estimate should be in the ballpark of analytical
            @test mc_pv ≈ analytical_pv rtol = 0.05
        end

        @testset "mean reversion" begin
            # With high mean reversion, long-term rate should approach b
            v_strong = ShortRate.Vasicek(2.0, 0.05, 0.01, Continuous(0.10))
            z_long = zero(v_strong, 30.0)
            @test abs(FinanceCore.rate(z_long) - 0.05) < 0.01
        end
    end

    @testset "CoxIngersollRoss" begin
        a, b, σ, r0 = 0.3, 0.05, 0.1, 0.03
        cir = ShortRate.CoxIngersollRoss(a, b, σ, Continuous(r0))

        @testset "discount (closed-form ZCB)" begin
            # Verify against CIR formula
            T = 5.0
            γ = sqrt(a^2 + 2σ^2)
            expγT = exp(γ * T)
            denom = (γ + a) * (expγT - 1) + 2γ
            B = 2(expγT - 1) / denom
            A = (2γ * exp((a + γ) * T / 2) / denom)^(2a * b / σ^2)
            expected = A * exp(-B * r0)
            @test discount(cir, T) ≈ expected

            @test discount(cir, 0.0) ≈ 1.0
            @test 0 < discount(cir, 10.0) < 1
        end

        @testset "zero, forward, par" begin
            @test zero(cir, 5.0) isa FinanceCore.Rate
            @test forward(cir, 2.0, 3.0) isa FinanceCore.Rate
            @test par(cir, 5.0) isa FinanceCore.Rate
        end

        @testset "present_value with fixed bond" begin
            bond = Bond.Fixed(0.04, Periodic(2), 10)
            pv_val = present_value(cir, bond)
            @test pv_val > 0
        end

        @testset "simulate" begin
            rng = Random.MersenneTwister(42)
            scenarios = simulate(cir; n_scenarios = 50, timestep = 1 / 12, horizon = 10.0, rng)
            @test length(scenarios) == 50

            # CIR rates should stay non-negative (with Feller condition: 2ab > σ²)
            # We just check the paths produce valid discounts
            @test discount(scenarios[1], 5.0) > 0
        end

        @testset "pv_mc" begin
            rng = Random.MersenneTwister(456)
            bond = Bond.Fixed(0.04, Periodic(2), 5)
            mc_pv = pv_mc(cir, bond; n_scenarios = 500, timestep = 1 / 12, rng)
            analytical_pv = present_value(cir, bond)
            @test mc_pv ≈ analytical_pv rtol = 0.1
        end
    end

    @testset "HullWhite" begin
        curve = Yield.Constant(0.03)
        hw = ShortRate.HullWhite(0.1, 0.01, curve)

        @testset "discount matches initial curve" begin
            # Hull-White should match the initial term structure
            for T in [1.0, 5.0, 10.0, 20.0]
                @test discount(hw, T) ≈ discount(curve, T)
            end
        end

        @testset "zero, forward, par" begin
            @test zero(hw, 5.0) isa FinanceCore.Rate
            @test forward(hw, 2.0, 3.0) isa FinanceCore.Rate
            @test par(hw, 5.0) isa FinanceCore.Rate
        end

        @testset "present_value" begin
            bond = Bond.Fixed(0.04, Periodic(2), 10)
            pv_hw = present_value(hw, bond)
            pv_curve = present_value(curve, bond)
            @test pv_hw ≈ pv_curve
        end

        @testset "simulate" begin
            rng = Random.MersenneTwister(42)
            scenarios = simulate(hw; n_scenarios = 50, timestep = 1 / 12, horizon = 10.0, rng)
            @test length(scenarios) == 50
            @test discount(scenarios[1], 0.0) ≈ 1.0
        end
    end

    @testset "fit" begin
        # Fit Vasicek to ZCB yields
        quotes = ZCBYield.([0.02, 0.025, 0.03], [1, 5, 10])
        v0 = ShortRate.Vasicek(0.1, 0.02, 0.01, Continuous(0.01))
        v_fit = fit(v0, quotes)

        @test v_fit isa ShortRate.Vasicek
        # The fitted model should reprice the quotes reasonably well
        for q in quotes
            pv_fit = present_value(v_fit, q.instrument)
            @test pv_fit ≈ q.price rtol = 0.01
        end
    end

    @testset "Hull-White ZCB options" begin
        # Reference: Hull "Options, Futures, and Other Derivatives" Technical Note 31
        # a=0.08, σ=0.01, flat 10% continuous, T=1, S=5, L=100, K=68
        hw_hull = ShortRate.HullWhite(0.08, 0.01, Yield.Constant(Continuous(0.10)))

        @testset "Hull textbook example" begin
            call_price = present_value(hw_hull, Option.ZCBCall(1.0, 5.0, 0.68))
            # Hull gives call ≈ 0.00439 per unit (0.439 per 100 face)
            @test call_price ≈ 0.00439 atol = 0.0005
        end

        @testset "tf-quant-finance reference (a=0.03, σ=0.02, flat 1%)" begin
            # Reference: google/tf-quant-finance zero_coupon_bond_option_test.py
            hw_tf = ShortRate.HullWhite(0.03, 0.02, Yield.Constant(Continuous(0.01)))

            # T=1, S=5, K = P(0,5)/P(0,1) (ATM forward)
            K1 = exp(-0.01 * 5) / exp(-0.01 * 1)
            call1 = present_value(hw_tf, Option.ZCBCall(1.0, 5.0, K1))
            @test call1 ≈ 0.02817777 atol = 1e-5

            # T=2, S=4, K = P(0,4)/P(0,2) (ATM forward)
            K2 = exp(-0.01 * 4) / exp(-0.01 * 2)
            call2 = present_value(hw_tf, Option.ZCBCall(2.0, 4.0, K2))
            @test call2 ≈ 0.02042677 atol = 1e-5
        end

        @testset "put-call parity" begin
            # Call - Put = P(0,S) - K·P(0,T) for ZCB options
            hw_pcp = ShortRate.HullWhite(0.1, 0.015, Yield.Constant(Continuous(0.05)))
            T, S, K = 1.0, 5.0, 0.75
            call = present_value(hw_pcp, Option.ZCBCall(T, S, K))
            put  = present_value(hw_pcp, Option.ZCBPut(T, S, K))
            P0S = discount(hw_pcp, S)
            P0T = discount(hw_pcp, T)
            @test call - put ≈ P0S - K * P0T atol = 1e-10
        end
    end

    @testset "Hull-White caps and floors" begin
        # Reference: google/tf-quant-finance cap_floor_test.py
        # a=0.03, σ=0.02, flat 1% continuous, quarterly, 1-year cap, notional=1 (per unit)
        hw_cap = ShortRate.HullWhite(0.03, 0.02, Yield.Constant(Continuous(0.01)))

        @testset "tf-quant-finance cap reference values" begin
            # Strike 1%, quarterly, 1-year cap
            # tf-quant-finance expected: 0.4072088281493774 for notional=100
            # Our Cap is per unit notional, so expected = 0.004072088...
            cap_1pct = present_value(hw_cap, Option.Cap(0.01, 4, 1.0))
            @test cap_1pct * 100 ≈ 0.40720883 atol = 0.001

            # Strike 2%
            cap_2pct = present_value(hw_cap, Option.Cap(0.02, 4, 1.0))
            @test cap_2pct * 100 ≈ 0.14283513 atol = 0.001

            # Strike 3%
            cap_3pct = present_value(hw_cap, Option.Cap(0.03, 4, 1.0))
            @test cap_3pct * 100 ≈ 0.03980642 atol = 0.001
        end

        @testset "cap-floor parity" begin
            # Cap(K) - Floor(K) = forward swap value
            # = P(0,τ) - P(0,Nτ) - K·τ·Σ_{i=2}^{N} P(0,iτ)
            hw_cf = ShortRate.HullWhite(0.1, 0.015, Yield.Constant(Continuous(0.05)))
            strike = 0.05
            freq = 2
            mat = 5.0
            τ = 1.0 / freq
            n_periods = round(Int, mat * freq)

            cap  = present_value(hw_cf, Option.Cap(strike, freq, mat))
            flr  = present_value(hw_cf, Option.Floor(strike, freq, mat))

            # Compute the forward swap value from the curve
            annuity = sum(discount(hw_cf, i * τ) for i in 2:n_periods)
            swap_value = discount(hw_cf, τ) - discount(hw_cf, n_periods * τ) - strike * τ * annuity
            @test cap - flr ≈ swap_value atol = 1e-8
        end

        @testset "cap increases with volatility" begin
            hw_lo = ShortRate.HullWhite(0.1, 0.005, Yield.Constant(Continuous(0.03)))
            hw_hi = ShortRate.HullWhite(0.1, 0.020, Yield.Constant(Continuous(0.03)))
            cap_lo = present_value(hw_lo, Option.Cap(0.03, 4, 2.0))
            cap_hi = present_value(hw_hi, Option.Cap(0.03, 4, 2.0))
            @test cap_hi > cap_lo
        end
    end

    @testset "Hull-White swaptions" begin
        # Reference: google/tf-quant-finance swaption_test.py
        # a=0.03, σ=0.02, flat 1% continuous, 1y into 1y quarterly, fixed rate 1.1%
        hw_sw = ShortRate.HullWhite(0.03, 0.02, Yield.Constant(Continuous(0.01)))

        @testset "tf-quant-finance swaption reference values" begin
            # Payer swaption: expiry=1, swap maturity=2, strike=0.011, quarterly
            # tf-quant-finance expected: 0.7163243383624043 for notional=100
            payer = Option.Swaption(1.0, 2.0, 0.011, 4; payer = true)
            payer_price = present_value(hw_sw, payer)
            @test payer_price * 100 ≈ 0.71632434 atol = 0.01

            # Receiver swaption
            # tf-quant-finance expected: 0.813482544626056 for notional=100
            receiver = Option.Swaption(1.0, 2.0, 0.011, 4; payer = false)
            receiver_price = present_value(hw_sw, receiver)
            @test receiver_price * 100 ≈ 0.81348254 atol = 0.01
        end

        @testset "put-call parity for swaptions" begin
            # Payer - Receiver = PV of forward swap
            hw_pcp = ShortRate.HullWhite(0.1, 0.015, Yield.Constant(Continuous(0.05)))
            expiry = 1.0
            swap_mat = 6.0
            strike = 0.05
            freq = 2
            τ = 1.0 / freq

            payer    = present_value(hw_pcp, Option.Swaption(expiry, swap_mat, strike, freq; payer = true))
            receiver = present_value(hw_pcp, Option.Swaption(expiry, swap_mat, strike, freq; payer = false))

            # Compute the forward swap value from the curve
            n_payments = round(Int, (swap_mat - expiry) * freq)
            payment_times = [expiry + i * τ for i in 1:n_payments]
            annuity = sum(discount(hw_pcp, Ti) for Ti in payment_times)
            fwd_swap = discount(hw_pcp, expiry) - discount(hw_pcp, swap_mat) - strike * τ * annuity
            @test payer - receiver ≈ fwd_swap atol = 1e-6
        end

        @testset "swaption increases with volatility" begin
            hw_lo = ShortRate.HullWhite(0.1, 0.005, Yield.Constant(Continuous(0.03)))
            hw_hi = ShortRate.HullWhite(0.1, 0.020, Yield.Constant(Continuous(0.03)))
            sw_lo = present_value(hw_lo, Option.Swaption(1.0, 6.0, 0.03, 2))
            sw_hi = present_value(hw_hi, Option.Swaption(1.0, 6.0, 0.03, 2))
            @test sw_hi > sw_lo
        end
    end

    @testset "Hull-White fit to swaptions" begin
        # Demonstrate that fit() works for Hull-White with derivative quotes
        hw0 = ShortRate.HullWhite(0.05, 0.01, Yield.Constant(Continuous(0.03)))

        # Generate "market" swaption prices from a known model
        hw_true = ShortRate.HullWhite(0.1, 0.015, Yield.Constant(Continuous(0.03)))
        instruments = [
            Option.Swaption(1.0, 6.0, 0.03, 2),
            Option.Swaption(2.0, 7.0, 0.03, 2),
            Option.Swaption(3.0, 8.0, 0.03, 2),
        ]
        quotes = [Quote(present_value(hw_true, inst), inst) for inst in instruments]

        hw_fit = fit(hw0, quotes)
        @test hw_fit isa ShortRate.HullWhite

        # Fitted model should reprice the quotes
        for q in quotes
            @test present_value(hw_fit, q.instrument) ≈ q.price rtol = 0.05
        end
    end

    @testset "RatePath as yield model" begin
        # Construct a RatePath manually and verify it works
        import DataInterpolations
        ts = [0.0, 1.0, 2.0, 3.0]
        cum = [0.0, 0.03, 0.065, 0.10]  # cumulative integral of rates
        interp = DataInterpolations.LinearInterpolation(
            cum, ts;
            extrapolation = DataInterpolations.ExtrapolationType.Extension
        )
        rp = RatePath(interp)

        @test discount(rp, 0.0) ≈ 1.0
        @test discount(rp, 1.0) ≈ exp(-0.03)
        @test discount(rp, 2.0) ≈ exp(-0.065)
        @test zero(rp, 2.0) isa FinanceCore.Rate
    end

    @testset "Degenerate parameters" begin
        @testset "σ = 0 (deterministic)" begin
            # Vasicek with σ=0: deterministic, r(t) = b + (r0-b)exp(-at)
            # Discount: exp(-(bT + (r0-b)(1-exp(-aT))/a))
            v0 = ShortRate.Vasicek(0.5, 0.04, 0.0, Continuous(0.02))
            T = 5.0
            a, b, r0 = 0.5, 0.04, 0.02
            expected = exp(-(b * T + (r0 - b) * (1 - exp(-a * T)) / a))
            @test discount(v0, T) ≈ expected rtol = 1e-8

            # CIR with σ=0: same deterministic formula
            cir0 = ShortRate.CoxIngersollRoss(0.5, 0.04, 0.0, Continuous(0.02))
            @test discount(cir0, T) ≈ expected rtol = 1e-8
        end

        @testset "a = 0 (no mean reversion)" begin
            # Vasicek with a=0: dr = σdW, r(t) = r0 + σW(t)
            # P(0,T) = E[exp(-∫r ds)] = exp(-r0*T + σ²T³/6)  (MGF of normal)
            v_a0 = ShortRate.Vasicek(0.0, 0.04, 0.01, Continuous(0.03))
            T = 5.0
            expected = exp(-0.03 * T + 0.01^2 * T^3 / 6)
            @test discount(v_a0, T) ≈ expected rtol = 1e-8

            # Larger σ makes the convexity correction more visible
            v_a0_bigσ = ShortRate.Vasicek(0.0, 0.04, 0.10, Continuous(0.03))
            expected_bigσ = exp(-0.03 * T + 0.10^2 * T^3 / 6)
            @test discount(v_a0_bigσ, T) ≈ expected_bigσ rtol = 1e-8
            # P > 1 because convexity adjustment dominates the short rate
            @test expected_bigσ > 1.0

            # Convergence: tiny a should match a=0 limit
            v_tiny = ShortRate.Vasicek(1e-8, 0.04, 0.10, Continuous(0.03))
            @test discount(v_tiny, T) ≈ expected_bigσ rtol = 1e-4

            # CIR with a=0: formula should still work (γ = σ√2)
            cir_a0 = ShortRate.CoxIngersollRoss(0.0, 0.04, 0.05, Continuous(0.03))
            @test discount(cir_a0, 1.0) > 0
            @test discount(cir_a0, 0.0) ≈ 1.0
        end

        @testset "large a (instant mean reversion)" begin
            # With very large a, long-term rate ≈ b
            v_big = ShortRate.Vasicek(50.0, 0.05, 0.01, Continuous(0.10))
            z = zero(v_big, 10.0)
            @test abs(FinanceCore.rate(z) - 0.05) < 0.005
        end

        @testset "negative b in Vasicek" begin
            # Vasicek supports negative long-term mean
            v_neg = ShortRate.Vasicek(0.5, -0.01, 0.01, Continuous(0.02))
            @test discount(v_neg, 5.0) > 0
            # discount factor can exceed 1 when rates go negative
            z = zero(v_neg, 30.0)
            @test FinanceCore.rate(z) < 0
        end

        @testset "very short maturity" begin
            v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))
            @test discount(v, 1e-15) ≈ 1.0 atol = 1e-10

            cir = ShortRate.CoxIngersollRoss(0.3, 0.05, 0.1, Continuous(0.03))
            @test discount(cir, 1e-15) ≈ 1.0 atol = 1e-10

            hw = ShortRate.HullWhite(0.1, 0.01, Yield.Constant(0.03))
            @test discount(hw, 1e-15) ≈ 1.0 atol = 1e-10
        end
    end

    @testset "Continuous(b) constructors" begin
        # Vasicek: Continuous(b) should be equivalent to passing the float
        v1 = ShortRate.Vasicek(0.5, 0.04, 0.01, Continuous(0.02))
        v2 = ShortRate.Vasicek(0.5, Continuous(0.04), 0.01, Continuous(0.02))
        @test discount(v1, 5.0) ≈ discount(v2, 5.0)

        # CIR: same
        c1 = ShortRate.CoxIngersollRoss(0.3, 0.05, 0.1, Continuous(0.03))
        c2 = ShortRate.CoxIngersollRoss(0.3, Continuous(0.05), 0.1, Continuous(0.03))
        @test discount(c1, 5.0) ≈ discount(c2, 5.0)
    end

    @testset "ZCB option validation" begin
        hw = ShortRate.HullWhite(0.1, 0.01, Yield.Constant(Continuous(0.05)))
        # S must be > T
        @test_throws ArgumentError present_value(hw, Option.ZCBCall(5.0, 3.0, 0.8))
        @test_throws ArgumentError present_value(hw, Option.ZCBPut(5.0, 3.0, 0.8))
        # S == T should also fail
        @test_throws ArgumentError present_value(hw, Option.ZCBCall(5.0, 5.0, 0.8))
    end

    @testset "HW drift consistency" begin
        # Under risk-neutral measure, E[discount(scenario, T)] should ≈ discount(curve, T)
        curve = Yield.Constant(Continuous(0.05))
        hw = ShortRate.HullWhite(0.1, 0.02, curve)
        rng = Random.MersenneTwister(42)
        scenarios = simulate(hw; n_scenarios = 5000, timestep = 1 / 12, horizon = 11.0, rng)
        for T in [1.0, 5.0, 10.0]
            expected_disc = sum(discount(sc, T) for sc in scenarios) / 5000
            curve_disc = discount(curve, T)
            @test expected_disc ≈ curve_disc rtol = 0.02
        end
    end

    @testset "MC convergence (multi-rep RMSE)" begin
        # Use multiple repetitions to robustly test that MC error ∝ 1/√N.
        # Single-run error comparisons are unreliable because the ratio of
        # two half-normal random variables has enormous variance.
        v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))
        bond = Bond.Fixed(0.05, Periodic(2), 5)
        analytical_pv = present_value(v, bond)

        M = 30  # independent repetitions
        N_small = 500
        N_large = 2000  # 4× more scenarios → RMSE should be ~2× smaller

        errors_small = Vector{Float64}(undef, M)
        errors_large = Vector{Float64}(undef, M)
        for rep in 1:M
            rng_s = Random.MersenneTwister(1000 + rep)
            mc_s = pv_mc(v, bond; n_scenarios = N_small, timestep = 1 / 12, rng = rng_s)
            errors_small[rep] = (mc_s - analytical_pv)^2

            rng_l = Random.MersenneTwister(2000 + rep)
            mc_l = pv_mc(v, bond; n_scenarios = N_large, timestep = 1 / 12, rng = rng_l)
            errors_large[rep] = (mc_l - analytical_pv)^2
        end
        rmse_small = sqrt(sum(errors_small) / M)
        rmse_large = sqrt(sum(errors_large) / M)
        ratio = rmse_small / rmse_large

        # Expect ratio ≈ 2.0 (= √4); allow generous range for finite-sample noise
        @test 1.2 < ratio < 3.5
    end

    # ──────────────────────────────────────────────────────────────────────────
    # Literature-validated ZCB prices
    # ──────────────────────────────────────────────────────────────────────────

    @testset "Literature reference: Vasicek ZCB prices" begin
        @testset "Reuvers (2020): a=0.86, b=0.08, σ=0.01, r₀=0.06" begin
            # Source: hannoreuvers.github.io/post/vasicek/ (MC-validated)
            v = ShortRate.Vasicek(0.86, 0.08, 0.01, Continuous(0.06))
            @test discount(v, 1.0) ≈ 0.935591823311 atol = 1e-8
            @test discount(v, 5.0) ≈ 0.686027543267 atol = 1e-8
            @test discount(v, 10.0) ≈ 0.460155726152 atol = 1e-8
            @test discount(v, 30.0) ≈ 0.093029933048 atol = 1e-8

            # Asymptotic long rate: R∞ = b − σ²/(2a²)
            R_inf = 0.08 - 0.01^2 / (2 * 0.86^2)
            z_100 = FinanceCore.rate(zero(v, 100.0))
            @test z_100 ≈ R_inf atol = 0.001
        end

        @testset "R-bloggers: a=0.3, b=0.10, σ=0.03, r₀=0.03" begin
            # Source: r-bloggers.com/2010/04/fun-with-the-vasicek-interest-rate-model/
            v = ShortRate.Vasicek(0.3, 0.10, 0.03, Continuous(0.03))
            @test discount(v, 1.0) ≈ 0.961362489229 atol = 1e-8
            @test discount(v, 5.0) ≈ 0.732195597561 atol = 1e-8
            @test discount(v, 10.0) ≈ 0.471590272742 atol = 1e-8
            @test discount(v, 30.0) ≈ 0.071240674775 atol = 1e-8
        end

        @testset "dpicone1: a=3.0, b=0.02, σ=0.15, r₀=0.05" begin
            # Source: github.com/dpicone1/Vasicek_CIR_HoLee_HullWhite_Models_Python
            v = ShortRate.Vasicek(3.0, 0.02, 0.15, Continuous(0.05))
            @test discount(v, 5.0) ≈ 0.900887404280 atol = 1e-8
            @test discount(v, 10.0) ≈ 0.820267313315 atol = 1e-8
        end
    end

    @testset "CIR ZCB analytical regression" begin
        # These values are computed from the CIR closed-form ZCB formula
        # (Cox, Ingersoll & Ross, 1985, Econometrica, Eq. 23). They serve
        # as regression tests across Feller-condition regimes.

        @testset "Feller satisfied: a=0.3, b=0.08, σ=0.15, r₀=0.05" begin
            # 2ab = 0.048 > σ² = 0.0225
            cir = ShortRate.CoxIngersollRoss(0.3, 0.08, 0.15, Continuous(0.05))
            @test discount(cir, 1.0) ≈ 0.947503053857 atol = 1e-8
            @test discount(cir, 5.0) ≈ 0.731755780942 atol = 1e-8
            @test discount(cir, 10.0) ≈ 0.514277221058 atol = 1e-8
            @test discount(cir, 30.0) ≈ 0.122202611804 atol = 1e-8
        end

        @testset "Higher vol: a=0.5, b=0.10, σ=0.20, r₀=0.05" begin
            # 2ab = 0.10 > σ² = 0.04
            cir = ShortRate.CoxIngersollRoss(0.5, 0.10, 0.20, Continuous(0.05))
            @test discount(cir, 1.0) ≈ 0.941394105702 atol = 1e-8
            @test discount(cir, 5.0) ≈ 0.673617793268 atol = 1e-8
            @test discount(cir, 10.0) ≈ 0.424642886353 atol = 1e-8
            @test discount(cir, 30.0) ≈ 0.066027931612 atol = 1e-8
        end

        @testset "Feller boundary: σ = √(2ab)" begin
            # At the exact Feller boundary: 2ab = σ², process can touch zero
            a, b = 0.5, 0.10
            σ_bdy = sqrt(2 * a * b)  # ≈ 0.3162
            cir = ShortRate.CoxIngersollRoss(a, b, σ_bdy, Continuous(0.05))
            @test discount(cir, 1.0) ≈ 0.941755369937 atol = 1e-8
            @test discount(cir, 5.0) ≈ 0.685261321235 atol = 1e-8
            @test discount(cir, 10.0) ≈ 0.447811827223 atol = 1e-8
        end

        @testset "Feller violated: σ=0.50 > √(2ab)=0.316" begin
            # 2ab = 0.10 < σ² = 0.25 — formula still valid, higher convexity
            cir = ShortRate.CoxIngersollRoss(0.5, 0.10, 0.50, Continuous(0.05))
            @test discount(cir, 1.0) ≈ 0.942631527602 atol = 1e-8
            @test discount(cir, 5.0) ≈ 0.708602161384 atol = 1e-8
            @test discount(cir, 10.0) ≈ 0.491497678238 atol = 1e-8

            # Higher vol → more convexity → higher bond prices
            cir_lo = ShortRate.CoxIngersollRoss(0.5, 0.10, 0.10, Continuous(0.05))
            @test discount(cir, 5.0) > discount(cir_lo, 5.0)
        end

        @testset "QuantLib near-deterministic: a=1.0, b=0.1, σ≈0, r₀=0.1" begin
            # Source: QuantLib test-suite/shortratemodels.cpp
            # testExtendedCoxIngersollRossDiscountFactor
            cir = ShortRate.CoxIngersollRoss(1.0, 0.1, 1e-4, Continuous(0.1))
            @test discount(cir, 1.0) ≈ exp(-0.1 * 1.0) atol = 1e-6
            @test discount(cir, 5.0) ≈ exp(-0.1 * 5.0) atol = 1e-5
            @test discount(cir, 10.0) ≈ exp(-0.1 * 10.0) atol = 1e-5
        end

        @testset "CIR asymptotic long rate" begin
            # R∞ = 2ab/(a + γ) where γ = √(a² + 2σ²)
            # Source: Cox, Ingersoll & Ross (1985), Eq. 23
            a, b, σ = 0.3, 0.08, 0.15
            γ = sqrt(a^2 + 2σ^2)
            R_inf = 2a * b / (a + γ)
            cir = ShortRate.CoxIngersollRoss(a, b, σ, Continuous(0.05))
            z_long = FinanceCore.rate(zero(cir, 500.0))
            @test z_long ≈ R_inf atol = 0.0005

            # Higher σ → lower R∞ (more convexity drag)
            γ_hi = sqrt(a^2 + 2 * 0.50^2)
            R_inf_hi = 2a * b / (a + γ_hi)
            @test R_inf_hi < R_inf
        end
    end

    # ──────────────────────────────────────────────────────────────────────────
    # Vasicek yield curve shapes
    # ──────────────────────────────────────────────────────────────────────────

    @testset "Vasicek yield curve shapes" begin
        @testset "Normal (upward sloping): r₀ < b" begin
            # r₀=0.03 < b=0.10, low σ → monotonically increasing toward R∞≈0.0992
            v = ShortRate.Vasicek(0.5, 0.10, 0.02, Continuous(0.03))
            tenors = [1.0, 5.0, 10.0, 20.0, 30.0]
            rates = [FinanceCore.rate(zero(v, T)) for T in tenors]
            for i in 1:length(rates)-1
                @test rates[i+1] > rates[i]  # monotonically increasing
            end
            @test rates[1] > 0.03  # above initial short rate
        end

        @testset "Inverted (downward sloping): r₀ > b" begin
            v = ShortRate.Vasicek(0.5, 0.03, 0.02, Continuous(0.10))
            tenors = [1.0, 5.0, 10.0, 20.0, 30.0]
            rates = [FinanceCore.rate(zero(v, T)) for T in tenors]
            for i in 1:length(rates)-1
                @test rates[i+1] < rates[i]  # monotonically decreasing
            end
            @test rates[1] < 0.10  # below initial short rate
        end

        @testset "Humped: r₀ < b but R∞ < r₀ (convexity drag)" begin
            # a=0.3, b=0.10, σ=0.08, r₀=0.08
            # R∞ = b − σ²/(2a²) = 0.10 − 0.0064/0.18 ≈ 0.0644
            v = ShortRate.Vasicek(0.3, 0.10, 0.08, Continuous(0.08))
            R_inf = 0.10 - 0.08^2 / (2 * 0.3^2)
            @test R_inf < 0.08  # R∞ < r₀ confirms hump is possible

            # Short end rises (mean reversion pulls toward b > r₀)
            z_short = FinanceCore.rate(zero(v, 1.0))
            @test z_short > 0.08

            # Long end falls below r₀ (convexity drag dominates)
            z_long = FinanceCore.rate(zero(v, 30.0))
            @test z_long < 0.08

            # Peak exists: mid-tenor rate > long-tenor rate
            z_mid = FinanceCore.rate(zero(v, 2.0))
            @test z_mid > z_long
        end
    end

    # ──────────────────────────────────────────────────────────────────────────
    # Hull-White additional literature references
    # ──────────────────────────────────────────────────────────────────────────

    @testset "Hull-White additional references" begin
        @testset "Hull textbook variant: a=0.08, σ=0.02, flat 5%, K=0.70" begin
            # Source: Hull "Options, Futures, and Other Derivatives" variant
            # Also verified in github.com/Royiswho/Use-Hull-White-Model-to-Valuate-European-Call-on-Bond
            hw = ShortRate.HullWhite(0.08, 0.02, Yield.Constant(Continuous(0.05)))
            call = present_value(hw, Option.ZCBCall(1.0, 5.0, 0.70))
            # Expected: 11.3077 per 100 face (0.113077 per unit)
            @test call ≈ 0.113077 atol = 0.0005
        end
    end

    # ──────────────────────────────────────────────────────────────────────────
    # Analytical moment matching (simulation quality)
    # ──────────────────────────────────────────────────────────────────────────

    @testset "Analytical moment matching" begin
        @testset "Vasicek r(T) moments" begin
            # Vasicek r(T)|r(0) ~ Normal(E[r(T)], Var[r(T)])
            # E[r(T)] = b + (r₀−b)exp(−aT)
            # Var[r(T)] = σ²(1−exp(−2aT))/(2a)
            a, b, σ, r0 = 0.3, 0.10, 0.03, 0.03
            T = 5.0
            E_rT = b + (r0 - b) * exp(-a * T)
            V_rT = σ^2 * (1 - exp(-2a * T)) / (2a)

            N = 10_000
            dt = 1 / 52  # weekly steps to reduce EM bias
            n_steps = ceil(Int, T / dt)
            sqrt_dt = sqrt(dt)
            rng = Random.MersenneTwister(42)
            rates = Vector{Float64}(undef, N)
            for i in 1:N
                r = Float64(r0)
                for _ in 1:n_steps
                    r = r + a * (b - r) * dt + σ * sqrt_dt * randn(rng)
                end
                rates[i] = r
            end
            sample_mean = sum(rates) / N
            sample_var = sum((r - sample_mean)^2 for r in rates) / (N - 1)

            # Mean: 4σ/√N tolerance (99.994% pass rate)
            @test sample_mean ≈ E_rT atol = 4 * sqrt(V_rT) / sqrt(N)
            # Variance: 5% relative tolerance
            @test sample_var ≈ V_rT rtol = 0.05
        end

        @testset "CIR r(T) moments" begin
            # CIR r(T)|r(0) has known moments:
            # E[r(T)] = b + (r₀−b)exp(−aT)
            # Var[r(T)] = r₀σ²exp(−aT)(1−exp(−aT))/a + bσ²(1−exp(−aT))²/(2a)
            a, b, σ, r0 = 0.3, 0.08, 0.15, 0.05
            T = 5.0
            E_rT = b + (r0 - b) * exp(-a * T)
            V_rT = r0 * σ^2 * exp(-a * T) * (1 - exp(-a * T)) / a +
                   b * σ^2 * (1 - exp(-a * T))^2 / (2a)

            N = 10_000
            dt = 1 / 52
            n_steps = ceil(Int, T / dt)
            sqrt_dt = sqrt(dt)
            rng = Random.MersenneTwister(42)
            rates = Vector{Float64}(undef, N)
            for i in 1:N
                r = Float64(r0)
                for _ in 1:n_steps
                    r_pos = max(r, 0.0)
                    r = max(r_pos + a * (b - r_pos) * dt + σ * sqrt(r_pos) * sqrt_dt * randn(rng), 0.0)
                end
                rates[i] = r
            end
            sample_mean = sum(rates) / N
            sample_var = sum((r - sample_mean)^2 for r in rates) / (N - 1)

            # Wider tolerances for CIR (full truncation scheme introduces small bias)
            @test sample_mean ≈ E_rT atol = 4 * sqrt(V_rT) / sqrt(N)
            @test sample_var ≈ V_rT rtol = 0.10
        end
    end

    # ──────────────────────────────────────────────────────────────────────────
    # Martingale tests: E^Q[exp(−∫r ds)] = P(0,T)
    # ──────────────────────────────────────────────────────────────────────────

    @testset "Martingale test E[D(T)] = P(0,T)" begin
        @testset "Vasicek" begin
            # Under risk-neutral measure, mean simulated discount factor
            # must equal the analytical ZCB price
            v = ShortRate.Vasicek(0.3, 0.10, 0.03, Continuous(0.03))
            rng = Random.MersenneTwister(42)
            scenarios = simulate(v; n_scenarios = 10_000, timestep = 1 / 12, horizon = 11.0, rng)
            for T in [1.0, 5.0, 10.0]
                mc_disc = sum(discount(sc, T) for sc in scenarios) / 10_000
                @test mc_disc ≈ discount(v, T) rtol = 0.015
            end
        end

        @testset "CIR" begin
            cir = ShortRate.CoxIngersollRoss(0.3, 0.08, 0.15, Continuous(0.05))
            rng = Random.MersenneTwister(42)
            scenarios = simulate(cir; n_scenarios = 10_000, timestep = 1 / 12, horizon = 11.0, rng)
            for T in [1.0, 5.0, 10.0]
                mc_disc = sum(discount(sc, T) for sc in scenarios) / 10_000
                @test mc_disc ≈ discount(cir, T) rtol = 0.02
            end
        end
    end
end
