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
end
