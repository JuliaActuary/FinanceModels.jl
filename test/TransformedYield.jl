@testset "TransformedYield" begin
    base = Yield.Constant(Continuous(0.05))

    @testset "parallel shift via Rate arithmetic" begin
        shifted = base + (z, t) -> z + Continuous(0.01)
        for t in [1.0, 5.0, 10.0, 30.0]
            @test zero(shifted, t) ≈ Continuous(0.06)
            @test discount(shifted, t) ≈ exp(-0.06 * t)
        end
    end

    @testset "Periodic rate shift" begin
        shifted = base + (z, t) -> z + Periodic(0.01, 1)
        # Periodic(0.01, 1) → Continuous(log(1.01))
        expected = 0.05 + log(1.01)
        for t in [1.0, 5.0, 10.0]
            @test zero(shifted, t) ≈ Continuous(expected)
            @test discount(shifted, t) ≈ exp(-expected * t)
        end
    end

    @testset "direct constructor" begin
        shifted = Yield.TransformedYield(base, (z, t) -> z + Continuous(0.01))
        @test zero(shifted, 5.0) ≈ Continuous(0.06)
        @test discount(shifted, 10.0) ≈ exp(-0.06 * 10.0)
    end

    @testset "commutative f + curve" begin
        f = (z, t) -> z + Continuous(0.01)
        lhs = base + f
        rhs = f + base
        for t in [1.0, 5.0, 10.0]
            @test zero(lhs, t) ≈ zero(rhs, t)
            @test discount(lhs, t) ≈ discount(rhs, t)
        end
    end

    @testset "identity transform" begin
        identity_curve = base + (z, t) -> z
        for t in [0.5, 1.0, 5.0, 20.0]
            @test discount(identity_curve, t) ≈ discount(base, t)
            @test zero(identity_curve, t) ≈ zero(base, t)
        end
    end

    @testset "tenor-dependent twist" begin
        twist = base + (z, t) -> z + Continuous(0.02 * max(0.0, 1.0 - t / 30.0))
        # At t=1, shift ≈ 0.02*(1 - 1/30); at t=15, shift = 0.01; at t≥30, shift = 0
        @test zero(twist, 1.0) ≈ Continuous(0.05 + 0.02 * (1.0 - 1.0 / 30.0))
        @test zero(twist, 15.0) ≈ Continuous(0.05 + 0.01)
        @test zero(twist, 30.0) ≈ Continuous(0.05)
        @test zero(twist, 50.0) ≈ Continuous(0.05)
    end

    @testset "Real return raises TypeError (strict Rate contract)" begin
        # Rules must return a Rate; a plain Real triggers TypeError at call time.
        shifted = base + (z, t) -> z.continuous_value + 0.01
        @test_throws TypeError zero(shifted, 5.0)
        @test_throws TypeError discount(shifted, 5.0)
    end

    @testset "composition with CompositeYield" begin
        spread_curve = Yield.Constant(Continuous(0.02))
        composite = base + spread_curve  # CompositeYield
        shifted = composite + (z, t) -> z + Continuous(0.005)
        @test zero(shifted, 5.0) ≈ Continuous(0.075)
    end

    @testset "NelsonSiegel base" begin
        ns = Yield.NelsonSiegel(1.5, 0.04, -0.02, 0.01)
        shifted = ns + (z, t) -> z + Continuous(0.005)
        for t in [1.0, 5.0, 10.0]
            z_base = zero(ns, t)
            z_shifted = zero(shifted, t)
            @test z_shifted ≈ Continuous(z_base.continuous_value + 0.005)
        end
    end

    @testset "ZeroRateCurve round-trip" begin
        zqs = ZCBYield.([0.04, 0.05, 0.06], [1.0, 5.0, 10.0])
        zrc = fit(Spline.Linear(), zqs, Fit.Bootstrap())
        shifted = zrc + (z, t) -> z + Continuous(0.01)
        for t in [1.0, 5.0, 10.0]
            z_base = zero(zrc, t)
            @test zero(shifted, t) ≈ Continuous(z_base.continuous_value + 0.01)
        end
    end

    @testset "negative base rate" begin
        neg = Yield.Constant(Continuous(-0.01))
        shifted = neg + (z, t) -> z + Continuous(0.02)
        @test zero(shifted, 5.0) ≈ Continuous(0.01)
        @test discount(shifted, 5.0) ≈ exp(-0.01 * 5.0)
    end

    @testset "edge cases" begin
        shifted = base + (z, t) -> z + Continuous(0.01)

        # t = 0 should give discount factor of 1
        @test discount(shifted, 0.0) == 1.0

        # zero(curve, 0.0) is NaN for most base curves (generic fallback: -log(1)/0 = 0/0)
        # This is pre-existing behavior, not specific to TransformedYield
        @test isnan(FinanceCore.rate(zero(shifted, 0.0)))

        # very small tenor
        @test isfinite(discount(shifted, 1e-10))

        # large tenor
        @test discount(shifted, 100.0) ≈ exp(-0.06 * 100.0)
    end

    @testset "forward rates" begin
        shifted = base + (z, t) -> z + Continuous(0.01)
        # For a flat base + flat shift, forward rate should equal the shifted zero rate
        fwd = forward(shifted, 5.0)
        @test fwd ≈ Continuous(0.06)
    end

    @testset "no dispatch conflict" begin
        # base + scalar → not TransformedYield (CompositeYield or Constant)
        @test !((base + 0.01) isa Yield.TransformedYield)
        # base + Rate → not TransformedYield
        @test !((base + Continuous(0.01)) isa Yield.TransformedYield)
        # base + Function → TransformedYield
        @test (base + (z, t) -> z) isa Yield.TransformedYield
    end

    @testset "TransformedYield alias" begin
        # TransformedYield is now a deprecation alias for TenorShift.
        @test Yield.TransformedYield === Yield.TenorShift
        ts = Yield.TenorShift(base, (z, t) -> z + Continuous(0.01))
        @test ts isa Yield.TransformedYield
        @test ts isa Yield.AbstractYieldShift
    end
end

@testset "ProjectedShift" begin
    base = Yield.Constant(Continuous(0.05))

    # A reusable year-curried rule: -150 bp parallel, phased in linearly over 10y.
    phase_in = (τ, z, _) -> z + Continuous(-0.015 * min(τ, 10) / 10)

    @testset "basic phase-in" begin
        # At τ=0 the shift is 0; curve matches base.
        c0 = Yield.ProjectedShift(base, phase_in, 0.0)
        for t in [1.0, 5.0, 10.0]
            @test zero(c0, t) ≈ zero(base, t)
            @test discount(c0, t) ≈ discount(base, t)
        end

        # At τ=3 the shift is -45 bp at every tenor.
        c3 = Yield.ProjectedShift(base, phase_in, 3.0)
        for t in [1.0, 5.0, 10.0, 30.0]
            @test zero(c3, t) ≈ Continuous(0.05 - 0.0045)
        end

        # At τ=10 (and beyond) the shift is fully phased in at -150 bp.
        for τ in [10.0, 15.0, 50.0]
            c = Yield.ProjectedShift(base, phase_in, τ)
            @test zero(c, 5.0) ≈ Continuous(0.05 - 0.015)
        end
    end

    @testset "equivalence with hand-rolled TenorShift closure" begin
        # For every (τ, t) grid point, ProjectedShift(base, rule, τ) must
        # match TenorShift(base, (z, u) -> rule(τ, z, u)).
        for τ in [0.0, 1.0, 5.0, 10.0, 25.0], t in [0.5, 1.0, 5.0, 10.0, 30.0]
            ps = Yield.ProjectedShift(base, phase_in, τ)
            ts = Yield.TenorShift(base, (z, u) -> phase_in(τ, z, u))
            @test zero(ps, t) ≈ zero(ts, t)
            @test discount(ps, t) ≈ discount(ts, t)
        end
    end

    @testset "both axes (τ and t) dependency" begin
        # Steepener that decays in tenor and fades across projection time.
        steepener_fade = (τ, z, t) -> z + Continuous(0.02 * max(0.0, 1.0 - t/30.0) * exp(-τ/10))

        # τ=0: full strength, tenor-decaying twist.
        c0 = Yield.ProjectedShift(base, steepener_fade, 0.0)
        @test zero(c0, 1.0) ≈ Continuous(0.05 + 0.02 * (1.0 - 1.0/30.0))
        @test zero(c0, 30.0) ≈ Continuous(0.05)  # twist dies at 30y

        # τ=10: shift is multiplied by exp(-1).
        c10 = Yield.ProjectedShift(base, steepener_fade, 10.0)
        @test zero(c10, 1.0) ≈ Continuous(0.05 + 0.02 * (1.0 - 1.0/30.0) * exp(-1.0))
    end

    @testset "Real return raises TypeError (strict Rate contract)" begin
        # Rules must return a Rate; a plain Real triggers TypeError at call time.
        rule_real = (τ, z, _) -> z.continuous_value + 0.01 * τ
        c = Yield.ProjectedShift(base, rule_real, 2.0)
        @test_throws TypeError zero(c, 5.0)
        @test_throws TypeError discount(c, 5.0)
    end

    @testset "Periodic-shaped Rate return is properly converted (no silent miscoercion)" begin
        # A rule returning a Rate in Periodic convention must be converted to
        # continuous, not silently misinterpreted (the footgun the strict
        # contract prevents).
        rule_periodic = (τ, z, _) -> z + Periodic(0.01 * τ, 1)  # +1% × τ in annual-effective
        c = Yield.ProjectedShift(base, rule_periodic, 2.0)
        # Expected: base 5% continuous, plus Periodic(0.02, 1) converted = log(1.02)
        expected = Continuous(0.05 + log(1.02))
        @test zero(c, 5.0) ≈ expected
        @test discount(c, 5.0) ≈ exp(-expected.continuous_value * 5.0)
    end

    @testset "edge cases" begin
        ps = Yield.ProjectedShift(base, phase_in, 5.0)

        # t = 0 → discount factor of 1, per AbstractYieldShift contract.
        @test discount(ps, 0.0) == 1.0

        # very small and very large tenor
        @test isfinite(discount(ps, 1e-10))
        @test discount(ps, 100.0) ≈ exp(-(0.05 - 0.0075) * 100.0)

        # negative τ is just numerically valid (min(τ, 10) clamps to τ when τ<10).
        # Caller-side semantics: passing a negative τ is the caller's choice.
        c_neg = Yield.ProjectedShift(base, phase_in, -2.0)
        @test zero(c_neg, 5.0) ≈ Continuous(0.05 - 0.015 * (-0.2))
    end

    @testset "composition with fitted curve" begin
        zqs = ZCBYield.([0.04, 0.05, 0.06], [1.0, 5.0, 10.0])
        zrc = fit(Spline.Linear(), zqs, Fit.Bootstrap())
        ps = Yield.ProjectedShift(zrc, phase_in, 5.0)  # -75 bp at every tenor

        for t in [1.0, 5.0, 10.0]
            z_base = zero(zrc, t)
            @test zero(ps, t) ≈ Continuous(z_base.continuous_value - 0.0075)
        end
    end

    @testset "composition with CompositeYield and ScaledYield" begin
        spread = Yield.Constant(Continuous(0.02))
        composite = base + spread  # CompositeYield: 7% combined
        ps = Yield.ProjectedShift(composite, phase_in, 3.0)
        @test zero(ps, 5.0) ≈ Continuous(0.07 - 0.0045)

        # Wrap a ProjectedShift in a ScaledYield (after-tax).
        ps_scaled = Yield.ProjectedShift(base, phase_in, 3.0) * 0.79
        @test discount(ps_scaled, 5.0) ≈ exp(-(0.05 - 0.0045) * 0.79 * 5.0)
    end

    @testset "TenorShift composes onto ProjectedShift" begin
        # ProjectedShift + (z, t) -> ... should produce a TenorShift wrapping the ProjectedShift.
        ps = Yield.ProjectedShift(base, phase_in, 5.0)  # -75 bp
        outer = ps + (z, t) -> z + Continuous(0.003)     # +30 bp on top
        @test outer isa Yield.TenorShift
        @test zero(outer, 5.0) ≈ Continuous(0.05 - 0.0075 + 0.003)
    end

    @testset "ProjectedShift composes onto ProjectedShift" begin
        # Two stackable rules: a τ-only phase-in and a τ-and-tenor twist.
        rule1 = (τ, z, _) -> z + Continuous(-0.015 * min(τ, 10) / 10)
        rule2 = (τ, z, t) -> z + Continuous(0.002 * t * τ / 10)

        @testset "same τ on both layers" begin
            τ = 5.0
            inner = Yield.ProjectedShift(base, rule1, τ)
            outer = Yield.ProjectedShift(inner, rule2, τ)
            @test outer isa Yield.ProjectedShift
            for t in [1.0, 5.0, 10.0]
                z_base = zero(base, t).continuous_value
                expected = z_base + (-0.015 * 5/10) + (0.002 * t * 5/10)
                @test zero(outer, t) ≈ Continuous(expected)
            end
        end

        @testset "different τ on each layer (each shift uses its own time field)" begin
            inner = Yield.ProjectedShift(base, rule1, 3.0)    # inner sees τ=3
            outer = Yield.ProjectedShift(inner, rule2, 10.0)  # outer sees τ=10
            for t in [1.0, 5.0, 10.0]
                z_base = zero(base, t).continuous_value
                expected = z_base + (-0.015 * 3/10) + (0.002 * t * 10/10)
                @test zero(outer, t) ≈ Continuous(expected)
            end
        end
    end

    @testset "forward rates" begin
        # For a flat base + flat-at-τ shift, forward rate should equal shifted zero rate.
        ps = Yield.ProjectedShift(base, phase_in, 10.0)  # -150 bp, fully phased
        @test forward(ps, 5.0) ≈ Continuous(0.05 - 0.015)
    end

    @testset "AbstractYieldShift supertype" begin
        ps = Yield.ProjectedShift(base, phase_in, 3.0)
        ts = Yield.TenorShift(base, (z, t) -> z + Continuous(0.01))
        @test ps isa Yield.AbstractYieldShift
        @test ts isa Yield.AbstractYieldShift
        @test Yield.AbstractYieldShift <: Yield.AbstractYieldModel
    end

    @testset "type inference" begin
        ps = Yield.ProjectedShift(base, phase_in, 3.0)
        @test (@inferred zero(ps, 5.0)) isa Rate{Float64, Continuous}
        @test (@inferred discount(ps, 5.0)) isa Float64
        @test (@inferred Yield.ProjectedShift(base, phase_in, 3.0)) isa Yield.ProjectedShift
    end

    @testset "ForwardDiff propagation through ProjectedShift" begin
        # AD must propagate through (a) the rule's parameters and (b) the .time
        # field itself. The latter is the novel capability this type enables:
        # projection-time Greeks of a stress.
        using ForwardDiff
        h = 1e-6  # central-difference step for cross-check

        # (a) ∂(discount)/∂(rule parameter): scenario-parameter sensitivity.
        f_slope(slope) = let rule = (τ, z, _) -> z + Continuous(slope * min(τ, 10) / 10)
            discount(Yield.ProjectedShift(base, rule, 5.0), 10.0)
        end
        ad_slope = ForwardDiff.derivative(f_slope, -0.015)
        fd_slope = (f_slope(-0.015 + h) - f_slope(-0.015 - h)) / (2h)
        @test isfinite(ad_slope)
        @test ad_slope ≈ fd_slope rtol = 1e-5

        # (b) ∂(discount)/∂τ: projection-time sensitivity — the new axis.
        # Differentiate at τ=5.0 (well inside the phase-in region; min is smooth here).
        f_tau(τ) = discount(Yield.ProjectedShift(base, phase_in, τ), 10.0)
        ad_tau = ForwardDiff.derivative(f_tau, 5.0)
        fd_tau = (f_tau(5.0 + h) - f_tau(5.0 - h)) / (2h)
        @test isfinite(ad_tau)
        @test ad_tau ≈ fd_tau rtol = 1e-5
    end
end
