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

    @testset "Real return treated as continuous" begin
        shifted = base + (z, t) -> z.continuous_value + 0.01
        for t in [1.0, 5.0, 10.0]
            @test zero(shifted, t) ≈ Continuous(0.06)
            @test discount(shifted, t) ≈ exp(-0.06 * t)
        end
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

    @testset "edge cases" begin
        shifted = base + (z, t) -> z + Continuous(0.01)

        # t = 0 should give discount factor of 1
        @test discount(shifted, 0.0) == 1.0

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
end
