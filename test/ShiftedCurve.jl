@testset "ShiftedCurve" begin
    base = Yield.Constant(Continuous(0.05))

    @testset "parallel shift" begin
        shifted = Yield.ShiftedCurve(base, t -> 0.01)
        for t in [1.0, 5.0, 10.0, 30.0]
            @test zero(shifted, t) ≈ Continuous(0.06)
            @test discount(shifted, t) ≈ exp(-0.06 * t)
        end
    end

    @testset "scalar convenience constructor" begin
        shifted = Yield.ShiftedCurve(base, 0.01)
        @test zero(shifted, 5.0) ≈ Continuous(0.06)
        @test discount(shifted, 10.0) ≈ exp(-0.06 * 10.0)
    end

    @testset "zero shift is identity" begin
        shifted = Yield.ShiftedCurve(base, 0.0)
        for t in [0.5, 1.0, 5.0, 20.0]
            @test discount(shifted, t) ≈ discount(base, t)
            @test zero(shifted, t) ≈ zero(base, t)
        end
    end

    @testset "tenor-dependent twist" begin
        twist = Yield.ShiftedCurve(base, t -> 0.02 * max(0.0, 1.0 - t / 30.0))
        # At t=0+, shift ≈ 0.02; at t=15, shift = 0.01; at t≥30, shift = 0
        @test zero(twist, 1.0) ≈ Continuous(0.05 + 0.02 * (1.0 - 1.0 / 30.0))
        @test zero(twist, 15.0) ≈ Continuous(0.05 + 0.01)
        @test zero(twist, 30.0) ≈ Continuous(0.05)
        @test zero(twist, 50.0) ≈ Continuous(0.05)
    end

    @testset "composition with CompositeYield" begin
        spread_curve = Yield.Constant(Continuous(0.02))
        composite = base + spread_curve  # CompositeYield
        shifted = Yield.ShiftedCurve(composite, 0.005)
        @test zero(shifted, 5.0) ≈ Continuous(0.075)
    end

    @testset "NelsonSiegel base" begin
        ns = Yield.NelsonSiegel(1.5, 0.04, -0.02, 0.01)
        shifted = Yield.ShiftedCurve(ns, 0.005)
        for t in [1.0, 5.0, 10.0]
            z_base = zero(ns, t)
            z_shifted = zero(shifted, t)
            @test z_shifted ≈ Continuous(z_base.continuous_value + 0.005)
        end
    end

    @testset "ZeroRateCurve round-trip" begin
        zqs = ZCBYield.([0.04, 0.05, 0.06], [1.0, 5.0, 10.0])
        zrc = fit(Spline.Linear(), zqs, Fit.Bootstrap())
        shifted = Yield.ShiftedCurve(zrc, 0.01)
        for t in [1.0, 5.0, 10.0]
            z_base = zero(zrc, t)
            @test zero(shifted, t) ≈ Continuous(z_base.continuous_value + 0.01)
        end
    end

    @testset "edge cases" begin
        shifted = Yield.ShiftedCurve(base, 0.01)

        # t = 0 should give discount factor of 1
        @test discount(shifted, 0.0) == 1.0

        # very small tenor
        @test isfinite(discount(shifted, 1e-10))

        # large tenor
        @test discount(shifted, 100.0) ≈ exp(-0.06 * 100.0)
    end

    @testset "forward rates" begin
        shifted = Yield.ShiftedCurve(base, t -> 0.01)
        # For a flat base + flat shift, forward rate should equal the shifted zero rate
        fwd = forward(shifted, 5.0)
        @test fwd ≈ Continuous(0.06)
    end
end
