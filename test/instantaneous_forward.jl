@testset "instantaneous_forward" begin
    # central finite difference of −log P(t), the definition of the instantaneous forward
    fd(m, t; h = 1.0e-6) = (log(discount(m, t - h)) - log(discount(m, t + h))) / (2h)

    @testset "Constant" begin
        c = Yield.Constant(Periodic(0.04, 2))
        f = Yield.instantaneous_forward(c, 5.0)
        @test f ≈ rate(convert(Continuous(), Periodic(0.04, 2)))
        @test f ≈ fd(c, 5.0) atol = 1.0e-8
        # flat curve: same forward at every tenor
        @test Yield.instantaneous_forward(c, 0.0) == Yield.instantaneous_forward(c, 30.0)
    end

    @testset "NelsonSiegel and NelsonSiegelSvensson" begin
        ns = Yield.NelsonSiegel(3.0, 0.04, -0.02, 0.03)
        for t in (0.5, 2.0, 7.3, 25.0)
            @test Yield.instantaneous_forward(ns, t) ≈ fd(ns, t) atol = 1.0e-6
        end
        # f(0) = β₀ + β₁, matching zero's t=0 limit
        @test Yield.instantaneous_forward(ns, 0.0) ≈ 0.04 - 0.02
        @test Yield.instantaneous_forward(ns, 0.0) ≈ FinanceCore.rate(zero(ns, 0.0))

        nss = Yield.NelsonSiegelSvensson(1.5, 3.0, 0.04, -0.02, 0.03, 0.015)
        for t in (0.5, 2.0, 7.3, 25.0)
            @test Yield.instantaneous_forward(nss, t) ≈ fd(nss, t) atol = 1.0e-6
        end
        @test Yield.instantaneous_forward(nss, 0.0) ≈ 0.04 - 0.02
    end

    @testset "Spline-backed ZeroRateCurve" begin
        zs = [0.02, 0.025, 0.03, 0.035]
        tenors = [1.0, 2.0, 5.0, 10.0]
        for sp in (Spline.PCHIP(), Spline.Cubic(), Spline.Linear())
            zrc = ZeroRateCurve(zs, tenors, sp)
            # off-knot points (Linear has kinks at knots where f is undefined)
            for t in (1.5, 3.7, 8.0)
                @test Yield.instantaneous_forward(zrc, t) ≈ fd(zrc, t) atol = 1.0e-5
            end
        end

        # MonotoneConvex-backed (default) delegates to the Hagan–West forward
        zrc_mc = ZeroRateCurve(zs, tenors)
        @test Yield.instantaneous_forward(zrc_mc, 3.0) ≈ fd(zrc_mc, 3.0) atol = 1.0e-5

        @test_throws DomainError Yield.instantaneous_forward(zrc_mc, -1.0)
    end

    @testset "wrapper models" begin
        a = Yield.Constant(0.03)
        b = Yield.Constant(0.01)
        ns = Yield.NelsonSiegel(3.0, 0.04, -0.02, 0.03)

        # CompositeYield: forwards add/subtract for the +/− compositions
        @test Yield.instantaneous_forward(a + b, 2.0) ≈
            Yield.instantaneous_forward(a, 2.0) + Yield.instantaneous_forward(b, 2.0)
        @test Yield.instantaneous_forward(ns - b, 2.0) ≈ fd(ns - b, 2.0) atol = 1.0e-6

        # ScaledYield: forward scales with the curve
        m = Yield.Constant(Continuous(0.05)) * 0.79
        @test Yield.instantaneous_forward(m, 2.0) ≈ 0.05 * 0.79
        @test Yield.instantaneous_forward(ns * 2.0, 4.0) ≈ fd(ns * 2.0, 4.0) atol = 1.0e-6

        # ForwardStarting: forward is the base curve's forward shifted by the start
        fs = Yield.ForwardStarting(ns, 1.0)
        @test Yield.instantaneous_forward(fs, 2.0) ≈ Yield.instantaneous_forward(ns, 3.0)
        @test Yield.instantaneous_forward(fs, 2.0) ≈ fd(fs, 2.0) atol = 1.0e-6
    end

    @testset "discount-native models error rather than approximate" begin
        sw = Yield.SmithWilson([1.0, 2.0], [0.1, 0.2]; ufr = 0.03, α = 0.1)
        @test_throws MethodError Yield.instantaneous_forward(sw, 1.0)
    end
end
