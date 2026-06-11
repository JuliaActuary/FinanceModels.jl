using ActuaryUtilities

@testset "ActuaryUtilities.jl integration tests" begin
    cfs = [5, 5, 105]
    times = [1, 2, 3]

    # scalar rates and Rate objects: ActuaryUtilities applies a periodic-rate
    # shock, so Modified duration = Macaulay / (1 + y)
    for d in [0.03, Periodic(0.03, 1)]
        @test present_value(d, cfs, times) ≈ 105.65722270978935
        @test duration(Macaulay(), d, cfs, times) ≈ 2.86350467067113
        @test duration(d, cfs, times) ≈ 2.7801016220108057
        @test convexity(d, cfs, times) ≈ 10.625805482685939
    end

    # yield *models*: the shock is applied to the continuous zero rate
    # (curve + Δ composes in continuous-zero space), so Modified ≡ Macaulay
    @testset "Yield.Constant uses the curve (continuous-shock) convention" begin
        d = Yield.Constant(0.03)
        V = present_value(d, cfs, times)
        @test V ≈ 105.65722270978935
        mac = duration(Macaulay(), d, cfs, times)
        @test mac ≈ 2.86350467067113
        @test duration(d, cfs, times) ≈ mac
        # ∂²V/∂Δ² with d(t) = (1+r)^(-t)·(1+Δ)^(-t) gives Σ cf·d·t·(t+1) / V
        expected_cvx = sum(cf * 1.03^-t * t * (t + 1) for (cf, t) in zip(cfs, times)) / V
        @test convexity(d, cfs, times) ≈ expected_cvx
    end
end
