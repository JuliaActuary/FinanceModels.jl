using ForwardDiff

@testset "Constant-forward interpolation and bootstrap" begin
    c = ZeroRateCurve([0.02, 0.04, 0.03], [1, 4, 10], Spline.ConstantForward())
    @test c isa Yield.ConstantForward
    @test discount(c, 0) == 1
    @test rate(zero(c, 0)) == 0.02
    @test discount(c, 2) ≈ exp(-(0.02 + (0.16 - 0.02) / 3))
    @test rate(forward(c, 1.1, 3.9)) ≈ (0.16 - 0.02) / 3
    @test rate(Yield.instantaneous_forward(c, 4)) ≈ (0.3 - 0.16) / 6
    @test ForwardDiff.derivative(t -> rate(zero(c, t)), 1.0) ≈ (0.16 - 0.02) / 3 - 0.02
    @test reconstruct(c) == c
    @test (@set c.rates = [0.03, 0.04, 0.03]) == reconstruct(c; rates = [0.03, 0.04, 0.03])
    @test_throws ArgumentError ZeroRateCurve([0.01, 0.02], [0, 1], Spline.ConstantForward())
    @test_throws ArgumentError reconstruct(c; extrapolation = :extension)
    @test_throws DomainError discount(c, -1)
    for policy in (:flat_forward, :flat_zero, :linear, Yield.FlatForwardAt(Continuous(0.025)))
        tail = reconstruct(c; extrapolation = policy)
        @test discount(tail, 10) == discount(c, 10)
        @test isfinite(discount(tail, 20))
    end
    @test rate(zero(Yield.ConstantForward([0.02], [3]), 100)) ≈ 0.02
    for frequency in (1, 2), rates in ([0.02, 0.035, 0.04], [-0.01, -0.005, 0.01], [0.0, 0.0, 0.0])
        qs = ParSwapYield.(rates, [2, 5, 20]; frequency)
        full = fit(Spline.ConstantForward(), qs, Fit.Bootstrap())
        for j in eachindex(qs)
            prefix = fit(Spline.ConstantForward(), qs[1:j], Fit.Bootstrap())
            for q in qs[1:j]
                @test pv(full, q.instrument) ≈ q.price atol = 1.0e-12
                @test pv(prefix, q.instrument) ≈ q.price atol = 1.0e-12
            end
        end
    end
    valuation(r) = discount(fit(Spline.ConstantForward(), ParSwapYield.(r, [1, 5, 20]; frequency = 2), Fit.Bootstrap()), 12)
    r = [0.02, 0.03, 0.04]
    fd = [(valuation(r + 1.0e-6 * (1:3 .== k)) - valuation(r - 1.0e-6 * (1:3 .== k))) / 2.0e-6 for k in 1:3]
    @test ForwardDiff.gradient(valuation, r) ≈ fd rtol = 1.0e-7
end

@testset "Exponential-forward extrapolation" begin
    base = Yield.ConstantForward([0.02, 0.03], [1, 20])
    make(u, l, a) = Yield.ExponentialForward(base; start = 10, initial_forward = Continuous(l), ufr = Continuous(u), α = a)
    c = make(0.035, 0.04, 0.11)
    @test discount(c, 0) == 1
    @test discount(c, 10) == discount(base, 10)
    @test rate(Yield.instantaneous_forward(c, 10)) == 0.04
    for seed in (-1.0, 0.0, 1.0)
        t = ForwardDiff.Dual{Nothing}(10.0, seed)
        @test ForwardDiff.value(rate(Yield.instantaneous_forward(c, t))) == 0.04
        @test ForwardDiff.partials(rate(Yield.instantaneous_forward(c, t)))[1] ≈ -0.11 * 0.005 * seed
    end
    @test ForwardDiff.derivative(t -> discount(c, 20 - t), 10.0) ≈ 0.04 * discount(c, 10)
    @test ForwardDiff.derivative(t -> rate(zero(c, t)), 0.0) == 0
    changed = @set c.α = 0.12
    @test changed.α == 0.12
    @test discount(changed, 60) ≈ discount(make(0.035, 0.04, 0.12), 60)
    moved = @set c.start = 5.0
    @test discount(moved, 5) == discount(base, 5)
    @test discount(moved, 60) ≈ discount(base, 5) * discount(moved, 5, 60)
    @test discount((@set c.base = Yield.Constant(Continuous(0.02))), 10) ≈ exp(-0.2)
    @test_throws ArgumentError (@set c.α = 0.0)
    @test_throws MethodError (@set c._log_start = 0.0)
    properties = Accessors.ConstructionBase.getproperties(c)
    rebuilt = Accessors.ConstructionBase.setproperties(c, properties)
    @test discount(rebuilt, 60) == discount(c, 60)
    @test rate(Yield.instantaneous_forward(c, 20)) ≈ 0.035 + 0.005 * exp(-1.1)
    @test discount(c, 10, 20) ≈ exp(-0.035 * 10 - 0.005 * (-expm1(-1.1)) / 0.11)
    @test discount(c, 0, 30) ≈ discount(c, 30)
    @test discount(c, 5, 30) * discount(c, 30, 5) ≈ 1
    @test discount(c, 10 + 1.0e-10) ≈ discount(base, 10) rtol = 1.0e-10
    @test rate(zero(c, Inf)) == 0.035
    @test discount(c, Inf) == 0
    @test rate(Yield.instantaneous_forward(c, 1000)) ≈ 0.035
    @test discount(make(-0.01, 0.02, 0.11), Inf) == Inf
    c0 = make(0.0, 0.02, 0.11)
    @test discount(c0, Inf) ≈ discount(base, 10) * exp(-0.02 / 0.11)
    @test discount(c0, 0, Inf) ≈ discount(c0, Inf)
    @test discount(c, Inf, Inf) == 1
    @test discount(make(0.03, 0.04, 1.0e-14), 10, 11) ≈ exp(-0.04)
    @test rate(forward(make(0.03, 0.04, nextfloat(0.0)), 10.0, 10.1)) ≈ 0.04
    setprecision(256) do
        high = Yield.ExponentialForward(base; start = 10, initial_forward = Continuous(big"0.04"), ufr = Continuous(big"0.035"), α = big"1e-6")
        expected = big"0.035" + big"0.005" * (-expm1(-big"1e-6")) / big"1e-6"
        @test -log(discount(high, 10, 11)) ≈ expected rtol = big"1e-70"
    end
    @test discount(c + Yield.Constant(Continuous(-0.035)), Inf) ≈ exp(-Yield.__log_tail(c).a0)
    @test isfinite(discount(c, 1.0e6, 1.0e6 + 1))
    g(x) = discount(make(x[1], x[2], x[3]), 60)
    x = [0.035, 0.04, 0.11]
    fd = [(g(x + 1.0e-6 * (1:3 .== k)) - g(x - 1.0e-6 * (1:3 .== k))) / 2.0e-6 for k in 1:3]
    @test ForwardDiff.gradient(g, x) ≈ fd rtol = 1.0e-7
    shifted_start(s) = discount((@set c.start = s), 60)
    @test ForwardDiff.derivative(shifted_start, 10.0) ≈ (shifted_start(10.000001) - shifted_start(9.999999)) / 0.000002 rtol = 1.0e-7
    @test_throws ArgumentError make(0.03, 0.04, 0)
    @test_throws DomainError discount(c, -1)
end
