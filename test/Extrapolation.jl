using ForwardDiff

@testset "shared curve tails" begin
    rates = [0.02, 0.025, 0.03, 0.035, 0.04, 0.042]
    times = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0]
    tn, zn, horizon = last(times), last(rates), 100.0
    fixed = Yield.FlatForwardAt(0.035)
    policies = (:flat_forward, :flat_zero, :linear, fixed)
    descriptors = (
        Spline.Linear(), Spline.Quadratic(), Spline.Cubic(),
        Spline.PCHIP(), Spline.Akima(), Spline.BSpline(3), Spline.MonotoneConvex(),
    )

    @testset "boundary and AD: $descriptor" for descriptor in descriptors
        base = Yield.build_model(descriptor, times, rates)
        # A second-order one-sided difference measures the INTERIOR zero slope.
        h = 1.0e-3
        left_slope = (
            3rate(zero(base, tn)) - 4rate(zero(base, tn - h)) +
                rate(zero(base, tn - 2h))
        ) / (2h)
        for policy in policies
            curve = Yield.build_model(descriptor, times, rates; extrapolation = policy)
            wrapped = ZeroRateCurve(rates, times, descriptor; extrapolation = policy)
            for t in (0.0, 0.5, 3.0, 17.0, tn)
                @test zero(curve, t) == zero(base, t)
                @test zero(wrapped, t) == zero(curve, t)
            end
            @test discount(curve, tn) == discount(base, tn)
            @test zero(wrapped, horizon) == zero(curve, horizon)
            if policy === :linear
                tail_slope = (rate(zero(curve, horizon)) - zn) / (horizon - tn)
                @test tail_slope ≈ left_slope atol = 1.0e-9
            elseif policy === :flat_zero
                @test rate(zero(curve, horizon)) == zn
            else
                f = policy isa Yield.FlatForwardAt ? policy.forward : zn + tn * left_slope
                @test rate(forward(curve, tn, horizon)) ≈ f atol = 3.0e-8
                @test rate(zero(curve, 1.0e5)) ≈ f + (zn - f) * tn / 1.0e5 atol = 3.0e-8
            end

            z(rs) = rate(zero(ZeroRateCurve(rs, times, descriptor; extrapolation = policy), horizon))
            ad = ForwardDiff.gradient(z, rates)
            bump = 1.0e-7
            fd = map(eachindex(rates)) do i
                up, down = copy(rates), copy(rates)
                up[i] += bump
                down[i] -= bump
                (z(up) - z(down)) / (2bump)
            end
            @test ad ≈ fd atol = 1.0e-7
            if policy === fixed
                # An independent forward assumption removes the endpoint-derivative
                # sensitivity; only the last-knot discount factor still varies.
                @test ad ≈ [zeros(length(rates) - 1); tn / horizon] atol = 1.0e-14
            end
            zdot = ForwardDiff.derivative(t -> rate(zero(curve, t)), horizon)
            @test zdot ≈ (
                rate(zero(curve, horizon + h)) -
                    rate(zero(curve, horizon - h))
            ) / (2h) atol = 1.0e-10
            if curve isa Yield.MonotoneConvex
                @test Yield.instantaneous_forward(curve, horizon) ≈
                    rate(zero(curve, horizon)) + horizon * zdot atol = 1.0e-14
                @test Yield.instantaneous_forward(curve, tn) == last(curve.f)
            end
        end
        # Different fixed assumptions change only the tail and remain differentiable.
        zf(f) = rate(
            zero(
                ZeroRateCurve(
                    rates, times, descriptor;
                    extrapolation = Yield.FlatForwardAt(f)
                ), horizon
            )
        )
        @test ForwardDiff.derivative(zf, fixed.forward) ≈ 1 - tn / horizon
    end

    @testset "policy reconstruction" begin
        CB = Accessors.ConstructionBase
        for curve in (Yield.MonotoneConvex(rates, times), ZeroRateCurve(rates, times, Spline.Cubic()))
            for policy in policies
                c = @set curve.extrapolation = policy
                @test c.extrapolation == policy
                for rebuilt in (
                        CB.setproperties(c, CB.getproperties(c)),
                        CB.constructorof(typeof(c))(CB.getfields(c)...),
                        Accessors.mapproperties(identity, c),
                    )
                    @test rebuilt.extrapolation == policy
                    @test zero(rebuilt, horizon) == zero(c, horizon)
                end
                bumped = @set c.rates[end - 1] = c.rates[end - 1] + 0.001
                @test bumped.extrapolation == policy
                if policy === fixed
                    @test rate(forward(bumped, tn, horizon)) ≈ fixed.forward
                elseif policy === :flat_forward
                    @test rate(forward(bumped, tn, horizon)) != rate(forward(c, tn, horizon))
                end
                if c isa Yield.MonotoneConvex
                    moved = @set c.times = times .* 2
                    @test moved.extrapolation == policy
                    @test_throws ArgumentError c._tail
                    @test_throws ArgumentError (@set c._tail = nothing)
                else
                    moved = @set c.tenors = times .* 2
                    @test moved.extrapolation == policy
                    @test c == CB.setproperties(c, CB.getproperties(c))
                end
                @test_throws ArgumentError (@set c.extrapolation = :unknown)
            end
        end
        sampled = ZeroRateCurve(
            Yield.Constant(Continuous(0.03)), times;
            extrapolation = fixed
        )
        @test sampled.extrapolation == fixed
        @test rate(forward(sampled, tn, horizon)) ≈ fixed.forward
        a = ZeroRateCurve(rates, times; extrapolation = Yield.FlatForwardAt(big"0.035"))
        b = ZeroRateCurve(rates, times; extrapolation = Yield.FlatForwardAt(big"0.035"))
        @test a == b && isequal(a, b) && hash(a) == hash(b)
        @test a != (@set a.extrapolation = Yield.FlatForwardAt(0.04))
        @test Yield.FlatForwardAt(-0.0) == Yield.FlatForwardAt(0.0)
        @test !isequal(Yield.FlatForwardAt(-0.0), Yield.FlatForwardAt(0.0))
    end

    @testset "fits retain the policy and native MonotoneConvex interface" begin
        qs = ZCBYield.(Continuous.(rates), times)
        for method in (Fit.Loss(abs2), Fit.Bootstrap())
            baseline = fit(Spline.MonotoneConvex(), qs, method)
            for policy in policies
                c = fit(Spline.MonotoneConvex(), qs, method; extrapolation = policy)
                @test c isa Yield.MonotoneConvex
                @test c.rates == baseline.rates && c.times == baseline.times
                @test c.extrapolation == policy
                @test isfinite(Yield.instantaneous_forward(c, horizon))
                @test zero(c, horizon) == zero(
                    Yield.MonotoneConvex(
                        c.rates, c.times;
                        extrapolation = policy
                    ), horizon
                )
            end
        end
        for descriptor in (Spline.Linear(), Spline.MonotoneConvex(), Yield.MonotoneConvex())
            c = fit(descriptor, qs; extrapolation = fixed)
            @test rate(forward(c, tn, horizon)) ≈ fixed.forward
        end
        linear = fit(Spline.Linear(), qs, Fit.Bootstrap(); extrapolation = fixed)
        @test rate(forward(linear, tn, horizon)) ≈ fixed.forward
        for policy in (:linear, fixed)
            for c in (
                    ZeroRateCurve(rates .+ 0.001, times; extrapolation = policy),
                    Yield.MonotoneConvex(rates .+ 0.001, times; extrapolation = policy),
                )
                fitted = fit(c, qs)
                @test fitted.extrapolation == policy
                @test fitted.rates ≈ rates atol = 1.0e-6
                @test zero(fitted, horizon) == zero(
                    Yield.MonotoneConvex(
                        fitted.rates, times;
                        extrapolation = policy
                    ), horizon
                )
            end
        end
    end

    @testset "validation, degenerate grids, and numeric types" begin
        for f in (Inf, -Inf, NaN)
            @test_throws ArgumentError Yield.FlatForwardAt(f)
        end
        @test_throws ArgumentError (@set fixed.forward = Inf)
        for policy in (:extension, :unknown)
            @test_throws ArgumentError Yield.MonotoneConvex(rates, times; extrapolation = policy)
            @test_throws ArgumentError fit(
                Yield.MonotoneConvex(),
                ZCBYield.(rates, times); extrapolation = policy
            )
        end
        for t in (0.0, 2.0), policy in policies
            c = Yield.MonotoneConvex([0.03], [t]; extrapolation = policy)
            @test rate(zero(c, 0.0)) == 0.03
            f = policy === fixed ? fixed.forward : 0.03
            @test rate(forward(c, t, 100.0)) ≈ f
        end
        for T in (Float32, BigFloat)
            c = Yield.MonotoneConvex(T.(rates), T.(times); extrapolation = Yield.FlatForwardAt(T(0.035)))
            @test rate(zero(c, T(horizon))) isa T
            @test Yield.instantaneous_forward(c, T(horizon)) isa T
        end
        negative = Yield.MonotoneConvex(rates, times; extrapolation = Yield.FlatForwardAt(-0.01))
        @test rate(forward(negative, tn, horizon)) ≈ -0.01
    end
end
