using ForwardDiff

# Exercise runtime Symbol selection, not just constant-folded keyword literals.
monotone_with_policy(r, t, policy::Symbol) = Yield.MonotoneConvex(r, t; extrapolation = policy)
zero_curve_with_policy(r, t, policy::Symbol) = ZeroRateCurve(r, t; extrapolation = policy)
spline_with_policy(d, r, t, policy::Symbol) = Yield.Spline(d, t, r; extrapolation = policy)
spline_zero_curve_with_policy(d, r, t, policy::Symbol) = ZeroRateCurve(r, t, d; extrapolation = policy)

@testset "shared curve tails" begin
    rates = [0.02, 0.025, 0.03, 0.035, 0.04, 0.042]
    times = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0]
    tn, zn, horizon = last(times), last(rates), 100.0
    fixed = Yield.FlatForwardAt(Continuous(0.035))
    policies = (:flat_forward, :flat_zero, :linear, fixed)
    descriptors = (
        Spline.Linear(), Spline.Quadratic(), Spline.Cubic(),
        Spline.PCHIP(), Spline.Akima(), Spline.BSpline(3), Spline.MonotoneConvex(),
    )

    @testset "concrete construction through the AD path" begin
        for T in (Float32, Float64), U in (Float32, Float64, BigFloat)
            r, t = T.(rates), U.(times)
            @test (@inferred Yield.MonotoneConvex(r, t)) isa Yield.MonotoneConvex
            @test (@inferred ZeroRateCurve(r, t)) isa Yield.MonotoneConvex
            for policy in (:flat_forward, :flat_zero, :linear)
                @test (@inferred monotone_with_policy(r, t, policy)).extrapolation === policy
                @test (@inferred zero_curve_with_policy(r, t, policy)).extrapolation === policy
            end
        end
        # DataInterpolations-backed curves: every policy Symbol, including `:extension`,
        # builds the same concrete type. (`PolynomialSpline` stores its order as a runtime
        # value, so its interpolant type is not inferable from the descriptor type.)
        for d in (Spline.PCHIP(), Spline.Akima(), Spline.BSpline(3))
            for T in (Float32, Float64), U in (Float64, BigFloat)
                r, t = T.(rates), U.(times)
                types = map((:flat_forward, :flat_zero, :linear, :extension)) do policy
                    c = @inferred spline_with_policy(d, r, t, policy)
                    z = @inferred spline_zero_curve_with_policy(d, r, t, policy)
                    @test z.extrapolation === policy
                    typeof(c)
                end
                @test allequal(types)
                @test isconcretetype(first(types))
            end
            # policies known at compile time
            @test (@inferred Yield.Spline(d, times, rates)) isa Yield.Spline
            @test (@inferred Yield.Spline(d, times, rates; extrapolation = :extension)) isa Yield.Spline
            @test (@inferred Yield.Spline(d, times, rates; extrapolation = :flat_zero)) isa Yield.Spline
            @test (@inferred Yield.Spline(d, times, rates; extrapolation = :linear)) isa Yield.Spline
            @test (@inferred Yield.Spline(d, times, rates; extrapolation = fixed)) isa Yield.Spline
        end
        function pv(r)
            c = @inferred ZeroRateCurve(r, times)
            return sum(discount(c, t) for t in range(1.0, 100.0; length = 80))
        end
        @test all(isfinite, ForwardDiff.gradient(pv, rates))
        # A zero-valued coefficient may still carry knot-rate derivatives. Both
        # reciprocal and linear terms must remain active for their selected policy.
        tm = times[end - 1]
        Δ = tn - tm
        for policy in (:flat_forward, :linear)
            z(r) = rate(zero(ZeroRateCurve(r, times, Spline.Linear(); extrapolation = policy), horizon))
            flat = fill(0.03, length(rates))
            ad = ForwardDiff.gradient(z, flat)
            if policy === :linear
                weight = (horizon - tn) / Δ
                @test ad ≈ [zeros(length(rates) - 2); -weight; 1 + weight]
            else
                # z(H) = f + (zₙ - f)tₙ/H with f = (zₙtₙ - zₘtₘ)/Δ, the last discrete forward
                tail = 1 - tn / horizon
                @test ad ≈ [zeros(length(rates) - 2); -tail * tm / Δ; tn / horizon + tail * tn / Δ]
            end
        end
    end

    @testset "boundary and AD: $descriptor" for descriptor in descriptors
        base = ZeroRateCurve(rates, times, descriptor)
        # A second-order one-sided difference measures the INTERIOR zero slope.
        h = 1.0e-3
        left_slope = (
            3rate(zero(base, tn)) - 4rate(zero(base, tn - h)) +
                rate(zero(base, tn - 2h))
        ) / (2h)
        for policy in policies
            curve = ZeroRateCurve(rates, times, descriptor; extrapolation = policy)
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
                # `:flat_forward` anchors DataInterpolations-backed curves on the last discrete
                # forward; MonotoneConvex keeps its boundary instantaneous forward.
                f = if policy isa Yield.FlatForwardAt
                    policy.forward
                elseif curve isa Yield.MonotoneConvex
                    zn + tn * left_slope
                else
                    (zn * tn - rates[end - 1] * times[end - 1]) / (tn - times[end - 1])
                end
                @test rate(forward(curve, tn, horizon)) ≈ f atol = 3.0e-8
                @test rate(zero(curve, 1.0e5)) ≈ f + (zn - f) * tn / 1.0e5 atol = 3.0e-8
                @test rate(zero(curve, Inf)) ≈ f atol = 3.0e-8
            end
            # discount factors are continuous at the last knot under every policy
            @test discount(curve, tn + 1.0e-9) ≈ discount(curve, tn) rtol = 1.0e-9

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
                @test Yield.instantaneous_forward(curve, tn) == last(curve._f)
            end
        end
        # Different fixed assumptions change only the tail and remain differentiable.
        zf(f) = rate(
            zero(
                ZeroRateCurve(
                    rates, times, descriptor;
                    extrapolation = Yield.FlatForwardAt(Continuous(f))
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
                    @test rebuilt == c && isequal(rebuilt, c) && hash(rebuilt) == hash(c)
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
                moved = @set c.tenors = times .* 2
                @test moved.extrapolation == policy && knot_tenors(moved) == times .* 2
                @test_throws ArgumentError CB.setproperties(c, (_tail = nothing,))
                @test_throws ArgumentError (@set c.extrapolation = :unknown)
            end
        end
        sampled = ZeroRateCurve(
            Yield.Constant(Continuous(0.03)), times;
            extrapolation = fixed
        )
        @test sampled.extrapolation == fixed
        @test rate(forward(sampled, tn, horizon)) ≈ fixed.forward
        a = ZeroRateCurve(rates, times; extrapolation = Yield.FlatForwardAt(Continuous(big"0.035")))
        b = ZeroRateCurve(rates, times; extrapolation = Yield.FlatForwardAt(Continuous(big"0.035")))
        @test a == b && isequal(a, b) && hash(a) == hash(b)
        @test a != (@set a.extrapolation = Yield.FlatForwardAt(Continuous(0.04)))
        @test Yield.FlatForwardAt(Continuous(-0.0)) == Yield.FlatForwardAt(Continuous(0.0))
        @test !isequal(Yield.FlatForwardAt(Continuous(-0.0)), Yield.FlatForwardAt(Continuous(0.0)))
        for policy in policies
            c = Yield.MonotoneConvex(rates, times; extrapolation = policy)
            @test c == Yield.MonotoneConvex(copy(rates), copy(times); extrapolation = policy)
            @test c != (@set c.rates[end] = c.rates[end] + 0.001)
            @test c != (@set c.tenors[end] = c.tenors[end] + 1)
            @test c != (@set c.extrapolation = Yield.FlatForwardAt(Continuous(0.045)))
            @test Set([c, deepcopy(c)]) == Set([c])
        end
        negative_zero = Yield.MonotoneConvex([-0.0], [1.0])
        positive_zero = Yield.MonotoneConvex([0.0], [1.0])
        @test negative_zero == positive_zero
        @test !isequal(negative_zero, positive_zero)
        @test length(Set([negative_zero, positive_zero])) == 2
    end

    @testset "explicit compounding conventions" begin
        @test Yield.FlatForwardAt(Continuous(0.035)) == fixed
        annual = Yield.FlatForwardAt(Periodic(0.035, 1))
        @test annual.forward ≈ log1p(0.035)
        @test annual != fixed
        @test annual.forward == rate(zero(Yield.Constant(0.035), 1.0))
        c = Yield.MonotoneConvex(rates, times; extrapolation = annual)
        @test rate(forward(c, tn, horizon)) ≈ log1p(0.035)
        @test ForwardDiff.derivative(f -> Yield.FlatForwardAt(Periodic(f, 2)).forward, 0.035) ≈ 1 / (1 + 0.035 / 2)
        for r in (Continuous(Inf), Continuous(NaN))
            @test_throws ArgumentError Yield.FlatForwardAt(r)
        end
        # A bare number has no stated convention (`Yield.Constant(0.035)` reads it as
        # annual effective), so it must be wrapped in a rate.
        for x in (0.035, 1, big"0.035", -0.01)
            err = try
                Yield.FlatForwardAt(x)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("Continuous(", sprint(showerror, err))
            @test occursin("Periodic(", sprint(showerror, err))
        end
        # the stored field is continuously compounded, and Accessors rebuilds it that way
        @test (@set fixed.forward = 0.04) == Yield.FlatForwardAt(Continuous(0.04))
    end

    @testset "fits retain the policy and native MonotoneConvex interface" begin
        qs = ZCBYield.(Continuous.(rates), times)
        # Bootstrap rejects MonotoneConvex under every tail policy and points to the loss fit.
        for policy in policies
            @test_throws "fit(Spline.MonotoneConvex(), quotes)" fit(
                Spline.MonotoneConvex(), qs, Fit.Bootstrap(); extrapolation = policy
            )
        end
        for method in (Fit.Loss(abs2),)
            baseline = fit(Spline.MonotoneConvex(), qs, method)
            for policy in policies
                c = fit(Spline.MonotoneConvex(), qs, method; extrapolation = policy)
                @test c isa Yield.MonotoneConvex
                @test c.rates == baseline.rates && c.tenors == baseline.tenors
                @test c.extrapolation == policy
                @test isfinite(Yield.instantaneous_forward(c, horizon))
                @test zero(c, horizon) == zero(
                    Yield.MonotoneConvex(
                        c.rates, c.tenors;
                        extrapolation = policy
                    ), horizon
                )
            end
        end
        for descriptor in (Spline.Linear(), Spline.MonotoneConvex())
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
        # Bounded tails retain their limits at infinity; linear forwards avoid
        # prematurely overflowing 2t at very large finite horizons.
        for policy in (:flat_forward, :flat_zero, fixed)
            c = Yield.MonotoneConvex(rates, times; extrapolation = policy)
            @test isfinite(rate(zero(c, Inf)))
            @test rate(zero(c, Inf)) == Yield.instantaneous_forward(c, Inf)
        end
        c = Yield.MonotoneConvex(rates, times; extrapolation = :linear)
        @test isfinite(Yield.instantaneous_forward(c, floatmax(Float64)))

        # Every policy at t = Inf. Bounded tails converge to their anchor; `:linear`
        # diverges with the sign of the boundary slope, or stays flat when the slope is
        # zero (previously 0 * Inf = NaN).
        flat_rates = fill(0.03, length(times))
        for d in (Spline.Linear(), Spline.Cubic(), Spline.PCHIP(), Spline.MonotoneConvex())
            for policy in policies
                c = ZeroRateCurve(rates, times, d; extrapolation = policy)
                expected = if policy === :flat_zero
                    zn
                elseif policy === :linear
                    Inf
                elseif policy isa Yield.FlatForwardAt
                    policy.forward
                elseif c isa Yield.MonotoneConvex
                    Yield.instantaneous_forward(c, tn)
                else
                    (zn * tn - rates[end - 1] * times[end - 1]) / (tn - times[end - 1])
                end
                @test rate(zero(c, Inf)) ≈ expected
                # (A one-knot MonotoneConvex has an exactly zero boundary slope; on a longer
                # flat grid its discrete forwards can differ from 0.03 by rounding.)
                flat = if d == Spline.MonotoneConvex()
                    Yield.MonotoneConvex([0.03], [tn]; extrapolation = policy)
                else
                    ZeroRateCurve(flat_rates, times, d; extrapolation = policy)
                end
                expected_flat = policy isa Yield.FlatForwardAt ? policy.forward : 0.03
                @test rate(zero(flat, Inf)) ≈ expected_flat
                @test !isnan(rate(zero(flat, Inf)))
                if flat isa Yield.MonotoneConvex
                    @test Yield.instantaneous_forward(flat, Inf) ≈ expected_flat
                end
            end
            if d != Spline.MonotoneConvex()
                # `:extension` keeps DataInterpolations' own continuation of the final
                # piece, even at t = Inf, where a sloped linear piece diverges and a cubic
                # piece can be undefined (Inf - Inf). This is the legacy behavior.
                c = ZeroRateCurve(rates, times, d; extrapolation = :extension)
                DI = FinanceModels.DataInterpolations
                E = DI.ExtrapolationType.Extension
                legacy = if d == Spline.Linear()
                    DI.LinearInterpolation(rates, times; extrapolation = E)
                elseif d == Spline.Cubic()
                    DI.CubicSpline(rates, times; extrapolation = E)
                else
                    DI.PCHIPInterpolation(rates, times; extrapolation = E)
                end
                @test isequal(rate(zero(c, Inf)), legacy(Inf))
                d == Spline.Linear() && @test rate(zero(c, Inf)) == Inf
            end
        end
        # `discount` at t = Inf is the tail's limit: 0 for a positive tail forward, Inf for a
        # negative one, and the discount factor at the last knot for a zero one (previously
        # exp(-0 * Inf) = NaN).
        for d in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            c3(e) = ZeroRateCurve([0.02, 0.025, 0.03], [1.0, 2.0, 3.0], d; extrapolation = e)
            @test discount(c3(Yield.FlatForwardAt(Continuous(0.0))), Inf) == discount(c3(:flat_zero), 3.0)
            @test discount(c3(Yield.FlatForwardAt(Continuous(0.01))), Inf) == 0.0
            @test discount(c3(Yield.FlatForwardAt(Continuous(-0.01))), Inf) == Inf
            @test discount(c3(:flat_forward), Inf) == 0.0
            @test discount(c3(:flat_zero), Inf) == 0.0
        end
        @test discount(ZeroRateCurve([0.02, 0.0], [1.0, 2.0], Spline.Linear(); extrapolation = :flat_zero), Inf) == 1.0
        @test_throws DomainError discount(ZeroRateCurve([0.02, 0.03], [1.0, 2.0]), -Inf)

        falling = Yield.Spline(Spline.Linear(), times, reverse(rates); extrapolation = :linear)
        @test rate(zero(falling, Inf)) == -Inf
        no_boundary = () -> error("must not inspect boundary")
        @test_throws ArgumentError Yield.__build_tail(:extension, tn, zn, no_boundary, no_boundary)
        @test_throws ArgumentError Yield.__build_tail(:unknown, tn, zn, no_boundary, no_boundary)
        for f in (Inf, -Inf, NaN)
            @test_throws ArgumentError Yield.FlatForwardAt(Continuous(f))
        end
        @test_throws ArgumentError (@set fixed.forward = Inf)
        for policy in (:extension, :unknown)
            @test_throws ArgumentError Yield.MonotoneConvex(rates, times; extrapolation = policy)
            @test_throws ArgumentError fit(
                Spline.MonotoneConvex(),
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
            c = Yield.MonotoneConvex(T.(rates), T.(times); extrapolation = Yield.FlatForwardAt(Continuous(T(0.035))))
            @test rate(zero(c, T(horizon))) isa T
            @test Yield.instantaneous_forward(c, T(horizon)) isa T
        end
        negative = Yield.MonotoneConvex(rates, times; extrapolation = Yield.FlatForwardAt(Continuous(-0.01)))
        @test rate(forward(negative, tn, horizon)) ≈ -0.01
    end
end
