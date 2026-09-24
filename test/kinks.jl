using ForwardDiff

# Shape-preserving interpolants are only piecewise smooth in their knot rates. These tests pin
# the derivatives at their kinks: MonotoneConvex reports the limit of a centered bump in each
# partial's direction; PCHIP, Akima, and fits that sit on a kink throw.
@testset "derivatives at interpolant kinks" begin
    E(n, i) = Float64.(1:n .== i)
    centred(f, z, v; h = 1.0e-7) = (f(z .+ h .* v) - f(z .- h .* v)) / 2h
    centred(f, z; h = 1.0e-7) = [centred(f, z, E(length(z), i); h) for i in eachindex(z)]
    ts = [1.0, 2.0, 3.0]

    @testset "MonotoneConvex: two adjacent discrete forwards equal (two sectors meet)" begin
        # fᵈ₁ = fᵈ₂ = 2%, fᵈ₃ = 4%: interval 2 has g0 = 0, where sectors (ii) and (iv) meet.
        z = [0.02, 0.02, 0.08 / 3]
        for t in (0.5, 1.5, 2.5)
            f = x -> discount(ZeroRateCurve(x, ts), t)
            g = ForwardDiff.gradient(f, z)
            @test g ≈ centred(f, z) rtol = 1.0e-6
            # the centered limit of two pieces is linear in the direction: sums are parallel
            @test sum(g) ≈ ForwardDiff.derivative(s -> f(z .+ s), 0.0) rtol = 1.0e-12
            # a combined direction seeded as one partial
            v = [1.0, 1.0, 0.0]
            @test ForwardDiff.derivative(s -> f(z .+ s .* v), 0.0) ≈ g[1] + g[2] rtol = 1.0e-12
            @test ForwardDiff.derivative(s -> f(z .+ s .* v), 0.0) ≈ centred(f, z, v) rtol = 1.0e-6
        end
        # the one-sided derivatives do differ at t = 1.5: this is a kink, not a smooth point
        f = x -> discount(ZeroRateCurve(x, ts), 1.5)
        h = 1.0e-7
        @test abs((f(z .+ h .* E(3, 2)) - f(z)) / h - (f(z) - f(z .- h .* E(3, 2))) / h) > 0.1
    end

    @testset "MonotoneConvex: flat intervals" begin
        z = fill(0.03, 3)
        for t in (0.5, 1.5, 2.5)
            f = x -> discount(ZeroRateCurve(x, ts), t)
            g = ForwardDiff.gradient(f, z)
            @test g ≈ centred(f, z) rtol = 1.0e-6 atol = 1.0e-9
            # a parallel shift keeps the curve flat: its derivative is exact, -t·D(t)
            @test ForwardDiff.derivative(s -> f(z .+ s), 0.0) ≈ -t * exp(-0.03 * t) rtol = 1.0e-12
            # chunk sizes evaluate the same partials
            for chunk in (1, 2, 3)
                cfg = ForwardDiff.GradientConfig(f, z, ForwardDiff.Chunk{chunk}())
                @test ForwardDiff.gradient(f, z, cfg) == g
            end
        end
        # Inside the curve a flat interval has no gradient: single-knot sensitivities do not add
        # up to the parallel one (documented).
        f = x -> discount(ZeroRateCurve(x, ts), 1.5)
        @test sum(ForwardDiff.gradient(f, z)) ≈ -1.3891838408 rtol = 1.0e-8
        @test ForwardDiff.derivative(s -> f(z .+ s), 0.0) ≈ -1.5 * exp(-0.045) rtol = 1.0e-12
        # a flat curve on a dense grid, where rounding leaves the forwards slightly unequal
        qs = collect(0.25:0.25:30.0)
        zq = fill(0.0437, length(qs))
        fq = x -> discount(ZeroRateCurve(x, qs), 15.1)
        @test ForwardDiff.gradient(fq, zq) ≈ centred(fq, zq) atol = 1.0e-8
    end

    @testset "MonotoneConvex: the positivity collar binds exactly" begin
        # Dyadic rates: fᵈ = (a, 3a, 4a) and the node forward at t = 1 is (a + 3a)/2 = 2a, the
        # collar bound 2·min(fᵈ₁, fᵈ₂), exactly.
        a = 1 / 64
        z = [a, 2a, 8a / 3]
        f, fᵈ = FinanceModels.Yield.__monotone_convex_fs(z, ts)
        @test f[2] == 2 * fᵈ[1] == 2a
        for t in (0.5, 1.5)
            v = x -> discount(ZeroRateCurve(x, ts), t)
            g = ForwardDiff.gradient(v, z)
            # The end intervals sit on sector boundaries where the second derivative jumps, so a
            # centered difference has an O(h) error; Richardson extrapolation removes it.
            @test g ≈ 2 .* centred(v, z; h = 5.0e-8) .- centred(v, z; h = 1.0e-7) rtol = 1.0e-6
            @test sum(g) ≈ ForwardDiff.derivative(s -> v(z .+ s), 0.0) rtol = 1.0e-12
        end
    end

    @testset "MonotoneConvex: near a kink the derivative is the local one" begin
        z = [0.03, 0.03 + 1.0e-4, 0.03]
        f = x -> discount(ZeroRateCurve(x, ts), 1.5)
        @test ForwardDiff.gradient(f, z) ≈ centred(f, z; h = 1.0e-8) rtol = 1.0e-5
    end

    @testset "MonotoneConvex: unresolvable ties throw" begin
        f = x -> discount(ZeroRateCurve(x, ts), 1.5)
        # a 0% curve: the collar's bound collapses to zero at every node
        @test_throws "positivity collar" ForwardDiff.gradient(f, zeros(3))
        # second derivatives do not exist at a kink
        @test_throws "second derivatives" ForwardDiff.hessian(f, fill(0.03, 3))
        # the value is the primal value
        d = f(ForwardDiff.Dual.(fill(0.03, 3), 1.0))
        @test ForwardDiff.value(d) == f(fill(0.03, 3))
    end

    @testset "PCHIP and Akima throw at kinks" begin
        a = 1 / 64
        cases = (
            (Spline.PCHIP(), fill(0.03, 3)),            # flat: zero secant slopes
            (Spline.PCHIP(), [0.02, 0.03, 0.03]),       # one flat segment
            (Spline.Akima(), fill(0.03, 3)),
            (Spline.Akima(), [a, 2a, 3a]),              # a straight run: equal secant slopes
        )
        for (sp, z) in cases
            f = x -> discount(ZeroRateCurve(x, ts, sp), 1.5)
            @test_throws "switches formula" ForwardDiff.gradient(f, z)
            # a direction that keeps the configuration (here a parallel shift) is smooth
            @test ForwardDiff.derivative(s -> f(z .+ s), 0.0) ≈ centred(f, z, ones(3)) rtol = 1.0e-6
        end
        # away from kinks, exact
        for sp in (Spline.PCHIP(), Spline.Akima())
            z = [0.02, 0.026, 0.029, 0.037, 0.038]
            f = x -> discount(ZeroRateCurve(x, [1.0, 2.0, 3.0, 4.0, 5.0], sp), 2.5)
            @test ForwardDiff.gradient(f, z) ≈ centred(f, z) rtol = 1.0e-6
        end
    end

    @testset "Akima: the fallback-slope cutoff" begin
        # Knot 3's slope weight is 3e and the largest weight about 0.019, so DataInterpolations
        # switches knot 3 to its fallback slope near e = 0.019e-9/3, and the curve's value jumps.
        t7 = collect(1.0:7.0)
        knots(e) = [0.02; 0.02 .+ cumsum([0.001, 0.001 + e, 0.01, 0.01 + 2e, 0.02, 0.022])]
        cutoff(e) = FinanceModels.Yield.__kink_quantities(Spline.Akima(), knots(e), t7)[9 + 3]
        lo, hi = 0.0, 1.0e-10
        for _ in 1:200
            mid = (lo + hi) / 2
            (mid == lo || mid == hi) && break
            cutoff(mid) > 0 ? (hi = mid) : (lo = mid)
        end
        f = x -> discount(ZeroRateCurve(x, t7, Spline.Akima()), 2.5)
        # The guard's quantity changes sign exactly where the installed DataInterpolations switches
        # formula: the value jumps between adjacent floating-point e. If this fails,
        # DataInterpolations changed the rule `__kink_quantities(::Sp.Akima, ...)` mirrors.
        @test abs(f(knots(hi)) - f(knots(lo))) > 1.0e-4
        @test_throws "switches formula" ForwardDiff.gradient(f, knots(lo))
        @test_throws "switches formula" ForwardDiff.gradient(f, knots(hi))
        # Either side, clear of the cutoff (and of knot 1's, at 1.5 times it), the derivative is
        # the local one: compare with a centered difference in high precision.
        for e in (0.8lo, 1.2hi)
            z = knots(e)
            ref = Float64.(centred(f, big.(z); h = big(1.0e-30)))
            @test ForwardDiff.gradient(f, z) ≈ ref rtol = 1.0e-5
        end
    end

    @testset "a first knot at t = 0" begin
        t0 = [0.0, 1.0, 2.0]
        for sp in (Spline.PCHIP(), Spline.Akima())
            z = [0.02, 0.026, 0.035]
            f = x -> discount(ZeroRateCurve(x, t0, sp), 0.5)
            @test ForwardDiff.gradient(f, z) ≈ centred(f, z) rtol = 1.0e-6
        end
        z = fill(0.03, 3)
        f = x -> discount(ZeroRateCurve(x, t0), 0.5)
        @test ForwardDiff.gradient(f, z) ≈ centred(f, z) atol = 1.0e-9
        qs = collect(0.0:0.25:30.0)
        zq = fill(0.0437, length(qs))
        fq = x -> discount(ZeroRateCurve(x, qs), 15.1)
        @test ForwardDiff.gradient(fq, zq) ≈ centred(fq, zq) atol = 1.0e-8
    end

    @testset "fits that sit on a kink throw" begin
        prices = exp.(-0.03 .* ts)
        for sp in (Spline.MonotoneConvex(), Spline.PCHIP(), Spline.Akima())
            f = p -> discount(fit(sp, ZCBPrice.(p, ts)), 1.5)
            @test_throws "switches shape" ForwardDiff.gradient(f, prices)
            @test f(prices) ≈ exp(-0.045) rtol = 1.0e-6   # the primal fit is fine
        end
        # a smooth interpolant differentiates the same flat quotes
        g = ForwardDiff.gradient(p -> discount(fit(Spline.Cubic(), ZCBPrice.(p, ts)), 1.5), prices)
        @test all(isfinite, g)
    end
end
