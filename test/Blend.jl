# Blends of curves in discount-factor, zero-rate and forward-rate space. Uses the curve types defined
# at the top of log_discount.jl (`__HalfFlatCurve`, `__NegativeDiscountCurve`).
import Random

@testset "Blend" begin
    flat(r) = Yield.Constant(Continuous(r))
    a, b = flat(0.02), flat(0.05)
    lin = ZeroRateCurve([0.02, 0.025, 0.03, 0.031], [1.0, 2.0, 5.0, 10.0], Spline.Linear())
    ns = Yield.NelsonSiegel(1.0, 0.05, -0.02, 0.01)
    DF, ZR = Yield.DiscountFactors(), Yield.ZeroRates()
    L(c, t) = Yield.__log_discount(c, t)
    grade(t) = clamp((8 - t) / 6, 0, 1)   # 1 until 2, 0 from 8

    @testset "discount factors: zero weights, overflow and underflow" begin
        # curves with L(1) = 0 against 1000 (the mixture underflows) or 700 (it doesn't): ∂D/∂w is
        # exact at the end weights, and ∂L/∂w where it is representable (-Inf where it overflows)
        for Lb in (1000.0, 700.0), w0 in (0.0, 0.3, 1.0)
            m(w) = Yield.Blend(flat(0.0), flat(Lb), w, DF)
            D, Db = w0 + (1 - w0) * exp(-Lb), exp(-Lb)
            @test ForwardDiff.derivative(w -> discount(m(w), 1.0), w0) == 1 - Db
            dL = ForwardDiff.derivative(w -> L(m(w), 1.0), w0)
            iszero(D) ? (@test dL == -Inf) : (@test dL ≈ -(1 - Db) / D rtol = 1.0e-14)
        end
        # and the second derivative, d²L/dw² = (D_a - D_b)²/D², at L = (0, 5)
        for w0 in (0.0, 0.3, 1.0)
            Db = exp(-5.0)
            d2 = ForwardDiff.derivative(u -> ForwardDiff.derivative(w -> L(Yield.Blend(flat(0.0), flat(5.0), w, DF), 1.0), u), w0)
            @test d2 ≈ (1 - Db)^2 / (w0 + (1 - w0) * Db)^2 rtol = 1.0e-13
        end
        # a curve whose discount factor overflows, with a tiny share; the reference is exact for these
        # inputs, and the error is the conditioning of an L of 720
        big_small = Yield.Blend((flat(-720.0), flat(0.0)), (1.0e-310, 1 - 1.0e-310), DF)
        @test discount(big_small, 1.0) ≈ 493.07009302638005 rtol = 8 * 720 * eps()
        # weighted discount factors that each underflow, with a representable sum
        tiny = Yield.Blend((flat(745.0), flat(745.0)), (0.5, 0.5), DF)
        @test discount(tiny, 1.0) == 5.0e-324 == discount(tiny, 0.0, 1.0)
        # identical curves are the curve itself, bit for bit
        for c in (lin, ns), t in (0.5, 3.0, 12.0)
            m = Yield.Blend((c, c), (0.3, 0.7), DF)
            @test discount(m, t) === discount(c, t)
            @test L(m, t) === L(c, t)
        end
        # an interval from 0 is `discount(b, t)`, also in its derivatives
        mix(w) = Yield.Blend(lin, ns, w, DF)
        @test discount(mix(0.4), 0.0, 7.0) === discount(mix(0.4), 7.0)
        @test ForwardDiff.derivative(w -> discount(mix(w), 0.0, 7.0), 0.4) === ForwardDiff.derivative(w -> discount(mix(w), 7.0), 0.4)
        @test discount(mix(0.4), 0.0) === 1.0
        @test L(mix(0.4), 0.0) === 0.0
    end

    @testset "the origin in mixed precision" begin
        # Float32 weights on Float64 curves, and Float64 weights on Float32 curves
        cs64 = (lin, ns, a)
        cs32 = (flat(0.02f0), flat(0.03f0), flat(0.05f0))
        for space in (DF, ZR, Yield.ForwardRates(), Yield.ForwardRates(1.0)), (cs, ws) in ((cs64, (0.1f0, 0.2f0, 0.7f0)), (cs32, (0.1, 0.2, 0.7)))
            m = Yield.Blend(cs, ws, space)
            @test discount(m, 0.0) == 1
            @test iszero(L(m, 0.0))
            @test rate(zero(m, 0.0)) ≈ rate(Yield.instantaneous_forward(m, 0.0)) rtol = 1.0e-7
        end
    end

    @testset "validation" begin
        @test_throws ArgumentError Yield.Blend((a, b), (1.0e16, -1.0e16), ZR)
        @test_throws ArgumentError Yield.Blend((a, b), (NaN, 1.0), ZR)
        @test_throws ArgumentError Yield.Blend((a, b), (Inf, -Inf), ZR)
        @test_throws ArgumentError Yield.Blend((a, b), (0.3, 0.6), ZR)
        @test_throws DimensionMismatch Yield.Blend((a, b), (1.0,), ZR)
        @test_throws ArgumentError Yield.Blend(a, b, 1.2, DF)
        @test_throws ArgumentError Yield.Blend((a, b, ns), (1.2, -0.1, -0.1), DF)
        for p in (0.0, -1.0, NaN, Inf)
            @test_throws ArgumentError Yield.ForwardRates(p)
        end
        @test_throws ArgumentError Yield.Blend(a, b, grade, Yield.ForwardRates())
        # weights that vary with tenor are checked where they are read
        bad = Yield.Blend((a, b), (t -> 0.5, t -> 0.6), ZR)
        @test_throws ArgumentError discount(bad, 1.0)
        @test_throws ArgumentError discount(Yield.Blend(a, b, t -> 1.5, DF), 1.0)
        # extrapolating zero-rate weights are allowed; a rounding error in the sum is too
        @test discount(Yield.Blend((a, b, lin), (0.1, 0.2, 0.7), ZR), 5.0) isa Float64
        # vectors are copied, and the stored ones can't be changed
        cs, ws = [a, b], [0.5, 0.5]
        m = Yield.Blend(cs, ws, DF)
        d = discount(m, 10.0)
        ws[1], cs[2] = 0.9, flat(0.09)
        @test discount(m, 10.0) === d
        @test_throws Exception (m.weights[1] = 0.9)
        # rebuilding through Accessors checks again
        m2 = Yield.Blend((a, b), (0.4, 0.6), ZR)
        @test_throws ArgumentError Accessors.@set m2.weights = (0.4, 0.4)
    end

    @testset "a single curve" begin
        h = __HalfFlatCurve()   # D(t) = 0.5·exp(-0.03t)
        mh = Yield.Blend((h,), (1.0,), DF)
        @test discount(mh, 0.0) == 0.5
        @test discount(mh, 0.0, 3.0) ≈ exp(-0.09) rtol = 1.0e-15
        @test zero(mh, 0.0) == zero(h, 0.0)
        @test discount(Yield.Blend((h,), (1.0,), ZR), 3.0) ≈ discount(h, 3.0) rtol = 1.0e-15
        # a forward blend is built from the curve's interval factors: it starts at 1
        fh = Yield.Blend((h,), (1.0,), Yield.ForwardRates())
        @test discount(fh, 0.0) == 1
        @test discount(fh, 3.0) ≈ discount(h, 0.0, 3.0) rtol = 1.0e-15
        # a log-native curve is itself in every space
        for space in (DF, ZR, Yield.ForwardRates(), Yield.ForwardRates(1.0)), t in (0.0, 0.7, 4.0, 12.5)
            @test discount(Yield.Blend((lin,), (1.0,), space), t) === discount(lin, t)
        end
    end

    @testset "a weight that varies with tenor, held constant" begin
        for space in (DF, ZR, Yield.ForwardRates(1.0)), t in (0.4, 3.0, 9.5)
            @test discount(Yield.Blend(lin, ns, t -> 0.3, space), t) ≈ discount(Yield.Blend(lin, ns, 0.3, space), t) rtol = 1.0e-14
        end
        # its limit at infinity is unknown
        @test isnan(discount(Yield.Blend(lin, ns, t -> 0.3, DF), Inf))
        @test discount(Yield.Blend(lin, ns, 0.3, DF), Inf) == 0.0
    end

    @testset "discount factors" begin
        m = Yield.Blend(lin, ns, 0.3, DF)
        cfs, times = [5.0, 5.0, 105.0], [1.0, 2.5, 7.0]
        @test pv(m, cfs, times) ≈ 0.3 * pv(lin, cfs, times) + 0.7 * pv(ns, cfs, times) rtol = 1.0e-14
        # a dominant curve with a tiny share, against a high-precision reference
        tiny_share = Yield.Blend((flat(0.01), flat(50.0)), (1.0e-10, 1 - 1.0e-10), DF)
        ref = setprecision(BigFloat, 256) do
            -log(big(1.0e-10) * exp(-big(0.01)) + (1 - big(1.0e-10)) * exp(-big(50.0)))
        end
        @test L(tiny_share, 1.0) ≈ ref rtol = 4eps()
        # both curves' discount factors underflow; the interval doesn't
        far = Yield.Blend(flat(0.9), flat(1.0), 0.5, DF)
        @test discount(far, 1000.0) == 0.0
        @test discount(far, 999.0, 1000.0) ≈ exp(-0.9) rtol = 1.0e-12
        # the expected discount factor over scenarios
        rng = Random.Xoshiro(1)
        paths = simulate(ShortRate.Vasicek(0.1, 0.03, 0.01, Continuous(0.03)); n_scenarios = 20, timestep = 0.5, horizon = 10.0, rng)
        mean_curve = Yield.Blend(collect(paths), fill(1 / 20, 20), DF)
        @test discount(mean_curve, 6.3) ≈ sum(p -> discount(p, 6.3), paths) / 20 rtol = 1.0e-14
        # rebased to s, the mixture's weights are each curve's share of the discount factor at s
        s, t = 2.0, 3.5
        p = (0.3 * discount(lin, s), 0.7 * discount(ns, s)) ./ discount(m, s)
        @test discount(Yield.ForwardStarting(m, s), t) ≈
            p[1] * discount(Yield.ForwardStarting(lin, s), t) + p[2] * discount(Yield.ForwardStarting(ns, s), t) rtol = 1.0e-14
        # the forward is each curve's, weighted by its share of the discount factor
        f = rate(Yield.instantaneous_forward(m, 4.0))
        q = 0.3 * discount(lin, 4.0) / discount(m, 4.0)
        @test f ≈ q * rate(Yield.instantaneous_forward(lin, 4.0)) + (1 - q) * rate(Yield.instantaneous_forward(ns, 4.0)) rtol = 1.0e-14
        # the gradient in N independent weights is (Dᵢ - D)/Σw
        g = ForwardDiff.gradient(w -> discount(Yield.Blend((lin, ns, a), (w[1], w[2], w[3]), DF), 4.0), [0.2, 0.3, 0.5])
        D = discount(Yield.Blend((lin, ns, a), (0.2, 0.3, 0.5), DF), 4.0)
        @test g ≈ [discount(lin, 4.0), discount(ns, 4.0), discount(a, 4.0)] .- D rtol = 1.0e-13
    end

    @testset "zero rates" begin
        for w in (0.3, -0.4, 1.5), (x, y) in ((a, b), (lin, ns))
            m, c = Yield.Blend(x, y, w, ZR), w * x + (1 - w) * y
            for t in (0.0, 0.5, 4.0, 12.0)
                @test discount(m, t) === discount(c, t)
                @test rate(zero(m, t)) === rate(zero(c, t))
                @test discount(m, 1.5, t) === discount(c, 1.5, t)
            end
            @test discount(m, Inf) === discount(c, Inf)
        end
        # graded: z(t) = w(t)·z₁(t) + (1 - w(t))·z₂(t)
        m = Yield.Blend(lin, ns, grade, ZR)
        for t in (0.5, 3.0, 9.0)
            @test rate(zero(m, t)) ≈ grade(t) * rate(zero(lin, t)) + (1 - grade(t)) * rate(zero(ns, t)) rtol = 1.0e-14
            @test discount(m, t) ≈ exp(-rate(zero(m, t)) * t) rtol = 1.0e-14
        end
        @test isnan(discount(m, Inf))
        @test discount(m, 5.0, 5.0) == 1
    end

    @testset "forward rates" begin
        m = Yield.Blend(lin, ns, grade, Yield.ForwardRates(1.0))
        # each period's forward is the curves' forwards, weighted at the period's start
        for k in 0:9
            @test rate(forward(m, k, k + 1)) ≈ grade(k) * rate(forward(lin, k, k + 1)) + (1 - grade(k)) * rate(forward(ns, k, k + 1)) rtol = 1.0e-12
        end
        # L on and off the grid, and at 0
        Lsum(t) = sum(k -> grade(k) * Yield.__log_interval(lin, k, min(k + 1, t)) + (1 - grade(k)) * Yield.__log_interval(ns, k, min(k + 1, t)), 0:(ceil(Int, t) - 1))
        @test L(m, 3.0) ≈ Lsum(3.0) rtol = 1.0e-14
        @test L(m, 3.6) ≈ Lsum(3.6) rtol = 1.0e-14
        @test L(m, 0.4) ≈ Lsum(0.4) rtol = 1.0e-14
        @test discount(m, 0.0) == 1
        # at a grid point the time derivative is the right-hand one, the forward of the period starting there
        f3 = grade(3) * rate(Yield.instantaneous_forward(lin, 3.0)) + (1 - grade(3)) * rate(Yield.instantaneous_forward(ns, 3.0))
        @test ForwardDiff.derivative(t -> L(m, t), 3.0) ≈ f3 rtol = 1.0e-14
        @test rate(Yield.instantaneous_forward(m, 3.0)) ≈ f3 rtol = 1.0e-14
        # intervals spanning periods, and reversed
        @test Yield.__log_interval(m, 1.5, 4.5) ≈ L(m, 4.5) - L(m, 1.5) rtol = 1.0e-13
        @test Yield.__log_interval(m, 4.5, 1.5) == -Yield.__log_interval(m, 1.5, 4.5)
        # a curve that is log-native takes its L once per grid point, with the same result as its
        # interval rule (here through `ForwardStarting` at 0, which isn't log-native)
        fs(c) = Yield.ForwardStarting(c, 0.0)
        for t in (0.4, 3.0, 9.5)
            @test discount(m, t) === discount(Yield.Blend(fs(lin), fs(ns), grade, Yield.ForwardRates(1.0)), t)
        end
        # a signed curve has interval factors but no log-discount of its own
        neg = __NegativeDiscountCurve()
        for t in (0.4, 3.7)
            numbers = discount(Yield.Blend((neg, a), (0.4, 0.6), Yield.ForwardRates()), t)
            @test numbers ≈ exp(-(0.4 * 0.03 + 0.6 * 0.02) * t) rtol = 1.0e-14
            @test discount(Yield.Blend((neg, a), (t -> 0.4, t -> 0.6), Yield.ForwardRates(1.0)), t) ≈ numbers rtol = 1.0e-14
        end
        # with number weights the period doesn't matter, and the blend is the zero-rate blend
        for p in (0.7, 1.0), t in (0.0, 2.5, 11.0)
            @test discount(Yield.Blend(lin, ns, 0.3, Yield.ForwardRates(p)), t) === discount(Yield.Blend(lin, ns, 0.3, ZR), t)
        end
        @test discount(Yield.Blend(lin, ns, 0.3, Yield.ForwardRates()), Inf) == 0.0
        @test isnan(discount(m, Inf))
    end

    @testset "limits at infinity" begin
        c4, c2, cm2 = flat(0.04), flat(0.02), flat(-0.02)
        mix(x, y) = Yield.Blend((x, y), (0.5, 0.5), DF)
        @test discount(mix(c4, c2), Inf) == 0.0
        @test discount(mix(c4, c2) - c2, Inf) ≈ 0.5 rtol = 1.0e-15
        @test discount(mix(c4, cm2), Inf) == Inf
        @test discount(mix(c4, cm2) + c2, Inf) ≈ 0.5 rtol = 1.0e-15
        # tied slopes combine their intercepts: L = t/16 - 1/32 and t/16
        knot = ZeroRateCurve([0.03125, 0.046875], [1.0, 2.0], Spline.Linear())
        @test discount(mix(knot, flat(0.0625)) - flat(0.0625), Inf) ≈ 0.5 * exp(0.03125) + 0.5 rtol = 1.0e-15
        # a curve whose intercept at infinity is unknown: alone it still decides the limit, and
        # against a curve with a lower slope it is negligible
        @test discount(Yield.Blend((ns,), (1.0,), DF), Inf) == 0.0
        @test discount(mix(ns, c2) - c2, Inf) ≈ 0.5 rtol = 1.0e-15
        # a curve with an unknown slope makes the limit unknown
        vas = ShortRate.Vasicek(0.1, 0.03, 0.01, Continuous(0.03))
        @test isnan(discount(mix(vas, c2), Inf))
        # growth faster than t of unknown order, against a quadratic tail: unknown
        u = Yield.TenorShift(Yield.Constant(Continuous(0.0)), (z, t) -> z + Continuous(t / 32))
        q = ZeroRateCurve([1 / 64, 2 / 64], [1.0, 2.0], Spline.Linear(); extrapolation = :linear)
        @test isnan(discount(mix(u, q) - q, Inf))
        # zero-rate blends combine tails linearly
        @test discount(Yield.Blend(c4, cm2, 0.25, ZR), Inf) == Inf
        @test discount(Yield.Blend(c4, cm2, 0.5, ZR), Inf) == 0.0
        # a tie that a bump would split, and a zero share on a curve that would dominate, have no derivative
        split(r) = discount(mix(flat(r), flat(0.03)) - flat(0.03), Inf)
        @test split(0.03) ≈ 1.0 rtol = 1.0e-15
        @test_throws "no derivative" ForwardDiff.derivative(split, 0.03)
        @test ForwardDiff.derivative(split, 0.04) == 0.0
        @test_throws "no derivative" ForwardDiff.derivative(w -> discount(Yield.Blend(flat(0.01), flat(0.03), w, DF), Inf), 0.0)
        @test ForwardDiff.derivative(w -> discount(Yield.Blend(flat(0.05), flat(0.03), w, DF), Inf), 0.0) == 0.0
    end

    @testset "number types and derivatives" begin
        for space in (DF, ZR, Yield.ForwardRates(), Yield.ForwardRates(1.0))
            @test discount(Yield.Blend(flat(0.02f0), flat(0.05f0), 0.3f0, space), 5.0f0) isa Float32
            mb = Yield.Blend(flat(big"0.02"), flat(big"0.05"), big"0.3", space)
            @test discount(mb, big"5.0") isa BigFloat
        end
        @test discount(Yield.Blend(flat(0.02f0), flat(0.05f0), t -> 0.3f0, Yield.ForwardRates(1.0f0)), 5.0f0) isa Float32
        # a curve parameter, against central differences, in each space
        for space in (DF, ZR, Yield.ForwardRates(1.0)), w in (0.4, t -> grade(t))
            f(r) = discount(Yield.Blend(ZeroRateCurve([r, 0.03], [1.0, 5.0], Spline.Linear()), ns, w, space), 3.0)
            h = 1.0e-6
            @test ForwardDiff.derivative(f, 0.02) ≈ (f(0.02 + h) - f(0.02 - h)) / 2h rtol = 1.0e-7
        end
        # a second derivative in time
        m = Yield.Blend(lin, ns, 0.3, DF)
        d1(t) = ForwardDiff.derivative(u -> discount(m, u), t)
        h = 1.0e-5
        @test ForwardDiff.derivative(d1, 2.7) ≈ (d1(2.7 + h) - d1(2.7 - h)) / 2h rtol = 1.0e-6
    end
end
