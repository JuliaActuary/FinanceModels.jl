using ForwardDiff

# A curve that defines only `discount`, with D(0) ≠ 1: intervals keep D(to)/D(from) semantics.
struct __HalfFlatCurve <: Yield.AbstractYieldModel end
FinanceCore.discount(::__HalfFlatCurve, t) = 0.5 * exp(-0.03 * t)

# A curve with negative discount factors (as a Smith-Wilson fit to arbitrary prices can have): its
# interval factors are still the exact ratio, but it has no log-discount.
struct __NegativeDiscountCurve <: Yield.AbstractYieldModel end
FinanceCore.discount(::__NegativeDiscountCurve, t) = -exp(-0.03 * t)

@testset "Interval factors from cumulative log-discounts" begin
    rates, tenors = [0.02, 0.03, 0.035, 0.04], [1.0, 2.0, 5.0, 10.0]
    zrc_lin = ZeroRateCurve(rates, tenors, Spline.Linear())
    ns = Yield.NelsonSiegel(1.0, 0.05, -0.02, 0.01)
    sw = Yield.SmithWilson([1.0, 5.0, 10.0], [1.0 0 0; 0 1.0 0; 0 0 1.0], [0.98, 0.86, 0.7]; ufr = 0.03, α = 0.1)

    @testset "far-tail intervals no longer underflow to NaN" begin
        # the long-run forward of this Nelson-Siegel curve is β₀ = 5%
        @test discount(ns, 20000.0, 20001.0) ≈ exp(-0.05) rtol = 1.0e-12
        @test discount(Yield.ForwardStarting(ns, 20000.0), 1.0) ≈ exp(-0.05) rtol = 1.0e-12
        # Smith-Wilson converges to its ultimate forward rate
        @test discount(sw, 30000.0, 30001.0) ≈ exp(-0.03) rtol = 1.0e-10
        @test isfinite(rate(zero(sw, 30000.0)))
    end

    @testset "the origin is preserved" begin
        @test discount(Yield.SmithWilson(ufr = 0.03, α = 0.1), 0.0, 1.0) == exp(-0.03)
        @test discount(sw, 0.0, 0.0) === 1.0
        @test discount(Yield.ForwardStarting(sw, 1.0) * 2.0, 0.0) == 1.0
    end

    @testset "intervals from 0 reproduce discount(c, t) bit for bit" begin
        curves = [
            Yield.Constant(0.04), Yield.Constant(Continuous(0.03)), zrc_lin,
            ZeroRateCurve(rates, tenors, Spline.Cubic()), ZeroRateCurve(rates, tenors), ns,
            Yield.NelsonSiegelSvensson(2.5, 3.0, 0.04, -0.02, 0.01, -0.005),
            Yield.CairnsPritchard(0.5, 1.5, 0.04, -0.02, -0.01),
            Yield.CairnsPritchardExtended(0.5, 1.5, 3.0, 0.04, -0.02, -0.01, 0.005),
            sw, zrc_lin + ns, ns - Yield.Constant(Continuous(0.01)), 2 * zrc_lin,
            Yield.TenorShift(zrc_lin, (z, t) -> z + Continuous(0.01)),
            Yield.ProjectedShift(ns, (τ, z, t) -> z + Continuous(0.001 * τ), 2.0),
            Yield.ForwardStarting(zrc_lin, 1.0),
            FinanceModels.ShortRate.Vasicek(0.1, 0.05, 0.01, Continuous(0.03)),
            FinanceModels.ShortRate.HullWhite(0.1, 0.01, zrc_lin),
        ]
        for c in curves, t in (0.0, 0.5, 1.0, 2.5, 7.0, 30.0)
            @test discount(c, 0.0, t) === discount(c, t)
        end
        # a discount-native curve without a closed-form log-discount agrees to rounding
        cir = FinanceModels.ShortRate.CoxIngersollRoss(0.1, 0.05, 0.01, Continuous(0.03))
        @test discount(cir, 0.0, 7.0) ≈ discount(cir, 7.0) rtol = 1.0e-14
    end

    @testset "interval algebra" begin
        a, b = zrc_lin, ns
        for (s, t) in ((0.5, 3.0), (1.0, 7.5), (2.5, 12.0))
            @test discount(a + b, s, t) ≈ discount(a, s, t) * discount(b, s, t) rtol = 1.0e-14
            @test discount(a - b, s, t) ≈ discount(a, s, t) / discount(b, s, t) rtol = 1.0e-14
            @test discount(2 * a, s, t) ≈ discount(a, s, t)^2 rtol = 1.0e-14
            @test discount(Yield.ForwardStarting(b, 1.5), s, t) ≈ discount(b, 1.5 + s, 1.5 + t) rtol = 1.0e-14
            # reversed intervals invert; accumulation is the reverse interval
            @test discount(a, t, s) ≈ inv(discount(a, s, t)) rtol = 1.0e-14
            @test accumulation(a, s, t) == discount(a, t, s)
            # the composition law
            m = (s + t) / 2
            @test discount(b, s, t) ≈ discount(b, s, m) * discount(b, m, t) rtol = 1.0e-14
            # the forward rate is the interval's continuously compounded rate
            @test rate(forward(b, s, t)) ≈ -log(discount(b, s, t)) / (t - s) rtol = 1.0e-13
        end
    end

    @testset "equal endpoints are the identity" begin
        for c in (zrc_lin, ns, sw, ZeroRateCurve(rates, tenors))
            @test discount(c, 3.0, 3.0) === 1.0
            @test accumulation(c, 3.0, 3.0) === 1.0
        end
        # L is infinite at both ends; the empty interval is still 1
        @test discount(zrc_lin, Inf, Inf) === 1.0
        @test discount(zrc_lin, 5.0, Inf) == 0.0
        @test discount(zrc_lin, Inf, 5.0) == Inf
    end

    @testset "curves that define only discount" begin
        c = __HalfFlatCurve()
        @test discount(c, 0.0, 1.0) ≈ exp(-0.03) rtol = 1.0e-14
        @test discount(c, 2.0, 5.0) ≈ exp(-0.09) rtol = 1.0e-14
        @test rate(forward(c, 0.0, 1.0)) ≈ 0.03 rtol = 1.0e-12
        @test discount(__NegativeDiscountCurve(), 1.0, 2.0) ≈ exp(-0.03) rtol = 1.0e-14
        @test_throws DomainError forward(__NegativeDiscountCurve(), 1.0, 2.0)
    end

    @testset "Smith-Wilson intervals are exact for either sign of the discount factor" begin
        # fitted to arbitrary prices, the discount factor turns negative (1 + s < 0) at some times
        u = [1.0, 5.0]
        neg = Yield.SmithWilson(u, [1.0 0.0; 0.0 1.0], [0.9, -0.5]; ufr = 0.03, α = 0.1)
        t_neg = only(t for t in 2.0:0.5:6.0 if discount(neg, t) < 0 && discount(neg, t - 0.5) > 0)
        for (s, t) in ((0.5, t_neg), (t_neg, 7.0), (0.5, 3.0))
            @test discount(neg, s, t) ≈ discount(neg, t) / discount(neg, s) rtol = 1.0e-13
        end
        @test discount(Yield.ForwardStarting(neg, 0.5), t_neg - 0.5) ≈ discount(neg, t_neg) / discount(neg, 0.5) rtol = 1.0e-13
        @test_throws DomainError zero(neg, t_neg)

        # Rebased where both absolute discount factors are negative, the relative factor is positive
        # and has a real log-discount; rebased across the sign change, it has none.
        f = Yield.ForwardStarting(neg, 5.0)
        @test discount(neg, 5.0) < 0 && discount(neg, 6.0) < 0 && discount(neg, 7.0) < 0
        d1, d2 = discount(neg, 6.0) / discount(neg, 5.0), discount(neg, 7.0) / discount(neg, 5.0)
        @test rate(zero(f, 1.0)) ≈ -log(d1) rtol = 1.0e-13
        @test rate(forward(f, 1.0, 2.0)) ≈ -log(d2 / d1) rtol = 1.0e-12
        @test discount(2 * f, 1.0) ≈ d1^2 rtol = 1.0e-13
        @test discount(2 * f, 0.0) == 1.0
        @test discount(f + Yield.Constant(Continuous(0.01)), 1.0) ≈ d1 * exp(-0.01) rtol = 1.0e-13
        @test_throws DomainError zero(Yield.ForwardStarting(neg, 0.5), t_neg - 0.5)
    end

    @testset "limits at infinity combine the components' tails" begin
        # Flat forwards of 4% and -2%: L is Inf and -Inf separately, and the sum tends to +Inf.
        a = ZeroRateCurve([0.04, 0.04], [1.0, 2.0], Spline.Linear())
        b = ZeroRateCurve([-0.02, -0.02], [1.0, 2.0], Spline.Linear())
        for c in (a + b, b + a, a - a / 2)
            @test discount(c, Inf) === 0.0
            @test discount(c, 1.0, Inf) === 0.0
        end
        @test discount(b - a, Inf) == Inf
        # Slopes that cancel leave the intercept. With knot rates 1/32 and 3/64 at 1 and 2 (exact in
        # binary), the flat-forward tail is L(t) = 0.0625·t − 0.03125, so subtracting a flat 6.25%
        # leaves exp(0.03125) at infinity.
        c = ZeroRateCurve([0.03125, 0.046875], [1.0, 2.0], Spline.Linear())
        flat = Yield.Constant(Continuous(0.0625))
        @test discount(c - flat, Inf) == exp(0.03125)
        @test discount(c - flat, 1.0, Inf) == 1.0   # L(1) = 0.03125 − 0.0625 as well
        @test discount(flat - c, Inf) == exp(-0.03125)
        # A rebased curve's tail is shifted: L(t + 3) − L(3) − L(t) → 0.0625·3 − L(3) = 0.03125.
        @test discount(Yield.ForwardStarting(c, 3.0) - c, Inf) ≈ exp(-0.03125) rtol = 1.0e-15
        # `:linear` tails grow as γ·t²; equal slopes γ = 1/4 cancel, leaving (z₁ − z₂)·t = t/8.
        l1 = ZeroRateCurve([0.25, 0.5], [1.0, 2.0], Spline.Linear(); extrapolation = :linear)
        l2 = ZeroRateCurve([0.125, 0.375], [1.0, 2.0], Spline.Linear(); extrapolation = :linear)
        @test discount(l1 - l2, Inf) == 0.0
        @test discount(l2 - l1, Inf) == Inf
        @test discount(l1 - l1, Inf) == 1.0
        # a zero multiple has L = 0 everywhere (0·Inf would be NaN)
        @test discount(0 * a, Inf) == 1.0
        # The combined deciding coefficient is classified like a single tail's: zero with partials
        # has no derivative, and a derivative is zero where the limit is 0.
        net = r -> discount(ZeroRateCurve([r, r], [1.0, 2.0], Spline.Linear()) - flat, Inf)
        @test net(0.0625) == 1.0
        @test_throws "no derivative" ForwardDiff.derivative(net, 0.0625)
        @test ForwardDiff.derivative(net, 0.07) == 0.0
        # A negative time still reaches the components, which reject it.
        @test_throws DomainError discount(a + b, -Inf)
        # Curves without a closed-form tail contribute their zero rate at infinity as the slope, as
        # FinanceModels 6 combined zero rates; an exact cancellation of such slopes is still 0·Inf.
        ns = Yield.NelsonSiegel(1.0, 0.05, -0.02, 0.01)
        @test discount(ns + b, Inf) == 0.0 == exp(-rate(zero(ns + b, Inf)) * Inf)
        ts = Yield.TenorShift(b, (z, t) -> z - Continuous(0.03))
        @test discount(ts + a, Inf) == Inf == exp(-rate(zero(ts + a, Inf)) * Inf)
        @test isnan(discount(ns - ns, Inf))
        # A fallback curve whose zero rate tends to ±Inf has L growing faster than t, at an order it
        # doesn't state. Here u has zero rate t/32 (L = t²/32) and q the `:linear` tail t/64
        # (L = t²/64), both exact in binary, so L(u − q) = t²/64 → +Inf; but the fallback can't
        # tell t²/32 from t^1.5, which t²/64 would outgrow. Against a known quadratic term of the
        # other sign the limit is unsupported: NaN, not the opposite answer.
        u = Yield.TenorShift(Yield.Constant(Continuous(0.0)), (z, t) -> z + Continuous(t / 32))
        u2 = Yield.TenorShift(Yield.Constant(Continuous(0.0)), (z, t) -> z + Continuous(t / 16))
        q = ZeroRateCurve([1 / 64, 2 / 64], [1.0, 2.0], Spline.Linear(); extrapolation = :linear)
        @test discount(u - q, 32.0) == exp(-16.0)   # finite times are exact
        @test isnan(discount(u - q, Inf))
        @test isnan(discount(q - u, Inf))
        @test isnan(discount(u - u2, Inf))   # two unknown growths of opposite sign
        # The unknown growth decides against linear or flat tails and same-sign quadratic ones,
        # and survives scaling and rebasing.
        @test discount(u + Yield.Constant(0.03), Inf) == 0.0
        @test discount(u - a, Inf) == 0.0
        @test discount(u + q, Inf) == 0.0
        @test discount(u + u2, Inf) == 0.0
        @test discount(2 * u, Inf) == 0.0
        @test discount(-1 * u, Inf) == Inf
        @test discount(Yield.ForwardStarting(u, 5.0), Inf) == 0.0
        @test isnan(discount(Yield.ForwardStarting(u - q, 5.0), Inf))
        # With L = t^1.5 + k·t²/64 the limit is 0 at k = 0 but Inf for every k < 0 (NaN here,
        # unsupported), so it has no derivative at k = 0; a zero quadratic term carrying partials
        # throws, as it does for a known tail. Away from 0 the limit is locally constant.
        v = Yield.TenorShift(Yield.Constant(Continuous(0.0)), (z, t) -> z + Continuous(sqrt(t)))
        lim(k) = discount(v + k * q, Inf)
        @test lim(0.0) == 0.0
        @test isnan(lim(-1.0))
        @test_throws "no derivative" ForwardDiff.derivative(lim, 0.0)
        @test ForwardDiff.derivative(lim, 1.0) == 0.0
    end

    @testset "endpoint derivatives" begin
        cd(f, x; h = 1.0e-6) = (f(x + h) - f(x - h)) / (2h)
        for c in (zrc_lin, ns, sw, zrc_lin + ns, Yield.ForwardStarting(ns, 2.0))
            f_from(x) = discount(c, x, 7.0)
            f_to(x) = discount(c, 1.3, x)
            @test ForwardDiff.derivative(f_from, 1.3) ≈ cd(f_from, 1.3) rtol = 1.0e-7
            @test ForwardDiff.derivative(f_to, 7.0) ≈ cd(f_to, 7.0) rtol = 1.0e-7
        end
        # coinciding primal endpoints with different partials: the derivative is not erased. It is
        # minus the instantaneous forward: z(t) + t·z′(t) = 0.035 for the linear zero rates at 2.5,
        # and β₀ + β₁e^{-q} + β₂·q·e^{-q} (q = t/τ) for Nelson-Siegel.
        @test ForwardDiff.derivative(h -> discount(zrc_lin, 2.5, 2.5 + h), 0.0) ≈ -0.035 rtol = 1.0e-12
        @test ForwardDiff.derivative(h -> discount(ns, 4.0, 4.0 + h), 0.0) ≈
            -(0.05 - 0.02 * exp(-4.0) + 0.01 * 4.0 * exp(-4.0)) rtol = 1.0e-12
        mc = ZeroRateCurve(rates, tenors)
        @test ForwardDiff.derivative(h -> discount(mc, 3.0, 3.0 + h), 0.0) ≈
            -Yield.instantaneous_forward(mc, 3.0) rtol = 1.0e-12
    end
end
