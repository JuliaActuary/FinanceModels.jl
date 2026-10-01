using ForwardDiff

# A curve that defines only `discount`, with D(0) ≠ 1: intervals keep D(to)/D(from) semantics.
struct __HalfFlatCurve <: Yield.AbstractYieldModel end
FinanceCore.discount(::__HalfFlatCurve, t) = 0.5 * exp(-0.03 * t)

# A curve with negative discount factors (as a Smith-Wilson fit to arbitrary prices can have): its
# interval factors are still the exact ratio, but it has no log-discount at a single time, only over
# an interval whose factor is positive.
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

    @testset "L(0) = 0 for every log-native curve" begin
        # An interval from an exact 0 skips L(from), which relies on this.
        path = RatePath(
            FinanceModels.DataInterpolations.LinearInterpolation(
                [0.0, 0.03, 0.065], [0.0, 1.0, 2.0];
                extrapolation = FinanceModels.DataInterpolations.ExtrapolationType.Extension
            )
        )
        curves = Any[
            Yield.Constant(0.04), Yield.Constant(Continuous(-0.01)),
            [
                ZeroRateCurve(rates, tenors, s; extrapolation = e)
                    for s in (Spline.Linear(), Spline.Cubic(), Spline.PCHIP())
                    for e in (:flat_forward, :flat_zero, :linear, :extension, Yield.FlatForwardAt(Continuous(0.05)))
            ]...,
            [ZeroRateCurve(rates, tenors; extrapolation = e) for e in (:flat_forward, :flat_zero, :linear)]...,
            ZeroRateCurve([-0.01, 0.0, 0.01, 0.02], tenors, Spline.Linear()),
            ns, Yield.NelsonSiegelSvensson(2.5, 3.0, 0.04, -0.02, 0.01, -0.005),
            Yield.CairnsPritchard(0.5, 1.5, 0.04, -0.02, -0.01),
            Yield.CairnsPritchardExtended(0.5, 1.5, 3.0, 0.04, -0.02, -0.01, 0.005),
            Yield.TenorShift(zrc_lin, (z, t) -> z + Continuous(0.01)),
            Yield.ProjectedShift(ns, (τ, z, t) -> z + Continuous(0.001 * τ), 2.0),
            FinanceModels.ShortRate.Vasicek(0.1, 0.05, 0.01, Continuous(0.03)), path,
        ]
        for c in curves
            @test Yield.__log_native(c)
            @test iszero(Yield.__log_discount(c, 0.0))
        end
        # A wrapper takes its components' intervals, so from 0 it is `discount(w, t)` exactly too.
        for w in (zrc_lin + ns, ns - Yield.Constant(Continuous(0.01)), 2 * zrc_lin, 2 * (zrc_lin + ns))
            @test discount(w, 0.0, 7.3) === discount(w, 7.3)
        end
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
        # a wrapper takes an infinite endpoint from its combined tail: the components' limits are
        # Inf and -Inf here, but the composite's forward is 2%
        w = Yield.Constant(Continuous(0.04)) + Yield.Constant(Continuous(-0.02))
        @test discount(w, 5.0, Inf) == 0.0
        @test discount(w, Inf, 5.0) == Inf
        @test discount(w, Inf, Inf) === 1.0
        @test discount(2 * Yield.Constant(Continuous(-0.02)), 5.0, Inf) == Inf
    end

    @testset "curves that define only discount" begin
        c = __HalfFlatCurve()
        @test discount(c, 0.0, 1.0) ≈ exp(-0.03) rtol = 1.0e-14
        @test discount(c, 2.0, 5.0) ≈ exp(-0.09) rtol = 1.0e-14
        @test rate(forward(c, 0.0, 1.0)) ≈ 0.03 rtol = 1.0e-12
        # Composed or scaled, the curve still has L(0) ≠ 0, so an interval from 0 keeps the ratio
        # and is continuous in its start.
        k = Yield.Constant(Continuous(0.01))
        for (w, expected) in (
                (c + k, exp(-0.04)), (k - c, exp(0.02)), (2 * c, exp(-0.06)),
                (zrc_lin + c, discount(zrc_lin, 1.0) * exp(-0.03)), (2 * (c + k), exp(-0.08)),
                (c + sw, exp(-0.03) * discount(sw, 0.0, 1.0)),
            )
            @test discount(w, 0.0, 1.0) ≈ expected rtol = 1.0e-14
            @test discount(w, 0.0, 1.0) ≈ discount(w, 1.0) / discount(w, 0.0) rtol = 1.0e-14
            @test discount(w, 1.0e-12, 1.0) ≈ discount(w, 0.0, 1.0) rtol = 1.0e-10
        end
        @test discount(__NegativeDiscountCurve(), 1.0, 2.0) ≈ exp(-0.03) rtol = 1.0e-14
        @test rate(forward(__NegativeDiscountCurve(), 1.0, 2.0)) ≈ 0.03 rtol = 1.0e-12
        @test_throws DomainError zero(__NegativeDiscountCurve(), 1.0)
    end

    @testset "wrappers keep their components' far-tail intervals" begin
        # Both of each wrapper's discount factors underflow here; the components' own intervals
        # (Smith–Wilson's ratio, a rebased curve's base interval) stay finite.
        k = Yield.Constant(Continuous(0.01))
        sw0 = Yield.SmithWilson(ufr = 0.03, α = 0.1)
        fs = Yield.ForwardStarting(Yield.Constant(Continuous(0.03)), 1.0)
        for (w, f) in (
                (sw0 + k, 0.04), (sw0 - k, 0.02), (2 * sw0, 0.06), (sw + k, 0.04), (2 * sw, 0.06),
                (fs + k, 0.04), (2 * fs, 0.06), (2 * (sw + k), 0.08), (zrc_lin + sw, 0.03 + rate(forward(zrc_lin, 30000.0, 30001.0))),
            )
            @test discount(w, 30000.0, 30001.0) ≈ exp(-f) rtol = 1.0e-12
        end
        # the primal view that `implied_quote` solves on keeps each curve's own interval (it took the
        # ratio of the underflowed factors, NaN, and was one ulp off for Smith–Wilson)
        for w in (zrc_lin + ns, sw0 + k, 2 * sw, fs + k, Yield.ForwardStarting(zrc_lin, 1.0), sw)
            for (a, b) in ((29000.0, 30000.0), (2.0, 7.0))
                @test discount(Yield.__PrimalCurve(w), a, b) === discount(w, a, b)
            end
        end
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
        # a composite interval takes Smith–Wilson's own signed interval
        @test discount(neg + Yield.Constant(Continuous(0.01)), 5.0, 7.0) ≈ discount(neg, 5.0, 7.0) * exp(-0.02) rtol = 1.0e-13
        d1, d2 = discount(neg, 6.0) / discount(neg, 5.0), discount(neg, 7.0) / discount(neg, 5.0)
        @test rate(zero(f, 1.0)) ≈ -log(d1) rtol = 1.0e-13
        @test rate(forward(f, 1.0, 2.0)) ≈ -log(d2 / d1) rtol = 1.0e-12
        # `forward` takes the curve's own interval, so it exists where both factors are negative (the
        # difference of the two log-discounts was a DomainError there), for the curve and in a
        # composite
        @test rate(forward(neg, 5.0, 7.0)) ≈ -log(discount(neg, 7.0) / discount(neg, 5.0)) / 2 rtol = 1.0e-12
        @test rate(forward(neg + Yield.Constant(Continuous(0.01)), 5.0, 7.0)) ≈ rate(forward(neg, 5.0, 7.0)) + 0.01 rtol = 1.0e-12
        @test_throws DomainError forward(neg, 0.5, t_neg)   # across the sign change there is no rate
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
        # A dual start at an exact 0 still carries its derivative: the instantaneous forward at 0
        # (the linear curve's flat short end, 2%) times the factor.
        @test ForwardDiff.derivative(x -> discount(zrc_lin, x, 7.0), 0.0) ≈ 0.02 * discount(zrc_lin, 7.0) rtol = 1.0e-12
        @test ForwardDiff.derivative(x -> discount(Yield.Constant(Continuous(0.03)), x, 7.0), 0.0) ≈
            0.03 * exp(-0.21) rtol = 1.0e-12
    end

    @testset "time derivatives at t = 0" begin
        # A dual time whose value is 0 skips the exact-zero shortcuts, so each curve's L and zero rate
        # must have no removable 0/0 there (MonotoneConvex and Nelson-Siegel(-Svensson) were NaN).
        # dD/dt at 0 is minus the zero rate at 0, which is the instantaneous forward there.
        mc = ZeroRateCurve(rates, tenors)
        # The Hagan-West boundary condition makes the forward flat at 0, unless the positivity
        # collar binds there, as it does on this steep curve.
        mc_steep = ZeroRateCurve([0.01, 0.05, 0.05, 0.05], tenors)
        mc_origin = ZeroRateCurve([0.02, 0.025, 0.03], [0.0, 1.0, 5.0])   # a knot at t = 0
        nss = Yield.NelsonSiegelSvensson(2.5, 3.0, 0.04, -0.02, 0.01, -0.005)
        f0 = Yield.instantaneous_forward(mc, 0.0)
        curves = (
            mc => f0, mc_steep => 0.0, mc_origin => 0.02,
            ns => 0.03, nss => 0.02, mc + ns => f0 + 0.03, 2 * mc => 2 * f0,
            Yield.TenorShift(mc, (z, t) -> z + Continuous(0.01)) => f0 + 0.01,
            Yield.ProjectedShift(mc, (τ, z, t) -> z + Continuous(0.001 * τ), 2.0) => f0 + 0.002,
            Yield.TenorShift(nss, (z, t) -> z + Continuous(0.01)) => 0.03,
        )
        # a second-order one-sided difference: times below 0 throw
        fd(f, h = 1.0e-4) = (-3 * f(0.0) + 4 * f(h) - f(2h)) / (2h)
        for (c, z0) in curves
            @test rate(zero(c, 0.0)) ≈ z0 rtol = 1.0e-12
            @test ForwardDiff.derivative(t -> discount(c, t), 0.0) ≈ -z0 rtol = 1.0e-12
            @test ForwardDiff.derivative(x -> discount(c, x, 7.0), 0.0) ≈ z0 * discount(c, 7.0) rtol = 1.0e-12
            # the second derivative and the zero rate's slope at 0 match finite differences
            d1(t) = ForwardDiff.derivative(u -> discount(c, u), t)
            @test ForwardDiff.derivative(d1, 0.0) ≈ fd(d1) rtol = 1.0e-6
            zr(t) = rate(zero(c, t))
            @test ForwardDiff.derivative(zr, 0.0) ≈ fd(zr) atol = 1.0e-9
        end
        @test ForwardDiff.derivative(t -> rate(zero(mc_steep, t)), 0.0) > 1.0e-3   # a slope is exercised
        # A single knot at t = 0 is a flat curve: a dual time at 0 compares above the knot and enters
        # the tail, whose zero rate α + β·tₙ/t formed 0/0 there.
        for spline in (Spline.MonotoneConvex(), Spline.Linear()),
                ex in (:flat_forward, :flat_zero, :linear, Yield.FlatForwardAt(Continuous(0.03)))
            c1 = ZeroRateCurve([0.03], [0.0], spline; extrapolation = ex)
            @test ForwardDiff.derivative(t -> discount(c1, t), 0.0) ≈ -0.03 rtol = 1.0e-14
            @test ForwardDiff.derivative(t -> rate(zero(c1, t)), 0.0) == 0
            @test ForwardDiff.derivative(s -> ForwardDiff.derivative(t -> discount(c1, t), s), 0.0) ≈ 0.03^2 rtol = 1.0e-14
            @test (discount(c1, 0.0), discount(c1, 2.0), discount(c1, Inf)) == (1.0, exp(-0.06), 0.0)
            @test rate(zero(c1, 0.0)) == rate(zero(c1, 1.0e-8)) == 0.03
            # Another forward would make the zero rate jump at the origin, from the knot's rate to it.
            @test_throws ArgumentError ZeroRateCurve([0.03], [0.0], spline; extrapolation = Yield.FlatForwardAt(Continuous(0.04)))
            # also when the tenor carries a derivative: the check reads its primal value
            @test_throws ArgumentError ZeroRateCurve(
                [0.03], [ForwardDiff.Dual(0.0, 1.0)], spline; extrapolation = Yield.FlatForwardAt(Continuous(0.04))
            )
        end
        # the check compares primal values: a knot rate carrying a partial still builds, and the
        # tail's forward, not the knot's rate, prices every positive time
        at(z) = ZeroRateCurve([z], [0.0]; extrapolation = Yield.FlatForwardAt(Continuous(0.03)))
        @test ForwardDiff.derivative(z -> discount(at(z), 2.0), 0.03) == 0
        # the rate partial survives: ∂²D/∂z∂t at 0 is -1 for the flat curve
        @test ForwardDiff.derivative(z -> ForwardDiff.derivative(t -> discount(ZeroRateCurve([z], [0.0]), t), 0.0), 0.03) ≈ -1 rtol = 1.0e-14
        # A curve without its own `zero` has only L(t)/t, 0/0 at 0: its discount has a derivative
        # there, but a zero-rate transformation of it throws rather than returning NaN.
        @test ForwardDiff.derivative(t -> discount(sw, t), 0.0) ≈ -rate(zero(sw, 1.0e-9)) rtol = 1.0e-6
        @test_throws DomainError ForwardDiff.derivative(t -> discount(Yield.TenorShift(sw, (z, t) -> z), t), 0.0)
        @test_throws DomainError ForwardDiff.derivative(t -> rate(zero(sw, t)), 0.0)
        # An exact 0 is the zero rate's limit there, the instantaneous forward (it was the 0/0 NaN)
        @test rate(zero(sw, 0.0)) == Yield.instantaneous_forward(sw, 0.0)
        @test rate(zero(sw, 0.0)) ≈ rate(zero(sw, 1.0e-7)) rtol = 1.0e-6
        # Float32 parameters keep Float32 zero rates, at a dual 0 (the decay's series) and elsewhere
        ns32 = Yield.NelsonSiegel(1.0f0, 0.05f0, -0.02f0, 0.01f0)
        nss32 = Yield.NelsonSiegelSvensson(2.5f0, 3.0f0, 0.04f0, -0.02f0, 0.01f0, -0.005f0)
        for c in (ns32, nss32), t in (0.0f0, 1.0f0)
            @test rate(zero(c, t)) isa Float32
            @test ForwardDiff.derivative(u -> rate(zero(c, u)), t) isa Float32
        end
    end

    @testset "instantaneous forwards, and zero rates at t = 0" begin
        quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
        vas = ShortRate.Vasicek(0.1, 0.03, 0.01, 0.04)
        cir = quiet(() -> ShortRate.CoxIngersollRoss(0.1, 0.03, 0.05, 0.04))
        cir_explosive = quiet(() -> ShortRate.CoxIngersollRoss(-0.1, 0.03, 0.05, 0.04))
        cubic = ZeroRateCurve([0.02, 0.025, 0.03, 0.031], [1.0, 2.0, 5.0, 10.0], Spline.Cubic())
        hw = ShortRate.HullWhite(0.1, 0.01, cubic)
        nss = Yield.NelsonSiegelSvensson(1.5, 3.0, 0.03, -0.01, 0.005, 0.002)
        curves = (
            Yield.Constant(Continuous(0.04)), ZeroRateCurve(rates, tenors, Spline.Linear()), cubic,
            ZeroRateCurve(rates, tenors), ns, nss, vas, cir, cir_explosive, hw, sw,
            Yield.ForwardStarting(cubic, 1.5), cubic + ns, 0.7 * cubic, vas - Yield.Constant(0.01),
        )
        # The closed forms equal the derivative of each curve's L (the generic method, which
        # Smith–Wilson uses), including at t = 0 and past the last knot
        L′(c, t) = ForwardDiff.derivative(s -> Yield.__log_discount(c, s), t)
        for c in curves, t in (0.0, 0.3, 1.0, 1.7, 4.0, 7.5, 12.0, 40.0)
            @test Yield.instantaneous_forward(c, t) ≈ L′(c, t) atol = 1.0e-14
        end
        # long-run forwards
        @test Yield.instantaneous_forward(vas, Inf) ≈ 0.03 - 0.01^2 / (2 * 0.1^2) rtol = 1.0e-14
        @test Yield.instantaneous_forward(cir, Inf) ≈ 2 * 0.1 * 0.03 / (sqrt(0.1^2 + 2 * 0.05^2) + 0.1) rtol = 1.0e-14
        @test Yield.instantaneous_forward(ns, Inf) == ns.β₀
        # CIR without mean reversion: the forward decays to 0 at τ = ∞; and a σ whose square underflows
        @test Yield.instantaneous_forward(ShortRate.CoxIngersollRoss(0.0, 0.03, 0.1, 0.04), Inf) == 0
        @test Yield.instantaneous_forward(ShortRate.CoxIngersollRoss(0.0, 0.03, 1.0e-200, 0.04), 10.0) ≈ 0.04 rtol = 1.0e-15
        # a < 0 with a tiny σ: 2γ/D′ is ~1e200 and e^{-γτ} ~1e-218, so B′ overflowed when squared first
        # (a 2048-bit reference gives 1.1399322250785743e178)
        @test Yield.instantaneous_forward(ShortRate.CoxIngersollRoss(-0.1, 0.0, 1.0e-100, 0.04), 5000.0) ≈ 1.1399322250785743e178 rtol = 1.0e-12
        # and through a = σ = 0: f = r(1 - σ²B²/2) at a = 0, so ∂²f/∂σ² = -rτ²
        fσ(σ) = Yield.instantaneous_forward(ShortRate.CoxIngersollRoss(0.0, 0.03, σ, 0.04), 10.0)
        @test ForwardDiff.derivative(s -> ForwardDiff.derivative(fσ, s), 0.0) ≈ -0.04 * 10.0^2 rtol = 1.0e-13

        # A curve without its own `zero` takes the limit at an exact 0, its short rate; it was NaN
        for c in (vas, cir, cir_explosive, hw, sw, Yield.ForwardStarting(cubic, 1.5), vas - Yield.Constant(0.01))
            z0 = rate(zero(c, 0.0))
            @test z0 == Yield.instantaneous_forward(c, 0.0)
            @test z0 ≈ rate(zero(c, 1.0e-7)) rtol = 1.0e-6
        end
        @test rate(zero(vas, 0.0)) == 0.04
        @test rate(zero(cir, 0.0)) ≈ 0.04 rtol = 1.0e-14
        @test rate(zero(hw, 0.0)) == 0.02   # its curve's flat short end
        # A shift's forward at 0 is its rule applied to the base's limit, so an identity shift of a curve
        # without its own zero rate keeps that curve's forward, at 0 too (it threw), and a shift
        # parameter keeps its derivative there
        for c in (vas, cir, sw), t in (0.0, 0.5, 7.5)
            @test Yield.instantaneous_forward(Yield.TenorShift(c, (z, t) -> z), t) ≈
                Yield.instantaneous_forward(c, t) rtol = 1.0e-12
        end
        shifted(h) = Yield.TenorShift(vas, (z, t) -> z + Continuous(h))
        @test ForwardDiff.derivative(h -> Yield.instantaneous_forward(shifted(h), 0.0), 0.01) == 1
        @test ForwardDiff.derivative(r -> rate(zero(ShortRate.Vasicek(0.1, 0.03, 0.01, r), 0.0)), 0.04) == 1

        # Nelson–Siegel near 0: at t = 1e-16, e^{-q} rounded to 1 and the zero rate jumped from 2% to
        # 2.5%. The short end now agrees with the closed form in high precision on both sides of
        # the series band, |q| < 0.1
        ns_reference(c, t) = setprecision(BigFloat, 512) do
            q = big(t) / big(c.τ₁)
            d = -expm1(-q) / q
            Float64(big(c.β₀) + big(c.β₁) * d + big(c.β₂) * (d - exp(-q)))
        end
        ns2 = Yield.NelsonSiegel(2.0, 0.03, -0.01, 0.005)
        for t in (1.0e-300, 1.0e-16, 1.0e-10, 1.0e-6, 1.0e-3, 0.1999, 0.2001, 1.0, 10.0)
            @test rate(zero(ns2, t)) ≈ ns_reference(ns2, t) rtol = 2.0e-15
        end
        @test rate(zero(ns2, 1.0e-16)) == rate(zero(ns2, 0.0))
        @test rate(zero(nss, 1.0e-16)) == rate(zero(nss, 0.0))
        # and the zero rate's slope there is the series', z′(0) = (f′(0))/2
        f′0 = ForwardDiff.derivative(t -> Yield.instantaneous_forward(ns2, t), 0.0)
        @test ForwardDiff.derivative(t -> rate(zero(ns2, t)), 0.0) ≈ f′0 / 2 rtol = 1.0e-14
        @test ForwardDiff.derivative(t -> rate(zero(ns2, t)), 1.0e-12) ≈ f′0 / 2 rtol = 1.0e-9
    end
end
