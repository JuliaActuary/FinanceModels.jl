@testset "Yield.Spline: shared knot grid (owned copies, validation)" begin
    rates = [0.02, 0.03, 0.035, 0.04]
    tenors = [1.0, 2.0, 5.0, 10.0]
    descriptors = (Spline.Linear(), Spline.Quadratic(), Spline.Cubic(), Spline.BSpline(3), Spline.PCHIP(), Spline.Akima())

    @testset "caller-input isolation: $d" for d in descriptors
        r = copy(rates); t = copy(tenors)
        c = Yield.Spline(d, t, r)
        d3 = discount(c, 3.0); z7 = zero(c, 7.0)
        r[1] = 0.2; t[1] = 0.5; r .= 0.0
        @test discount(c, 3.0) == d3
        @test zero(c, 7.0) == z7
        @test discount(c, 3.0) == discount(Yield.Spline(d, tenors, rates), 3.0)
    end

    @testset "validation: $d" for d in descriptors
        @test_throws ArgumentError Yield.Spline(d, [2.0, 1.0, 3.0], [0.1, 0.2, 0.3])        # unsorted
        @test_throws ArgumentError Yield.Spline(d, [1.0, 1.0, 2.0], [0.1, 0.2, 0.3])        # duplicate
        @test_throws ArgumentError Yield.Spline(d, [-1.0, 1.0, 2.0], [0.1, 0.2, 0.3])       # negative
        @test_throws ArgumentError Yield.Spline(d, [1.0, 2.0, 3.0], [0.1, NaN, 0.3])        # non-finite rate
        @test_throws ArgumentError Yield.Spline(d, [1.0, 2.0, Inf], [0.1, 0.2, 0.3])        # non-finite tenor
        @test_throws ArgumentError Yield.Spline(d, [1.0, 2.0], [0.1, 0.2, 0.3])             # length
        @test_throws ArgumentError Yield.Spline(d, Float64[], Float64[])                     # empty
        k = FinanceModels.Yield.__min_knots(d)
        k > 1 && @test_throws ArgumentError Yield.Spline(d, collect(1.0:(k - 1)), fill(0.02, k - 1))
        ck = Yield.Spline(d, collect(1.0:k), fill(0.02, k))
        @test 0 < discount(ck, 0.5 + k / 2) <= 1
        # same errors through `ZeroRateCurve` and `reconstruct`
        @test_throws ArgumentError ZeroRateCurve([0.1, 0.2, 0.3], [2.0, 1.0, 3.0], d)
        @test_throws ArgumentError reconstruct(ck; tenors = [2.0, 1.0, 3.0], rates = [0.1, 0.2, 0.3])
    end

    @testset "promotion: $d" for d in descriptors
        c = Yield.Spline(d, 1:4, [2, 3, 4, 5] ./ 100)           # range tenors, Int-derived rates
        @test discount(c, 2.5) ≈ discount(Yield.Spline(d, [1.0, 2.0, 3.0, 4.0], [0.02, 0.03, 0.04, 0.05]), 2.5)
        ct = Yield.Spline(d, (1, 2, 3, 4), (0.02, 0.03, 0.04, 0.05))  # tuples
        @test discount(ct, 2.5) ≈ discount(c, 2.5)
        @test_throws MethodError Yield.Spline(d, [1.0, 2.0, 3.0], ["a", "b", "c"])   # non-numeric: float(String)
    end

    @testset "fit paths: grid validated up front, results validated" begin
        t = [1.0, 2.0, 5.0, 10.0]
        target = [0.02, 0.025, 0.03, 0.032]
        qs = ZCBYield.(Continuous.(target), t)
        for d in (Spline.Linear(), Spline.Cubic())
            fl = fit(d, qs, Fit.Loss(x -> x^2))
            @test maximum(abs, present_value(fl, q.instrument) - q.price for q in qs) < 1.0e-8
        end
        fb = fit(Spline.Linear(), qs, Fit.Bootstrap())
        @test maximum(abs, present_value(fb, q.instrument) - q.price for q in qs) < 1.0e-8
        @test_throws ArgumentError fit(Spline.Cubic(), qs, Fit.Bootstrap())
        dup = [qs[1], qs[1], qs[3]]
        @test_throws ArgumentError fit(Spline.Linear(), dup, Fit.Loss(x -> x^2))
        @test_throws ArgumentError fit(Spline.Linear(), dup, Fit.Bootstrap())
        # too few quotes for the interpolant is reported before any optimisation runs
        @test_throws ArgumentError fit(Spline.PCHIP(), qs[1:2], Fit.Loss(x -> x^2))

        selected = fit(Spline.Linear(), qs; extrapolation = :flat_zero)
        @test rate(zero(selected, 100.0)) ≈ rate(zero(selected, last(t))) atol = 1.0e-10
        selected_bootstrap = fit(
            Spline.Linear(), qs, Fit.Bootstrap();
            extrapolation = :extension
        )
        @test rate(zero(selected_bootstrap, 100.0)) !=
            rate(zero(selected_bootstrap, last(t)))

        selected_mc = fit(Spline.MonotoneConvex(), qs; extrapolation = :flat_zero)
        @test selected_mc isa Yield.MonotoneConvex
        @test selected_mc.extrapolation === :flat_zero
        @test rate(zero(selected_mc, 100.0)) ≈ last(selected_mc.rates) atol = 1.0e-10
        @test_throws ArgumentError fit(
            Spline.MonotoneConvex(), qs;
            extrapolation = :extension
        )
    end

    @testset "KnotGrid (internal): owned copies, independent promotion, unchecked trial form" begin
        KG = FinanceModels.Yield.KnotGrid
        r = [0.02, 0.03]; t = [1, 2]
        g = KG(r, t, Spline.Linear())
        @test g.rates !== r && g.rates == r && g.tenors == [1.0, 2.0]
        @test eltype(g.tenors) == Float64 && eltype(g.rates) == Float64
        gb = KG((0.02f0, 0.03f0), (big"1", big"2"), Spline.Linear())
        @test eltype(gb.rates) == Float32 && eltype(gb.tenors) == BigFloat   # promoted independently
        @test_throws ArgumentError KG([0.02, 0.03], [2.0, 1.0], Spline.Linear())
        @test_throws ArgumentError KG([0.02, 0.03], [1.0, 2.0], Spline.PCHIP())   # min knots
        # without a spline the grid validates and copies, too
        g2 = KG(r, t)
        @test g2.rates !== r && g2.rates == r && g2.tenors == [1.0, 2.0]
        @test_throws ArgumentError KG([NaN, 0.0], [1.0, 2.0])
        @test_throws ArgumentError KG([0.02, 0.03], [2.0, 1.0])
        # untyped inputs promote from their values; empty ones reach the grid's own length check
        # (on every Julia version, rather than a native empty-reduction error)
        ga = KG(Any[0.02, 0.03f0], Any[1, 2.0])
        @test eltype(ga.rates) == Float64 && eltype(ga.tenors) == Float64
        @test_throws "at least one knot is required (got 0)" KG([], [])
        @test_throws "requires at least 1 knots (got 0)" KG(Any[], (), Spline.Linear())
        # only the internal unchecked form skips copying and validation (optimizer trial curves)
        raw = KG(FinanceModels.Yield.Unchecked(), [NaN, 0.0], [2.0, 1.0])
        @test isnan(raw.rates[1]) && raw.tenors == [2.0, 1.0]
    end

    @testset "curves built from a KnotGrid cannot alias caller vectors" begin
        KG = FinanceModels.Yield.KnotGrid
        # Regression: the two-argument grid used to keep the caller's vector without
        # validation, so mutating it changed `mc.rates` but not the cached forwards.
        rr = [0.02, 0.03, 0.035]
        mc = Yield.MonotoneConvex(KG(rr, [1.0, 2, 5]))
        d3 = discount(mc, 3.0)
        rr[2] = 0.08
        @test mc.rates[2] == 0.03
        @test discount(mc, 3.0) == d3
        @test discount(mc, 3.0) == discount(Yield.MonotoneConvex([0.02, 0.03, 0.035], [1.0, 2.0, 5.0]), 3.0)
        rr = [0.02, 0.03, 0.035]
        sp = FinanceModels.Yield.__build(Spline.Linear(), KG(rr, [1.0, 2, 5]))
        d3 = discount(sp, 3.0)
        rr[2] = 0.08
        @test discount(sp, 3.0) == d3
        # a curve built over a grid checks its own minimum knot count
        @test_throws ArgumentError FinanceModels.Yield.__build(Spline.PCHIP(), KG([0.02, 0.03], [1.0, 2.0]))
        # one knot is a flat curve
        @test rate(zero(FinanceModels.Yield.__build(Spline.Linear(), KG([0.02], [1.0])), 5.0)) == 0.02
    end

    @testset "flat-forward long-end extrapolation" begin
        # The default tail holds the last discrete forward. It neither continues the final
        # polynomial piece (DataInterpolations.Extension made quadratic/cubic curves explode
        # at 100y) nor depends on the interpolant's endpoint derivative.
        long_rates = [0.02, 0.025, 0.03, 0.035, 0.04, 0.042]
        long_tenors = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0]
        tₘ, tₙ = long_tenors[(end - 1):end]
        zₘ, zₙ = long_rates[(end - 1):end]
        f = (zₙ * tₙ - zₘ * tₘ) / (tₙ - tₘ)
        @test f ≈ 0.046

        for d in descriptors
            # Both public knot-curve construction paths use the same policy.
            for curve in (
                    Yield.Spline(d, long_tenors, long_rates),
                    ZeroRateCurve(long_rates, long_tenors, d),
                )
                for t in (tₙ + 1.0e-6, 31.0, 45.0, 100.0, 1.0e4)
                    @test rate(zero(curve, t)) ≈ f + (zₙ - f) * tₙ / t atol = 1.0e-14
                end
                @test rate(zero(curve, Inf)) ≈ f atol = 1.0e-14
                # Every later interval has the anchor as its continuously compounded forward.
                @test rate(forward(curve, tₙ, 2tₙ)) ≈ f atol = 1.0e-14
                @test rate(forward(curve, 2tₙ, 100.0)) ≈ f atol = 1.0e-14
                # Discount factors are continuous at the last knot; the forward may jump there.
                @test discount(curve, tₙ) ≈ exp(-zₙ * tₙ) rtol = 1.0e-14
                for ε in (1.0e-6, 1.0e-9)
                    @test discount(curve, tₙ + ε) ≈ discount(curve, tₙ) rtol = 2ε
                    @test discount(curve, tₙ - ε) ≈ discount(curve, tₙ) rtol = 2ε
                end
            end
        end
    end

    @testset "flat-forward tail is anchored on the last discrete forward" begin
        # An upward-sloping 2% → 4% curve. The spline's endpoint instantaneous forward
        # zₙ + tₙz′(tₙ⁻) is about 2.95% for Quadratic, 3.58% for Cubic and −6.81% for
        # BSpline(3); the tail must instead use the last discrete forward, 4.25%.
        rates = [0.02, 0.025, 0.03, 0.035, 0.04]
        tenors = [1.0, 2.0, 5.0, 10.0, 30.0]
        f = (0.04 * 30 - 0.035 * 10) / (30 - 10)
        @test f ≈ 0.0425
        DI = FinanceModels.DataInterpolations
        E = DI.ExtrapolationType.Extension
        endpoint = Dict(
            Spline.Quadratic() => DI.QuadraticSpline(rates, tenors; extrapolation = E),
            Spline.Cubic() => DI.CubicSpline(rates, tenors; extrapolation = E),
            Spline.BSpline(3) => DI.BSplineInterpolation(rates, tenors, 3, :Average; extrapolation = E),
        )
        for (d, itp) in endpoint
            c = Yield.Spline(d, tenors, rates)
            @test rate(forward(c, 30.0, 31.0)) ≈ f atol = 1.0e-14
            @test rate(forward(c, 30.0, 100.0)) ≈ f atol = 1.0e-14
            @test rate(zero(c, 100.0)) ≈ f + (0.04 - f) * 30 / 100 atol = 1.0e-14
            @test all(t -> rate(forward(c, t, t + 1)) > 0, 30.0:10.0:500.0)
            # the endpoint-derivative anchor the tail no longer uses
            old = 0.04 + 30 * DI.derivative(itp, 30.0)
            @test abs(old - f) > 0.005
        end
        @test 0.04 + 30 * DI.derivative(endpoint[Spline.BSpline(3)], 30.0) < 0
        # a single knot implies no forward other than its own zero rate
        @test Yield.__last_discrete_forward([0.03], [2.0]) == 0.03
        # MonotoneConvex keeps its native boundary instantaneous forward, so its forward
        # curve is continuous at the last knot.
        mc = Yield.MonotoneConvex(rates, tenors)
        fₙ = Yield.instantaneous_forward(mc, 30.0)
        @test fₙ == last(mc._f)
        @test rate(forward(mc, 30.0, 100.0)) ≈ fₙ atol = 1.0e-14
        @test abs(fₙ - f) > 1.0e-4
    end

    @testset "selectable long-end extrapolation" begin
        long_rates = [0.02, 0.025, 0.03, 0.035, 0.04, 0.042]
        long_tenors = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0]
        tₙ, horizon = last(long_tenors), 100.0
        d = Spline.Cubic()

        default = Yield.Spline(d, long_tenors, long_rates)
        explicit_default = Yield.Spline(
            d, long_tenors, long_rates;
            extrapolation = :flat_forward
        )
        @test rate(zero(default, horizon)) == rate(zero(explicit_default, horizon))

        flat_zero = Yield.Spline(d, long_tenors, long_rates; extrapolation = :flat_zero)
        @test rate(zero(flat_zero, horizon)) == last(long_rates)

        linear = Yield.Spline(d, long_tenors, long_rates; extrapolation = :linear)
        z31 = rate(zero(linear, tₙ + 1))
        @test rate(zero(linear, horizon)) ≈
            last(long_rates) + (horizon - tₙ) * (z31 - last(long_rates)) atol = 1.0e-14

        extension = Yield.Spline(d, long_tenors, long_rates; extrapolation = :extension)
        DI = FinanceModels.DataInterpolations
        legacy = DI.CubicSpline(
            long_rates, long_tenors;
            extrapolation = DI.ExtrapolationType.Extension, cache_parameters = true
        )
        @test rate(zero(extension, horizon)) ≈ legacy(horizon) atol = 1.0e-14

        # ZeroRateCurve stores the selection so Accessors/fit reconstruction can preserve it.
        for method in (:flat_forward, :flat_zero, :linear, :extension)
            zrc = ZeroRateCurve(long_rates, long_tenors, d; extrapolation = method)
            @test zrc.extrapolation === method
            @test rate(zero(zrc, horizon)) ≈
                rate(
                zero(
                    Yield.Spline(
                        d, long_tenors, long_rates;
                        extrapolation = method
                    ), horizon
                )
            ) atol = 1.0e-14
        end

        # MonotoneConvex supports the financial boundary policies, but has no
        # DataInterpolations polynomial for `:extension` to continue.
        mc_flat = ZeroRateCurve(long_rates, long_tenors; extrapolation = :flat_zero)
        @test rate(zero(mc_flat, horizon)) ≈ last(long_rates) atol = 1.0e-14
        @test ZeroRateCurve(long_rates, long_tenors; extrapolation = :linear).extrapolation === :linear

        @test_throws ArgumentError Yield.Spline(
            d, long_tenors, long_rates;
            extrapolation = :unknown
        )
    end
end
