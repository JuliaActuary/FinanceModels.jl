@testset "ZeroRateCurve" begin
    rates = [0.02, 0.03, 0.035, 0.04]
    tenors = [1.0, 2.0, 5.0, 10.0]

    @testset "MonotoneConvex (default)" begin
        zrc = ZeroRateCurve(rates, tenors)

        @testset "t=0 returns 1.0" begin
            @test discount(zrc, 0.0) == 1.0
        end

        @testset "exact tenor points" begin
            for (r, t) in zip(rates, tenors)
                @test discount(zrc, t) ≈ exp(-r * t)
            end
        end

        @testset "callable interface" begin
            @test zrc(1.0) ≈ exp(-0.02 * 1.0)
            @test zrc(0.0) == 1.0
        end

        @testset "inherited methods" begin
            # zero rate extraction
            z = zero(zrc, 2.0)
            @test FinanceCore.rate(z) ≈ 0.03 atol = 1e-10
        end
    end

    @testset "Linear" begin
        zrc = ZeroRateCurve(rates, tenors, Spline.Linear())

        @testset "exact tenor points" begin
            for (r, t) in zip(rates, tenors)
                @test discount(zrc, t) ≈ exp(-r * t)
            end
        end

        @testset "interpolation between tenors" begin
            # t=3.5 is between tenors 2.0 (r=0.03) and 5.0 (r=0.035)
            t = 3.5
            w = (3.5 - 2.0) / (5.0 - 2.0)
            r_interp = 0.03 + w * (0.035 - 0.03)
            @test discount(zrc, t) ≈ exp(-r_interp * t) atol = 1e-10
        end
    end

    @testset "from AbstractYieldModel" begin
        @testset "round-trip with Constant" begin
            c = Yield.Constant(0.05)
            tenors = [1.0, 2.0, 5.0, 10.0]
            zrc = ZeroRateCurve(c, tenors)
            for t in tenors
                @test discount(zrc, t) ≈ discount(c, t) atol = 1e-10
            end
        end

        @testset "round-trip with NelsonSiegel" begin
            ns = Yield.NelsonSiegel(1.0, 0.04, -0.02, 0.01)
            tenors = [1.0, 2.0, 5.0, 10.0, 20.0]
            zrc = ZeroRateCurve(ns, tenors)
            for t in tenors
                @test discount(zrc, t) ≈ discount(ns, t) atol = 1e-10
            end
        end

        @testset "explicit spline kwarg" begin
            c = Yield.Constant(0.04)
            tenors = [1.0, 5.0, 10.0]
            zrc = ZeroRateCurve(c, tenors; spline=Spline.Linear())
            for t in tenors
                @test discount(zrc, t) ≈ discount(c, t) atol = 1e-10
            end
        end

        @testset "unsorted tenors are sorted automatically" begin
            c = Yield.Constant(0.04)
            zrc_sorted = ZeroRateCurve(c, [1.0, 5.0, 10.0])
            zrc_unsorted = ZeroRateCurve(c, [10.0, 1.0, 5.0])
            for t in [1.0, 3.0, 5.0, 7.0, 10.0]
                @test discount(zrc_sorted, t) ≈ discount(zrc_unsorted, t) atol = 1e-10
            end
        end

        @testset "error on non-positive tenors" begin
            c = Yield.Constant(0.05)
            @test_throws ArgumentError ZeroRateCurve(c, [0.0, 1.0, 2.0])
            @test_throws ArgumentError ZeroRateCurve(c, [-1.0, 1.0, 2.0])
        end
    end

    @testset "Cubic" begin
        cubic_rates = [0.02, 0.025, 0.03, 0.035, 0.04]
        cubic_tenors = [1.0, 2.0, 5.0, 7.0, 10.0]
        zrc = ZeroRateCurve(cubic_rates, cubic_tenors, Spline.Cubic())

        @testset "exact tenor points" begin
            for (r, t) in zip(cubic_rates, cubic_tenors)
                @test discount(zrc, t) ≈ exp(-r * t) atol = 1e-8
            end
        end

        @testset "two-point case matches linear" begin
            r2 = [0.03, 0.05]
            t2 = [1.0, 5.0]
            zrc_lin = ZeroRateCurve(r2, t2, Spline.Linear())
            zrc_cub = ZeroRateCurve(r2, t2, Spline.Cubic())
            @test discount(zrc_lin, 3.0) ≈ discount(zrc_cub, 3.0) atol = 1e-6
        end
    end

    @testset "eager-build: ForwardDiff pass-through" begin
        # The eager build runs once with the input rate type; for Dual-typed
        # rates, the model itself becomes Dual-typed and propagates through
        # discount. At an exact knot t = tenors[k], discount = exp(-rates[k] * t),
        # so ∂discount/∂rates[k] = -t · discount.
        using ForwardDiff
        tenors_ad = [1.0, 2.0, 5.0, 10.0]
        rates_ad  = [0.02, 0.03, 0.035, 0.04]
        for spl in (Spline.Linear(), Spline.MonotoneConvex())
            f(r) = discount(ZeroRateCurve([r, rates_ad[2:end]...], tenors_ad, spl), 1.0)
            ad = ForwardDiff.derivative(f, 0.02)
            @test ad ≈ -1.0 * exp(-0.02 * 1.0) atol = 1e-12
        end
    end

    @testset "eager-build: structural equality preserved" begin
        # Two ZRCs built from `==`-equal inputs must compare `==` despite the
        # new `_model` field carrying different prebuilt interpolant instances.
        r = [0.02, 0.03, 0.04]; t = [1.0, 2.0, 5.0]
        @test ZeroRateCurve(r, t, Spline.Linear()) == ZeroRateCurve(copy(r), copy(t), Spline.Linear())
        @test hash(ZeroRateCurve(r, t, Spline.Linear())) == hash(ZeroRateCurve(copy(r), copy(t), Spline.Linear()))
        # Different rates must compare unequal
        @test ZeroRateCurve(r, t, Spline.Linear()) != ZeroRateCurve(r .+ 0.01, t, Spline.Linear())
        # Different splines must compare unequal — relies on Sp.SplineCurve
        # subtypes implementing equality correctly (singleton splines under
        # `Sp.Linear()`, `Sp.Cubic()` etc. are `===` to other instances of the
        # same type).
        @test ZeroRateCurve(r, t, Spline.Linear()) != ZeroRateCurve(r, t, Spline.Cubic())
        @test hash(ZeroRateCurve(r, t, Spline.Linear())) != hash(ZeroRateCurve(r, t, Spline.Cubic()))
    end
end
