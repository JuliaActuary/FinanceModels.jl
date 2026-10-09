# A coupon, margin, basis spread or interest-rate strike may be a `Periodic` rate of the contract's
# frequency, whose nominal rate is used. Another frequency throws, and a `Continuous` rate has no
# method.
@testset "typed coupons, margins, spreads and strikes" begin
    same_cashflows(a, b) = length(a) == length(b) &&
        all(x.time == y.time && abs(x.amount - y.amount) <= 1.0e-15 for (x, y) in zip(a, b))
    curve = ZeroRateCurve([0.02, 0.03, 0.035, 0.04], [1.0, 2.0, 5.0, 10.0])

    @testset "Bond.Fixed" begin
        # regular schedules, short first stubs, a sub-period bond, a monthly and a negative coupon
        for (c, f, T) in ((0.05, 2, 10), (0.03625, 2, 7.25), (0.045, 4, 0.2), (0.01, 12, 3), (-0.002, 1, 3.5))
            typed = Bond.Fixed(Periodic(c, f), Periodic(f), T)
            numeric = Bond.Fixed(c, Periodic(f), T)
            @test typed isa typeof(numeric)
            @test typed.coupon_rate ≈ c rtol = 1.0e-15
            @test same_cashflows(collect(typed), collect(numeric))
            @test pv(curve, typed) ≈ pv(curve, numeric) rtol = 1.0e-14
        end
        # numbers keep their exact value; a typed coupon can differ in the last bits
        @test Bond.Fixed(0.03625, Periodic(2), 7).coupon_rate === 0.03625
        @test Bond.Fixed(Periodic(0.03625, 2), Periodic(2), 7).coupon_rate == rate(Periodic(0.03625, 2))
        @test_throws ArgumentError Bond.Fixed(Periodic(0.05, 1), Periodic(2), 10)
        @test_throws "frequency 1, but the contract's frequency is 2" Bond.Fixed(Periodic(0.05, 1), Periodic(2), 10)
        @test_throws MethodError Bond.Fixed(Continuous(0.05), Periodic(2), 10)
        # a coupon is a nominal rate, while a par yield on a stub schedule solves for its coupon
        @test Bond.Fixed(Periodic(0.04, 2), Periodic(2), 2.25).coupon_rate ≈ 0.04 rtol = 1.0e-15
        @test !isapprox(ParYield(Periodic(0.04, 2), 2.25).instrument.coupon_rate, 0.04; rtol = 1.0e-6)
        @test_throws ArgumentError ParYield(Periodic(0.04, 1), 10; frequency = 2)
    end

    @testset "Bond.Floating" begin
        projected(b) = collect(Projection(b, Dict(:R => curve), CashflowProjection()))
        # regular, short-stub and sub-period schedules, positive and negative margins
        for (m, f, T) in ((0.002, 4, 5.0), (-0.0015, 4, 2.3), (0.001, 2, 0.25))
            typed = Bond.Floating(Periodic(m, f), Periodic(f), T, :R)
            numeric = Bond.Floating(m, Periodic(f), T, :R)
            @test typed isa typeof(numeric)
            @test same_cashflows(projected(typed), projected(numeric))
            @test pv(Models(curve; index = curve), typed) ≈ pv(Models(curve; index = curve), numeric) rtol = 1.0e-14
        end
        @test_throws ArgumentError Bond.Floating(Periodic(0.002, 2), Periodic(4), 5.0, :R)
        @test_throws MethodError Bond.Floating(Continuous(0.002), Periodic(4), 5.0, :R)
    end

    @testset "FX.ParBasisSwap" begin
        eurusd = FX.Pair(:EUR, :USD)
        estr = fit(Spline.Linear(), ZCBYield.([0.028, 0.03, 0.031, 0.032], [1.0, 2.0, 5.0, 10.0]), Fit.Bootstrap())
        for (b, T) in ((-0.0015, 3.0), (-0.0015, 2.3), (0.001, 0.2))
            typed = FX.ParBasisSwap(eurusd, Periodic(b, 4), T; reference = estr)
            numeric = FX.ParBasisSwap(eurusd, b, T; reference = estr)
            @test same_cashflows(typed.instrument.cashflows, numeric.instrument.cashflows)
        end
        # basis-swap calibration is unchanged
        usd = Yield.Constant(Continuous(0.045))
        tenors = [1.0, 2.0, 5.0]
        quotes(spreads) = [FX.ParBasisSwap(eurusd, s, T; reference = estr) for (s, T) in zip(spreads, tenors)]
        spreads = [-0.001, -0.0015, -0.002]
        numeric = fit(FX.Forwards(eurusd, 1.1, usd, Spline.Linear()), quotes(spreads), Fit.Bootstrap())
        typed = fit(FX.Forwards(eurusd, 1.1, usd, Spline.Linear()), quotes(Periodic.(spreads, 4)), Fit.Bootstrap())
        @test rate.(knot_rates(typed.foreign)) ≈ rate.(knot_rates(numeric.foreign)) rtol = 1.0e-12
        @test_throws ArgumentError FX.ParBasisSwap(eurusd, Periodic(-0.0015, 1), 3.0; reference = estr)
        @test_throws MethodError FX.ParBasisSwap(eurusd, Continuous(-0.0015), 3.0; reference = estr)
    end

    @testset "Option.Cap, Option.Floor and Option.Swaption" begin
        hw = ShortRate.HullWhite(0.1, 0.01, curve)
        # the frequency is an integer or a `Periodic`
        for freq in (4, Periodic(4)), K in (0.03, 0.045)
            @test pv(hw, Option.Cap(Periodic(K, 4), freq, 5.0)) ≈ pv(hw, Option.Cap(K, freq, 5.0)) rtol = 1.0e-12
            @test pv(hw, Option.Floor(Periodic(K, 4), freq, 5.0)) ≈ pv(hw, Option.Floor(K, freq, 5.0)) rtol = 1.0e-12
            @test_throws ArgumentError Option.Cap(Periodic(K, 2), freq, 5.0)
            @test_throws ArgumentError Option.Floor(Periodic(K, 2), freq, 5.0)
        end
        for (freq, n, payer) in ((1, 1, true), (Periodic(2), 2, false))
            typed = Option.Swaption(1.0, 6.0, Periodic(0.035, n), freq; payer)
            numeric = Option.Swaption(1.0, 6.0, 0.035, freq; payer)
            @test typed isa typeof(numeric) && typed.payer == payer
            @test pv(hw, typed) ≈ pv(hw, numeric) rtol = 1.0e-12
            @test_throws ArgumentError Option.Swaption(1.0, 6.0, Periodic(0.035, 12), freq; payer)
        end
        @test_throws MethodError Option.Cap(Continuous(0.03), 4, 5.0)
        @test_throws MethodError Option.Floor(Continuous(0.03), 4, 5.0)
        @test_throws MethodError Option.Swaption(1.0, 6.0, Continuous(0.035), 1)
    end
end
