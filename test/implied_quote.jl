using ForwardDiff

@testset "implied_quote" begin
    flat = Yield.Constant(0.04)

    @testset "flat-curve closed forms" begin
        @test implied_quote(flat, ZCBYield, 5.0) ≈ 0.04 atol = 1.0e-14
        @test implied_quote(flat, ZCBPrice, 2.0) ≈ 1.04^-2 atol = 1.0e-14
        @test implied_quote(flat, (r, t) -> ParYield(r, t; frequency = 1), 5.0) ≈ 0.04 atol = 1.0e-14
        # Semiannual par yield of an annual-effective 4% curve.
        @test implied_quote(flat, CMTYield, 5.0) ≈ 2 * (sqrt(1.04) - 1) atol = 1.0e-14
        @test implied_quote(flat, OISYield, 0.5) ≈ 0.04 atol = 1.0e-14
        continuous = Yield.Constant(Continuous(0.03))
        @test implied_quote(continuous, ZCBYield, 2.0) ≈ exp(0.03) - 1 atol = 1.0e-14
    end

    @testset "round trip on fitted curves" begin
        tenors = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
        rates = [0.028, 0.03, 0.032, 0.035, 0.037, 0.04]
        for family in (ZCBYield, CMTYield, OISYield, (r, t) -> ParYield(r, t; frequency = 1))
            curve = fit(Spline.Linear(), family.(rates, tenors), Fit.Bootstrap())
            @test implied_quote.(curve, family, tenors) ≈ rates atol = 1.0e-12
        end
    end

    @testset "agrees with par" begin
        curve = fit(Spline.Linear(), CMTYield.([0.03, 0.035, 0.04], [1.0, 3.0, 7.0]), Fit.Bootstrap())
        for t in (2.0, 2.3, 5.0), f in (1, 2, 4)
            @test implied_quote(curve, (r, T) -> ParYield(r, T; frequency = f), t) ≈ FinanceCore.rate(par(curve, t; frequency = f)) atol = 1.0e-9
        end
    end

    @testset "exact first-order derivatives through the curve" begin
        # ZCBYield on a continuous curve: r = exp(x) - 1, dr/dx = exp(x).
        dr = ForwardDiff.derivative(x -> implied_quote(Yield.Constant(Continuous(x)), ZCBYield, 2.0), 0.03)
        @test dr ≈ exp(0.03) rtol = 1.0e-14
        # Against central differences of a parallel continuous-zero bump.
        curve = fit(Spline.Linear(), CMTYield.([0.03, 0.032, 0.035, 0.037], [1.0, 2.0, 3.0, 5.0]), Fit.Bootstrap())
        bumped(s) = Yield.TenorShift(curve, (z, t) -> Continuous(s) + z)
        for t in (1.0, 2.0, 5.0)
            ad = ForwardDiff.derivative(s -> implied_quote(bumped(s), CMTYield, t), 0.0)
            h = 1.0e-6
            fd = (implied_quote(bumped(h), CMTYield, t) - implied_quote(bumped(-h), CMTYield, t)) / (2h)
            @test ad ≈ fd rtol = 1.0e-8
        end
        # The value is the primal root exactly.
        d = implied_quote(bumped(ForwardDiff.Dual(0.0, 1.0)), CMTYield, 5.0)
        @test ForwardDiff.value(d) == implied_quote(curve, CMTYield, 5.0)
        # Chunked gradients and an outer tag created before the call.
        g = ForwardDiff.gradient(s -> implied_quote(bumped(s[1] + s[2]), CMTYield, 5.0), [0.0, 0.0])
        @test g[1] ≈ g[2]
    end

    @testset "the derivative does not depend on the notional" begin
        curve(z) = Yield.Constant(Continuous(z))
        unit = ForwardDiff.derivative(z -> implied_quote(curve(z), ZCBPrice, 2.0), 0.03)
        @test unit ≈ -2 * exp(-0.06) rtol = 1.0e-12
        for n in (1.0e-10, 1.0, 1.0e10)
            family(q, t) = Quote(n * q, Cashflow(n, t))
            @test ForwardDiff.derivative(z -> implied_quote(curve(z), family, 2.0), 0.03) ≈ unit rtol = 1.0e-10
        end
    end

    @testset "nonlinear quotes do not depend on the notional" begin
        # An annual-effective zero-coupon yield and a deposit-style quote, both scaled by a
        # notional n: the implied rate is expm1(z) and its derivative exp(z) at every n. The
        # solve used to stop on an absolute residual, so at n = 1e-10 the rate was about 1e-5
        # off and at n = 1e-14 it was the starting guess.
        curve(z) = Yield.Constant(Continuous(z))
        families = (
            n -> ((q, t) -> Quote(n * (1 + q)^-t, Cashflow(n, t))),
            n -> ((q, t) -> Quote(n, Cashflow(n * (1 + q)^t, t))),
        )
        for make in families, n in (1.0e-14, 1.0e-10, 1.0, 1.0e10)
            family = make(n)
            @test implied_quote(curve(0.035), family, 0.1) ≈ expm1(0.035) rtol = 1.0e-12
            @test ForwardDiff.derivative(z -> implied_quote(curve(z), family, 0.1), 0.035) ≈ exp(0.035) rtol = 1.0e-12
        end
    end

    @testset "loud errors" begin
        @test_throws ArgumentError ForwardDiff.derivative(
            x -> ForwardDiff.derivative(y -> implied_quote(Yield.Constant(Continuous(x + y)), ZCBYield, 2.0), 0.0), 0.03
        )
        @test_throws ArgumentError ForwardDiff.derivative(x -> implied_quote(flat, (q, t) -> ZCBYield(q * x, t), 2.0), 1.0)
        # A quote whose value does not depend on the rate has no implied rate.
        @test_throws Exception implied_quote(flat, (r, t) -> Quote(0.5, Cashflow(1.0, t)), 2.0)
    end
end
