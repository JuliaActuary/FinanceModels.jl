# A test contract whose reported maturity precedes its last payment.
struct PaysAfterMaturity{C} <: FinanceCore.AbstractContract
    inner::C
    stated_maturity::Float64
end
FinanceCore.maturity(c::PaysAfterMaturity) = c.stated_maturity
FinanceCore.present_value(model, c::PaysAfterMaturity) = FinanceCore.present_value(model, c.inner)

@testset "Bootstrap strategy boundary" begin
    # Uneven coupon dates expose changed earlier segments; zero-coupon quotes at
    # knots alone cannot detect that a later solve has undone an earlier fit.
    qs = ParYield.([0.02, 0.08, 0.01, 0.06], [1.0, 2.5, 5.0, 10.0])
    for spline in (Spline.Linear(), Spline.BSpline(1))
        curve = fit(spline, reverse(qs), Fit.Bootstrap())
        for i in eachindex(qs)
            prefix = fit(spline, qs[1:i], Fit.Bootstrap())
            for q in qs[1:i]
                @test present_value(prefix, q.instrument) ≈ q.price atol = 1.0e-12
                @test present_value(curve, q.instrument) ≈ q.price atol = 1.0e-12
            end
            for t in range(0, maturity(qs[i]); length = 13)
                @test discount(prefix, t) ≈ discount(curve, t) atol = 1.0e-12
            end
        end
        single = fit(spline, [ZCBPrice(0.95, 2.0)], Fit.Bootstrap())
        @test discount(single, 2.0) ≈ 0.95
        @test zero(single, 0.0) ≈ zero(single, 2.0)
        @test_throws ArgumentError fit(spline, [], Fit.Bootstrap())
        for t in (0.0, -1.0, Inf, NaN)
            @test_throws ArgumentError fit(spline, [ZCBPrice(0.95, t)], Fit.Bootstrap())
        end
    end

    for spline in (
            Spline.Quadratic(), Spline.Cubic(),
            Spline.BSpline(2), Spline.BSpline(3), Spline.PCHIP(), Spline.Akima(),
        )
        @test_throws "Fit.Loss" fit(spline, qs, Fit.Bootstrap())
        # Reject the strategy before even consuming the quotes or entering a solver.
        quotes = (error("quotes should not be consumed") for _ in 1:4)
        @test_throws "Fit.Loss" fit(spline, quotes, Fit.Bootstrap())
    end
    # Monotone convex curves have their own full-curve fit.
    @test_throws "fit(Spline.MonotoneConvex(), quotes)" fit(Spline.MonotoneConvex(), qs, Fit.Bootstrap())

    # A contract paying after its stated maturity is priced off the curve beyond
    # its knot, so later knots move it. The returned curve must not silently
    # misprice it.
    late = Quote(0.9, PaysAfterMaturity(Cashflow(1.0, 3.0), 1.5))
    @test_throws "could not reprice" fit(Spline.Linear(), [ZCBPrice(0.97, 1.0), late, ZCBPrice(0.93, 2.0)], Fit.Bootstrap())

    # Scaling every quote by a notional leaves the curve unchanged: the knot rates are the
    # zero rates the quotes were priced from, at every notional.
    z = [0.03, 0.035, 0.04]
    for n in (1.0e-12, 1.0, 1.0e10)
        scaled = [
            Quote(n * exp(-z[1]), Cashflow(n, 1.0)),
            Quote(
                n * (0.05 * exp(-z[1]) + 1.05 * exp(-2 * z[2])),
                FinanceCore.Composite(Cashflow(0.05n, 1.0), Cashflow(1.05n, 2.0)),
            ),
            Quote(n * exp(-3 * z[3]), Cashflow(n, 3.0)),
        ]
        curve = fit(Spline.Linear(), scaled, Fit.Bootstrap())
        @test knot_rates(curve) ≈ z rtol = 1.0e-12
        @test knot_tenors(curve) == [1.0, 2.0, 3.0]
    end

    # The documented full-grid alternative must actually fit the quote set.
    # A flat optimizer seed previously stalled PCHIP and Akima at the seed curve.
    quotes = ZCBPrice.([0.98, 0.94, 0.88, 0.75], [1.0, 2.5, 5.0, 10.0])
    for spline in (Spline.PCHIP(), Spline.Akima(), Spline.Cubic(), Spline.BSpline(3), Spline.MonotoneConvex())
        curve = fit(spline, quotes, Fit.Loss(x -> x^2))
        @test maximum(abs(present_value(curve, q.instrument) - q.price) for q in quotes) < 1.0e-7
    end
end

@testset "primal root search" begin
    # shared by each bootstrap step and `__implicit_root`
    solve = FinanceModels.__solve_primal_root
    # the secant method finds the root from the start
    @test solve(x -> x^3 - 8, 1.0, (-1.0, 3.0)) ≈ 2.0
    # a secant that diverges (as on a cube root) falls back to the bracketed search
    @test_throws FinanceModels.Roots.ConvergenceFailed FinanceModels.Roots.find_zero(x -> cbrt(x - 0.3), 0.9, FinanceModels.Roots.Order1())
    @test solve(x -> cbrt(x - 0.3), 0.9, (-1.0, 1.0)) ≈ 0.3
    # an error raised while evaluating the function surfaces, and is not retried
    calls = Ref(0)
    @test_throws DomainError solve(x -> (calls[] += 1; x < 0 ? throw(DomainError(x)) : x - 2), -1.0, (0.0, 3.0))
    @test calls[] == 1
    # when both searches fail, the bracketed search's error surfaces
    @test_throws ArgumentError solve(x -> x^2 + 1, 1.0, (-1.0, 1.0))
end
