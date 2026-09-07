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
            Spline.Quadratic(), Spline.Cubic(), Spline.PolynomialSpline(4),
            Spline.BSpline(2), Spline.BSpline(3), Spline.PCHIP(), Spline.Akima(), Spline.MonotoneConvex(),
        )
        @test_throws "Fit.Loss" fit(spline, qs, Fit.Bootstrap())
        # Reject the strategy before even consuming the quotes or entering a solver.
        quotes = (error("quotes should not be consumed") for _ in 1:4)
        @test_throws "Fit.Loss" fit(spline, quotes, Fit.Bootstrap())
    end

    # The documented full-grid alternative must actually fit the quote set.
    # A flat optimizer seed previously stalled PCHIP and Akima at the seed curve.
    quotes = ZCBPrice.([0.98, 0.94, 0.88, 0.75], [1.0, 2.5, 5.0, 10.0])
    for spline in (Spline.PCHIP(), Spline.Akima(), Spline.Cubic(), Spline.BSpline(3), Spline.MonotoneConvex())
        curve = fit(spline, quotes, Fit.Loss(x -> x^2))
        @test maximum(abs(present_value(curve, q.instrument) - q.price) for q in quotes) < 1.0e-7
    end
end
