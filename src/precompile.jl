using PrecompileTools

@setup_workload begin
    rates = [0.01, 0.02, 0.03]
    times = [0.0, 1.0, 2.0, 3.0]
    curve_rates = [first(rates); rates]
    q_rate = ZCBYield(rates)
    q_cmt = CMTYield.([0.04, 0.045, 0.05], [1, 2, 3])          # coupon bonds, integer tenors
    q_par = ParYield.([0.04, 0.045, 0.05], [1.0, 2.0, 3.0])    # coupon bonds, float tenors
    knot_tenors3 = [1.0, 2.0, 3.0]
    cfs = Cashflow.([5.0, 5.0, 105.0], knot_tenors3)

    @compile_workload begin
        # Cover the public quote constructors separately because q_rate was
        # deliberately created in setup (and therefore is not traced).
        ZCBYield(rates)
        ZCBPrice(0.9, 1.0)
        ParYield(0.9, 1.0)
        CMTYield(0.9, 1.0)

        # One representative curve is enough to compile the bootstrap/root-find
        # machinery. The remaining interpolants only need their constructors
        # traced; running a complete bootstrap for each one repeats the same
        # end-to-end root-finding workflow during package precompilation.
        model_rate = fit(Spline.Linear(), q_rate, Fit.Bootstrap())
        Yield.Spline(Spline.Quadratic(), times, curve_rates)
        Yield.Spline(Spline.Cubic(), times, curve_rates)
        Yield.Spline(Spline.BSpline(3), times, curve_rates)

        # Retain one exact spline fit and one generic optic-based fit, which
        # minimizes with Optimization, so both core solve/AD paths remain warm
        # without repeating them for every interpolation implementation.
        fit(Spline.Linear(), q_rate, Fit.Exact())
        fit(Yield.NelsonSiegelSvensson(), q_rate)

        present_value(model_rate, Cashflow(1.0, 1.0))
        present_value(model_rate, [Bond.Fixed(0.05, Periodic(1), 2.0), Cashflow(1.0, 1.0)])
        present_value(Models(model_rate; index = model_rate), Bond.Floating(0.0, Periodic(1), 2.0, :index))
        first(q_rate).instrument |> collect

        # Coupon-bond quotes (CMTYield/ParYield) are the most common market inputs.
        fit(Spline.Linear(), q_cmt, Fit.Exact())
        fit(Spline.Linear(), q_cmt, Fit.Bootstrap())
        fit(Spline.Linear(), q_par, Fit.Bootstrap())
        # Whole-curve fits with the other common interpolants.
        fit(Spline.Cubic(), q_rate, Fit.Loss())
        fit(Spline.Cubic(), q_rate, Fit.Exact())
        fit(Spline.MonotoneConvex(), q_rate, Fit.Loss())
        fit(Spline.MonotoneConvex(), q_rate, Fit.Exact())
        # A least-squares refit of coupon bonds, maturities repeated, on knots of its own.
        fit(ZeroRateCurve(rates, knot_tenors3, Spline.Linear()), [q_cmt; q_par], Fit.LeastSquares())

        # Direct construction and the curve queries on each curve family.
        for c in (
                model_rate,
                ZeroRateCurve(rates, knot_tenors3),
                ZeroRateCurve(rates, knot_tenors3, Spline.Linear()),
                Yield.Spline(Spline.Cubic(), times, curve_rates),
                Yield.Constant(0.03),
            )
            discount(c, 2.5); accumulation(c, 2.5); zero(c, 2.5); forward(c, 1.0, 2.0)
            par(c, 2.0); present_value(c, cfs); present_value(c, Bond.Fixed(0.05, Periodic(2), 3.0))
        end
    end
end
