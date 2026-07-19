using PrecompileTools

@setup_workload begin
    rates = [0.01, 0.02, 0.03]
    times = [0.0, 1.0, 2.0, 3.0]
    curve_rates = [first(rates); rates]
    q_rate = ZCBYield(rates)

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

        # Retain one spline fit and one generic optic-based fit so both core
        # Optimization/AD dispatch paths remain warm without repeating them for
        # every interpolation implementation.
        fit(Spline.Linear(), q_rate)
        fit(Yield.NelsonSiegelSvensson(), q_rate)

        present_value(model_rate, Cashflow(1.0, 1.0))
        first(q_rate).instrument |> collect
    end
end
