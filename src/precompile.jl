using PrecompileTools    # this is a small dependency


@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.


    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)


        q_rate = ZCBYield([0.01, 0.02, 0.03])
        ZCBPrice(0.9, 1.0)
        ParYield(0.9, 1.0)
        CMTYield(0.9, 1.0)

        # bootstrap a linear spline yield model
        model_rate = fit(Spline.Linear(), q_rate, Fit.Bootstrap())
        fit(Spline.Quadratic(), q_rate, Fit.Bootstrap())
        fit(Spline.Cubic(), q_rate, Fit.Bootstrap())
        model_rate = fit(Spline.Linear(), q_rate)
        fit(Spline.Quadratic(), q_rate)
        fit(Spline.Cubic(), q_rate)
        fit(Yield.NelsonSiegelSvensson(), q_rate)


        present_value(model_rate, Cashflow(1.0, 1.0))

        first(q_rate).instrument |> collect


    end
end
