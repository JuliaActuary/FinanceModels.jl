@testset "v8 fitting defaults" begin
    qs = ZCBPrice.([0.97, 0.94, 0.9, 0.86], [1.0, 2.0, 3.0, 5.0])
    spline = Spline.Cubic()
    aligned = ZeroRateCurve(fill(0.02, 4), [1.0, 2.0, 3.0, 5.0], spline; extrapolation = :linear)
    shifted = ZeroRateCurve(fill(0.02, 4), [0.5, 2.0, 4.0, 6.0], spline)
    fewer = ZeroRateCurve(fill(0.02, 3), [1.0, 3.0, 5.0], Spline.Linear())

    @testset "the model and maturities select the method" begin
        @test isequal(fit(spline, qs), fit(spline, qs, Fit.Exact()))
        @test isequal(fit(aligned, reverse(qs)), fit(aligned, qs, Fit.Exact()))
        @test fit(aligned, qs).extrapolation === :linear
        # Equal counts are insufficient: the knot locations must also match.
        @test isequal(fit(shifted, qs), fit(shifted, qs, Fit.LeastSquares()))
        @test isequal(fit(fewer, qs), fit(fewer, qs, Fit.LeastSquares()))
        repeated = [qs; ZCBPrice(0.939, 2.0)]
        @test isequal(fit(fewer, repeated), fit(fewer, repeated, Fit.LeastSquares()))
        for (model, method) in ((aligned, Fit.Exact()), (fewer, Fit.LeastSquares()))
            @test isequal(fit(model, Iterators.Stateful(reverse(qs))), fit(model, reverse(qs), method))
            @test isequal(fit(model, qs; optimizer = nothing, solve_kwargs = nothing), fit(model, qs, method))
        end
        @test isequal(fit(spline, qs; optimizer = nothing, solve_kwargs = nothing), fit(spline, qs, Fit.Exact()))
        @test fit(spline, qs; extrapolation = :flat_zero).extrapolation === :flat_zero
        # A non-square default refit differentiates the least-squares optimality condition.
        rates(p, method...) = rate.(knot_rates(fit(fewer, ZCBPrice.(p, maturity.(qs)), method...)))
        prices = [q.price for q in qs]
        @test ForwardDiff.jacobian(p -> rates(p), prices) ==
            ForwardDiff.jacobian(p -> rates(p, Fit.LeastSquares()), prices)
    end

    pair = FX.Pair(:EUR, :USD)
    domestic = Yield.Constant(Continuous(0.05))
    truth = FX.Forwards(pair, 1.1, domestic, Yield.Constant(Continuous(0.03)))
    fxqs = [FX.Outright(pair, forward(truth, t), t) for t in maturity.(qs)]
    fxmodels = [FX.Forwards(pair, 1.1, domestic, foreign) for foreign in (spline, aligned, fewer)]
    @testset "FX follows its foreign curve" begin
        for (model, method) in zip(fxmodels, (Fit.Exact(), Fit.Exact(), Fit.LeastSquares()))
            @test isequal(fit(model, fxqs).foreign, fit(model, fxqs, method).foreign)
            @test isequal(fit(model, Iterators.Stateful(fxqs)).foreign, fit(model, fxqs, method).foreign)
            @test isequal(fit(model, fxqs; optimizer = nothing, solve_kwargs = nothing).foreign, fit(model, fxqs, method).foreign)
        end
    end

    @testset "settings require an explicit loss method" begin
        settings = (
            (; optimizer = FailingFitOptimizer()),
            (; solve_kwargs = (; g_tol = 1.0e-12)),
            (; solve_kwargs = (;)),
            (; optimizer = nothing, solve_kwargs = (;)),
        )
        for (model, quotes) in [(spline, qs), (aligned, qs), (shifted, qs), (fewer, qs), [(m, fxqs) for m in fxmodels]...], kw in settings
            err = try
                fit(model, quotes; kw...)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("fit(model, quotes, Fit.Loss();", sprint(showerror, err))
        end
        @test_throws FitConvergenceError fit(spline, qs, Fit.Loss(); optimizer = FailingFitOptimizer())
        @test_throws FitConvergenceError fit(aligned, qs, Fit.Loss(); solve_kwargs = (; maxiters = 1))
        @test isequal(fit(fewer, qs, Fit.Loss(); solve_kwargs = (;)), fit(fewer, qs, Fit.Loss()))
        # Generic parametric models retain their optimizer API and default method.
        model = Yield.Constant(0.02)
        optimizer = FinanceModels.default_optimizer(model)
        @test isequal(fit(model, qs), fit(model, qs, Fit.Loss()))
        @test isequal(
            fit(model, qs; optimizer, solve_kwargs = (; g_tol = 1.0e-12)),
            fit(model, qs, Fit.Loss(); optimizer, solve_kwargs = (; g_tol = 1.0e-12)),
        )
    end

    @testset "exact defaults reject an unrepricable quote" begin
        bad = ZCBPrice.([0.97, -0.1], [1.0, 2.0])
        @test_throws FitConvergenceError fit(Spline.Linear(), bad)
        @test_throws FitConvergenceError fit(ZeroRateCurve([0.02, 0.02], [1.0, 2.0], Spline.Linear()), bad)
    end
end
