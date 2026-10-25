using ForwardDiff
using FinanceModels.IntervalSets: (..)

# a curve with one rate up to t = 3 and another after it, to fit with optics into a tuple or a vector
struct TwoRateCurve{T} <: FinanceModels.Yield.AbstractYieldModel
    r::T
end
FinanceCore.discount(m::TwoRateCurve, t) = t <= 3 ? exp(-m.r[1] * t) : exp(-m.r[1] * 3 - m.r[2] * (t - 3))

@testset "fit variables" begin
    reprice(m, qs) = maximum(abs(present_value(m, q.instrument) - q.price) for q in qs)
    OJL = FinanceModels.OptimizationOptimJL
    zcb = ZCBPrice([0.9, 0.8, 0.7, 0.6])
    c0 = Yield.Constant(Continuous(0.05))
    c_fit = fit(c0, zcb)
    two = ZCBPrice.(exp.(-[0.03, 0.06, 0.09, 0.03 * 3 + 0.05 * 2]), [1.0, 2.0, 3.0, 5.0])

    @testset "an interval maps the optimizer's [0, 1] onto it" begin
        # the Uniform cdf and quantile AccessibleModels used, so fits did not change
        vars = (@optic(_.rate.continuous_value) => -1.0 .. 1.0,)
        counts = map(v -> FinanceModels.__value_count(getall(c0, FinanceModels.__variable_optic(v))), vars)
        @test FinanceModels.__unit_start(c0, vars) == [0.525]
        @test FinanceModels.__unit_bounds(vars, counts) == ([0.0], [1.0])
        @test rate(FinanceModels.__set_variables(c0, vars, counts, [0.25], 0).rate) == -0.5
        # a start outside the interval is clamped to its end, then moved to the midpoint
        @test rate(fit(Yield.Constant(Continuous(3.0)), zcb).rate) ≈ rate(c_fit.rate) atol = 1.0e-7
    end

    @testset "the rebuilt model's type is inferred" begin
        for m0 in (Yield.NelsonSiegelSvensson(), ShortRate.Vasicek(0.1, 0.02, 0.01, Continuous(0.01)), TwoRateCurve((0.01, 0.01)))
            vars = m0 isa TwoRateCurve ? (@optic(_.r[∗]) => -1.0 .. 1.0,) : FinanceModels.default_variables(m0)
            counts = map(v -> FinanceModels.__value_count(getall(m0, FinanceModels.__variable_optic(v))), vars)
            u = FinanceModels.__unit_start(m0, vars)
            @inferred FinanceModels.__set_variables(m0, vars, counts, ForwardDiff.Dual{:t}.(u, 1.0), 0)
        end
    end

    @testset "unbounded variables" begin
        lens = @optic(_.rate.continuous_value)
        # an optic alone or in a 1-tuple
        @test isequal(fit(c0, zcb; variables = ((lens,),)), fit(c0, zcb; variables = (lens,)))
        @test rate(fit(c0, zcb; variables = (lens,)).rate) ≈ rate(c_fit.rate) atol = 1.0e-7
        # with every variable unbounded, the optimizer runs without bounds, so Newton() works
        @test rate(fit(c0, zcb; variables = (lens,), optimizer = OJL.Newton()).rate) ≈ rate(c_fit.rate) atol = 1.0e-7
        # a bounded variable runs the optimizer inside Fminbox, which does not take Newton()
        @test_throws ArgumentError fit(c0, zcb; optimizer = OJL.Newton())
        # bounded and unbounded variables together
        m = fit(TwoRateCurve((0.01, 0.01)), two; variables = (@optic(_.r[1]) => 0.0 .. 0.1, @optic(_.r[2])))
        @test reprice(m, two) < 1.0e-8
        # through an FX model's foreign curve
        eurusd = FX.Pair(:EUR, :USD)
        usd = Yield.Constant(Continuous(0.05))
        truth = FX.Forwards(eurusd, 1.1, usd, Yield.Constant(Continuous(0.03)))
        qfx = [FX.Outright(eurusd, forward(truth, t), t) for t in 0.5:0.5:5.0]
        m0 = FX.Forwards(eurusd, 1.1, usd, Yield.Constant(Continuous(0.01)))
        m = fit(m0, qfx; variables = ((@optic(_.foreign.rate.continuous_value),),))
        @test rate(m.foreign.rate) ≈ 0.03 atol = 1.0e-7
    end

    @testset "optics that target several values" begin
        mt = fit(TwoRateCurve((0.01, 0.01)), two; variables = (@optic(_.r[∗]) => -1.0 .. 1.0,))
        @test mt.r isa NTuple{2, Float64}
        @test all(isapprox.(mt.r, (0.03, 0.05); atol = 1.0e-8))
        mv = fit(TwoRateCurve([0.01, 0.01]), two; variables = (@optic(_.r |> Elements()) => -1.0 .. 1.0,))
        @test mv.r isa Vector{Float64}
        @test mv.r ≈ collect(mt.r) atol = 1.0e-12
        @test fit(TwoRateCurve([0.01, 0.01]), two; variables = ((@optic(_.r |> Elements()),),)).r ≈ [0.03, 0.05] atol = 1.0e-7
    end

    @testset "a spec that is not a finite interval or an optic throws" begin
        @test_throws MethodError fit(c0, zcb; variables = (@optic(_.rate.continuous_value) => 0:1,))
    end
end
