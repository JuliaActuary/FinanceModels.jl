using ForwardDiff

@testset "Credit hazard curves" begin
    @testset "ConstantHazard" begin
        h = Credit.ConstantHazard(0.02)
        @test Credit.hazard(h, 3.0) == 0.02
        @test Credit.cumhazard(h, 5.0) ≈ 0.1
        @test Credit.cumhazard(h, 2.0, 5.0) ≈ 0.06
        @test Credit.survival(h, 5.0) ≈ exp(-0.1)
        @test Credit.survival(h, 2.0, 5.0) ≈ exp(-0.06)
        @test Credit.default_probability(h, 1.0) ≈ 1 - exp(-0.02)
        # small probabilities keep their digits
        tiny = Credit.ConstantHazard(1.0e-12)
        @test Credit.default_probability(tiny, 1.0) ≈ 1.0e-12 rtol = 1.0e-10
        @test Credit.survival(Credit.ConstantHazard(0.0), Inf) == 1.0
        @test Credit.cumhazard(h, Inf) == Inf
        @test_throws DomainError Credit.ConstantHazard(-0.01)
        @test_throws DomainError Credit.ConstantHazard(Inf)
        @test_throws DomainError Credit.survival(h, -1.0)
    end

    @testset "PiecewiseHazard" begin
        h = Credit.PiecewiseHazard([0.01, 0.02, 0.03], [1.0, 5.0, 10.0])
        # right-continuous: at a tenor the intensity is the next interval's
        @test Credit.hazard(h, 0.0) == 0.01
        @test Credit.hazard(h, 1.0) == 0.02
        @test Credit.hazard(h, prevfloat(1.0)) == 0.01
        @test Credit.hazard(h, 10.0) == 0.03
        @test Credit.hazard(h, 50.0) == 0.03   # flat beyond the last tenor
        @test Credit.cumhazard(h, 0.5) ≈ 0.005
        @test Credit.cumhazard(h, 7.0) ≈ 0.01 + 0.02 * 4 + 0.03 * 2
        @test Credit.cumhazard(h, 12.0) ≈ 0.01 + 0.08 + 0.15 + 0.06
        @test Credit.survival(h, 7.0) ≈ exp(-0.15)
        # the interval form agrees with differences of the cumulative form...
        for (a, b) in ((0.0, 0.5), (0.5, 3.0), (2.0, 7.0), (0.3, 12.0), (11.0, 13.0), (5.0, 5.0))
            @test Credit.cumhazard(h, a, b) ≈ Credit.cumhazard(h, b) - Credit.cumhazard(h, a) atol = 1.0e-15
        end
        # ...and keeps its relative precision over a short interval late in the curve
        t1 = 30.0 + 1.0e-9
        Δ = t1 - 30.0   # the interval that t1 represents
        @test Credit.default_probability(h, 30.0, t1) ≈ -expm1(-0.03 * Δ) rtol = 1.0e-14
        # one tenor is a constant hazard
        one_tenor = Credit.PiecewiseHazard([0.02], [3.0])
        @test Credit.survival(one_tenor, 10.0) ≈ Credit.survival(Credit.ConstantHazard(0.02), 10.0)
        # a zero final intensity survives to infinity
        @test Credit.cumhazard(Credit.PiecewiseHazard([0.02, 0.0], [1.0, 2.0]), Inf) ≈ 0.02
        # the inputs are copied
        rates = [0.01, 0.02]
        h2 = Credit.PiecewiseHazard(rates, [1.0, 2.0])
        rates[1] = 1.0
        @test Credit.cumhazard(h2, 1.0) ≈ 0.01
        @test_throws ArgumentError Credit.PiecewiseHazard([0.01, 0.02], [2.0, 1.0])
        @test_throws ArgumentError Credit.PiecewiseHazard([0.01], [0.0])
        @test_throws ArgumentError Credit.PiecewiseHazard(Float64[], Float64[])
        @test_throws DimensionMismatch Credit.PiecewiseHazard([0.01, 0.02], [1.0])
        @test_throws DomainError Credit.PiecewiseHazard([0.01, -0.02], [1.0, 2.0])
        @test_throws DomainError Credit.cumhazard(h, -1.0)
    end

    @testset "derivatives" begin
        survival_at(λ) = Credit.survival(Credit.PiecewiseHazard(λ, [1.0, 5.0, 10.0]), 7.0)
        λ = [0.01, 0.02, 0.03]
        @test ForwardDiff.gradient(survival_at, λ) ≈ -survival_at(λ) .* [1.0, 4.0, 2.0]
        @test ForwardDiff.derivative(l -> Credit.default_probability(Credit.ConstantHazard(l), 2.0), 0.0) ≈ 2.0
    end

    @testset "Assumptions" begin
        a = Credit.Assumptions(Credit.ConstantHazard(0.02); recovery = 0.4)
        @test a.recovery == 0.4
        @test a.hazard isa Credit.ConstantHazard
        @test_throws DomainError Credit.Assumptions(Credit.ConstantHazard(0.02); recovery = 1.2)
        @test_throws DomainError Credit.Assumptions(Credit.ConstantHazard(0.02); recovery = -0.1)
        @test_throws UndefKeywordError Credit.Assumptions(Credit.ConstantHazard(0.02))
    end
end
