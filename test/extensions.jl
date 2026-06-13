@testset "Makie Extension" begin
    using CairoMakie
    cfs = [Cashflow(t + randn(), t) for t in 1:10]

    @test stem(cfs) isa Makie.FigureAxisPlot

    b = Bond.Fixed(0.1, Periodic(2), 20)
    proj = Projection(b, NullModel(), CashflowProjection())
    # a stem plot:
    @test stem(proj) isa Makie.FigureAxisPlot

    @test stem(b) isa Makie.FigureAxisPlot
end

@testset "UnicodePlots Extension" begin
    using UnicodePlots
    c = Yield.Constant(0.04)
    # rich display renders a plot...
    plain = sprint(show, MIME("text/plain"), c)
    @test occursin("Yield Curve", plain)
    # ...but plain `print`/string interpolation must NOT (it used to render a
    # 60-char-wide plot into every log line and interpolated string)
    @test !occursin("Yield Curve", string(c))
    @test !occursin("Yield Curve", sprint(show, c))
end
