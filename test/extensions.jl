@testset "Makie Extension" begin
    using CairoMakie
    cfs = [Cashflow(t + randn(), t) for t in 1:10]

    @test stem(cfs) isa Makie.FigureAxisPlot

    b = Bond.Fixed(0.10, Periodic(2), 20)
    proj = Projection(b, NullModel(), CashflowProjection())
    # a stem plot:
    @test stem(proj) isa Makie.FigureAxisPlot

    @test stem(b) isa Makie.FigureAxisPlot
end