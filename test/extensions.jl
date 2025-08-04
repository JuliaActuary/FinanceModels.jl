@testset "Makie Extension" begin
    using CairoMakie
    cfs = [Cashflow(t + randn(), t) for t in 1:10]

    @test_broken stem(cfs) isa Makie.FigureAxisPlot

    b = Bond.Fixed(0.1, Periodic(2), 20)
    proj = Projection(b, NullModel(), CashflowProjection())
    # a stem plot:
    @test_broken stem(proj) isa Makie.FigureAxisPlot

    @test_broken stem(b) isa Makie.FigureAxisPlot
end
