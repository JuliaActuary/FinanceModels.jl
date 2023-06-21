#Issue #117
@testset "DecFP" begin
    import DecFP

    @test FinanceModels.Periodic(1 / DecFP.Dec64(1 / 6)) == FinanceModels.Periodic(6)
    mats = convert.(DecFP.Dec64, [1 / 12, 2 / 12, 3 / 12, 6 / 12, 1, 2, 3, 5, 7, 10, 20, 30])
    rates = convert.(DecFP.Dec64, [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100)
    y = fit(Spline.Linear(), CMTYield.(rates, mats), Fit.Bootstrap())
    @test y isa FinanceModels.Yield.AbstractYieldModel
end