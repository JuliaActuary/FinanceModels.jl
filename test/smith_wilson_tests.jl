@testset "SmithWilson" begin
    ufr = 0.03
    alpha = 0.1

    # A trivial Qb vector (=0) should result in a flat yield curve
    ufr_curve = SmithWilsonYield(ufr, alpha, [5.0, 7.0], [0.0, 0.0])
    @test discount(ufr_curve, 10.0) == exp(-ufr * 10.0)
end
