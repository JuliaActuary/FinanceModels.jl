
# Cashflows
@testset "Cashflows" begin
    c = Cashflow(1.0, 1.0)

    p = Projection(c)

    @test collect(p) == [c]
    @test collect(c) == [c]
end

# fixed bonds
@testset "Fixed Bonds" begin
    c = Bond.Fixed(0.05, Periodic(1), 3.0)
    p = Projection(c)

    @test collect(p) == [Cashflow(0.05, 1.0), Cashflow(0.05, 2.0), Cashflow(1.05, 3.0)]
    @test collect(c) == [Cashflow(0.05, 1.0), Cashflow(0.05, 2.0), Cashflow(1.05, 3.0)]
    @test collect(Bond.Fixed(0.05, Periodic(1), 2.5)) == [Cashflow(0.05, 0.5), Cashflow(0.05, 1.5), Cashflow(1.05, 2.5)]
    @test collect(Bond.Fixed(0.05, Periodic(1), 1)) == [Cashflow(1.05, 1.0)]

    @test pv(Yield.Constant(0.05), Bond.Fixed(0.05, Periodic(1), 3.0)) ≈ 1.0
end

@testset "Floating Bonds" begin
    p = Projection(
        Bond.Floating(0.02, Periodic(1), 3.0, "SOFR"),
        Dict("SOFR" => Yield.Constant(0.05)),
        CashflowProjection(),
    )

    @test all(collect(p) .≈ [Cashflow(0.07, 1.0), Cashflow(0.07, 2.0), Cashflow(1.07, 3.0)])
end

@testset "Composite Contracts" begin
    a = Bond.Fixed(0.05, Periodic(1), 3.0)
    b = Bond.Fixed(0.1, Periodic(4), 3.0)
    c = FinanceCore.Composite(a, b)

    p = Projection(c)
    @test collect(p) == [collect(a); collect(b)]
end

@testset "Fit Models" begin

    @testset "Yield Models" begin
        qs = [
            Quote(1.0, Bond.Fixed(0.05, Periodic(1), 3.0)),
            Quote(1.0, Bond.Fixed(0.07, Periodic(1), 3.0)),
        ]
        @test fit(Yield.Constant(), qs, Fit.Loss(x -> abs2(x))).rate ≈ Yield.Constant(0.0602).rate atol = 1e-4

    end

end
@testset "Yield Models" begin
    y = Yield.Constant(0.05)
    @test Yield.discount(y, 5) ≈ 1 / (1.05)^5
end