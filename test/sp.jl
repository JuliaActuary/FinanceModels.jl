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
    @test pv(Yield.Constant(0.05), p) ≈ 1.0
end

@testset "Quote IRR" begin
    bond = Bond.Fixed(0.05, Periodic(1), 3.0)
    bond_quote = Quote(1.0, bond)
    @test rate(irr(bond_quote)) ≈ 0.05 atol = 1.0e-10

    price = pv(Yield.Constant(0.04), bond)
    bond_quote_at_price = Quote(price, bond)
    @test rate(irr(bond_quote_at_price)) ≈ 0.04 atol = 1.0e-10

    cashflow_quote = Quote(0.95, Cashflow(1.0, 1.0))
    @test rate(irr(cashflow_quote)) ≈ (1 / 0.95 - 1) atol = 1.0e-10

    semi_bond = Bond.Fixed(0.05, Periodic(2), 2.0)
    semi_bond_quote = Quote(1.0, semi_bond)
    expected = irr([Cashflow(-semi_bond_quote.price, 0.0); collect(semi_bond)])
    @test rate(irr(semi_bond_quote)) ≈ rate(expected) atol = 1.0e-10
    @test irr(semi_bond_quote).compounding == Periodic(1)
end

@testset "AbstractArray of contracts" begin
    c = Bond.Fixed(0.05, Periodic(1), 3.0)
    p = Projection([c, c])
    @test length(collect(p)) == 6
    @test pv(Yield.Constant(0.0), p) ≈ 2.3
    p = Projection([c, c |> Map(-)])
    @test length(collect(p)) == 6
    @test pv(Yield.Constant(0.0), p) ≈ 0.0


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


    # a swap value with the same curve used to parameterize should have
    # zero value at inception
    curve = Yield.Constant(0.05)
    swap = InterestRateSwap(curve, 10)
    proj = Projection(swap, Dict("OIS" => curve), CashflowProjection())
    @test pv(0.04, proj |> collect) ≈ 0.0 atol = 1.0e-12
end

@testset "Transduced Contracts" begin

    @test (Bond.Fixed(0.05, Periodic(1), 3) |> Map(-) |> collect) == [
        Cashflow{Float64, Float64}(-0.05, 1.0),
        Cashflow{Float64, Float64}(-0.05, 2.0),
        Cashflow{Float64, Float64}(-1.05, 3.0),
    ]

    @test (Bond.Fixed(0.05, Periodic(1), 3) |> Map(-) |> Map(x -> x * 2) |> collect) == [
        Cashflow{Float64, Float64}(-0.1, 1.0),
        Cashflow{Float64, Float64}(-0.1, 2.0),
        Cashflow{Float64, Float64}(-2.1, 3.0),
    ]


end

@testset "Fit Models" begin

    @testset "Yield Models" begin
        qs = [
            Quote(1.0, Bond.Fixed(0.05, Periodic(1), 3.0)),
            Quote(1.0, Bond.Fixed(0.07, Periodic(1), 3.0)),
        ]
        @test fit(Yield.Constant(), qs, Fit.Loss(x -> abs2(x))).rate ≈ Yield.Constant(0.0602).rate atol = 1.0e-4

    end

end
@testset "Yield Models" begin
    y = Yield.Constant(0.05)
    @test Yield.discount(y, 5) ≈ 1 / (1.05)^5
end
