@testset "Equity Models" begin
    m = Equity.BlackScholesMerton(0.01, 0.02, 0.15)

    a = Option.EuroCall(CommonEquity(), 1.0, 1.0)
    b = Option.EuroCall(CommonEquity(), 1.0, 2.0)

    @test pv(m, a) ≈ 0.05410094201902403

    qs = [
        Quote(0.0541, a),
        Quote(0.072636, b),
    ]
    m = Equity.BlackScholesMerton(0.01, 0.02, Volatility.Constant())
    fit(m, qs)
    @test fit(m, qs).σ ≈ 0.15 atol = 1e-4
end