@testset "Derivatives" begin

    @testset "Euro Options" begin
        # tested against https://option-price.com/index.php
        params = (S=1.0, K=1.0, τ=1, r=0.05, σ=0.25, q=0.0)

        @test FinanceModels.eurocall(; params...) ≈ 0.12336 atol = 1e-5
        @test FinanceModels.europut(; params...) ≈ 0.07459 atol = 1e-5

        params = (S=1.0, K=1.0, τ=1, r=0.05, σ=0.25, q=0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.105493 atol = 1e-5
        @test FinanceModels.europut(; params...) ≈ 0.086277 atol = 1e-5

        params = (S=1.0, K=0.5, τ=1, r=0.05, σ=0.25, q=0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.49494 atol = 1e-5
        @test FinanceModels.europut(; params...) ≈ 0.00011 atol = 1e-5

        params = (S=1.0, K=0.5, τ=1, r=0.05, σ=0.25, q=0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.49494 atol = 1e-5
        @test FinanceModels.europut(; params...) ≈ 0.00011 atol = 1e-5

        params = (S=1.0, K=0.5, τ=0, r=0.05, σ=0.25, q=0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.5 atol = 1e-5
        @test FinanceModels.europut(; params...) ≈ 0.0 atol = 1e-5

        params = (S=1.0, K=1.5, τ=0, r=0.05, σ=0.25, q=0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.0 atol = 1e-5
        @test FinanceModels.europut(; params...) ≈ 0.5 atol = 1e-5

    end
end



@testset "Equity Models" begin
    m = Equity.BlackScholesMerton(0.01, 0.02, 0.15)

    a = Option.EuroCall(CommonEquity(), 1.0, 1.0)
    b = Option.EuroCall(CommonEquity(), 1.0, 2.0)

    @test pv(m, a) ≈ 0.05410094201902403

    ap = Option.EuroPut(CommonEquity(), 1.0, 1.0)

    # put call parity with continuous dividend yield
    # https://web.ma.utexas.edu/users/mcudina/m339d_lecture6.pdf, pg 2
    @test pv(m, a) - pv(m, ap) ≈ 1.0 * exp(-0.02) - exp(-0.01) * 1.0


    # test fitting volatility
    qs = [
        Quote(0.0541, a),
        Quote(0.072636, b),
    ]
    m = Equity.BlackScholesMerton(0.01, 0.02, Volatility.Constant())
    fit(m, qs)
    @test fit(m, qs).σ ≈ 0.15 atol = 1e-4

    # test constructor with Rates
    @test Equity.BlackScholesMerton(Continuous(0.01), Continuous(0.02), 0.15) isa Equity.BlackScholesMerton
    m = Equity.BlackScholesMerton(Periodic(0.01, 1), Periodic(0.02, 1), 0.15)
    # periodic rates should convert to Continuous equivalent that's higher than the Continuous rate
    @test m.r < 0.01
    @test m.q < 0.02


end