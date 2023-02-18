
@testset "Forward Starting" begin
    maturity = [0.5, 1.0, 1.5, 2.0]
    zeros = [5.0, 5.8, 6.4, 6.8] ./ 100
    obs = ZCBYield.(zeros, maturity)
    c = curve(obs)

    fwd = Yields.ForwardStarting(c, 1.0)
    @test discount(fwd, 0) ≈ 1
    @test discount(fwd, 0.5) ≈ discount(c, 1, 1.5)
    @test discount(fwd, 1) ≈ discount(c, 1, 2)
    @test accumulation(fwd, 1) ≈ accumulation(c, 1, 2)

    @test zero(fwd,1) ≈ forward(c,1,2)
    @test zero(fwd,1) ≈ Continuous()(forward(c,1,2))
end