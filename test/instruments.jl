# a Zero Cooupon Bond Quote
@testset "intruments" begin
    @test ZCBPrice(0.9,1) == Quote(0.9,Cashflow(1.,1.))
    @test ZCBYield(0.1,1) == Quote(1/1.1,Cashflow(1.,1.))

    @test Cashflow(1,1) + Cashflow(1,1) == Cashflow(2,1)
    @test Cashflow(1,1) + Cashflow(1,2) == Yields.Composite(Cashflow(1,1), Cashflow(1,2))
end