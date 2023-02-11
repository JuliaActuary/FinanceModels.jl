# a Zero Cooupon Bond Quote
@testset "instruments" begin
    @test ZCBPrice(0.9,1) == Quote(0.9,Cashflow(1.,1.))
    @test ZCBYield(0.1,1) == Quote(1/1.1,Cashflow(1.,1.))

    @test Cashflow(1,1) + Cashflow(1,1) == Cashflow(2,1)
    @test Cashflow(1,1) + Cashflow(1,2) == Yields.Composite(Cashflow(1,1), Cashflow(1,2))

    @test ParYield(0.1,1) == Quote(1.,Bond(0.1,Periodic(2),1.))

    @test CMTYield(0.0,0.5) == Quote(1.,Bond(0.00,Periodic(1),0.5))
end