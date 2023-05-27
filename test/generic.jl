# test extending type and that various generic methods are defined 
# to fulfill the API of AbstractYieldCurve
@testset "generic methods and type extensions" begin
    struct MyYield <: FinanceModels.AbstractYieldCurve
        rate
    end

    FinanceModels.discount(c::MyYield,to) = exp(-c.rate * to)
    # Base.zero(c::MyYield,to) = Continuous(c.rate)

    my = MyYield(0.05)
    @test zero(my,1) ≈ Continuous(0.05)
    @test FinanceModels.forward(my,1) ≈ Continuous(0.05)
    @test FinanceModels.forward(my,1,2) ≈ Continuous(0.05)
    @test discount(my,1,2) ≈ exp(-0.05*1)
    @test accumulation(my,1,2) ≈ exp(0.05*1)
    @test accumulation(my,1) ≈ exp(0.05*1)
    @test FinanceModels.__ratetype(my) == FinanceModels.Rate{Float64,Continuous}
    @test FinanceModels.FinanceCore.CompoundingFrequency(my) == Continuous()



    @test FinanceModels.par(my,1) |> FinanceModels.rate > 0 
end