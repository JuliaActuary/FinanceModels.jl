@testset "Rate type" begin

    default = FinanceModels.Rate{Float64,typeof(FinanceModels.DEFAULT_COMPOUNDING)}
    
    y = FinanceModels.Constant(0.05)
    @test FinanceModels.__ratetype(y) == FinanceModels.Rate{Float64, Periodic}
    @test FinanceModels.__ratetype(typeof(Periodic(0.05,2))) == FinanceModels.Rate{Float64, Periodic}
    @test FinanceModels.FinanceCore.CompoundingFrequency(y) == Periodic(1)
    
    y = FinanceModels.Constant(Continuous(0.05))
    @test FinanceModels.__ratetype(y) == FinanceModels.Rate{Float64, Continuous}
    @test FinanceModels.__ratetype(typeof(Continuous(0.05))) == FinanceModels.Rate{Float64, Continuous}
    @test FinanceModels.FinanceCore.CompoundingFrequency(y) == Continuous()
    
    y = FinanceModels.Step([0.02,0.05], [1,2])
    @test FinanceModels.__ratetype(y) == FinanceModels.Rate{Float64, Periodic}
    @test FinanceModels.FinanceCore.CompoundingFrequency(y) == Periodic(1)
    
    y = FinanceModels.Forward( [0.01,0.02,0.03] )
    @test FinanceModels.__ratetype(y) == default
    
    rates =[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
    mats = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]
	
    y = FinanceModels.CMT(rates,mats)
    @test FinanceModels.__ratetype(y) == default
    @test FinanceModels.FinanceCore.CompoundingFrequency(y) == Continuous()
    
    
    combination = y + y
    @test FinanceModels.__ratetype(combination) == default


end

@testset "type coercion" begin
    @test FinanceModels.__coerce_rate(0.05, Periodic(1)) == Periodic(0.05,1)
    @test FinanceModels.__coerce_rate(Periodic(0.05,12), Periodic(1)) == Periodic(0.05,12)
end

#Issue #117
@testset "DecFP" begin
    import DecFP

    @test FinanceModels.Periodic(1/DecFP.Dec64(1/6)) == FinanceModels.Periodic(6)
    mats = convert.(DecFP.Dec64,[1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30])
    rates = convert.(DecFP.Dec64,[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100)
    y = FinanceModels.CMT(rates,mats)
    @test y isa FinanceModels.AbstractYieldCurve
end