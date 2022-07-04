@testset "Rate type" begin

    default = Yields.Rate{Float64,typeof(Yields.DEFAULT_COMPOUNDING)}
    
    y = Yields.Constant(0.05)
    @test Yields.__ratetype(y) == Yields.Rate{Float64, Periodic}
    @test Yields.__ratetype(typeof(Periodic(0.05,2))) == Yields.Rate{Float64, Periodic}
    @test Yields.CompoundingFrequency(y) == Periodic(1)
    
    y = Yields.Constant(Continuous(0.05))
    @test Yields.__ratetype(y) == Yields.Rate{Float64, Continuous}
    @test Yields.__ratetype(typeof(Continuous(0.05))) == Yields.Rate{Float64, Continuous}
    @test Yields.CompoundingFrequency(y) == Continuous()
    
    y = Yields.Step([0.02,0.05], [1,2])
    @test Yields.__ratetype(y) == Yields.Rate{Float64, Periodic}
    @test Yields.CompoundingFrequency(y) == Periodic(1)
    
    y = Yields.Forward( [0.01,0.02,0.03] )
    @test Yields.__ratetype(y) == default
    
    rates =[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
    mats = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]
	
    y = Yields.CMT(rates,mats)
    @test Yields.__ratetype(y) == default
    @test Yields.CompoundingFrequency(y) == Continuous()
    
    
    combination = y + y
    @test Yields.__ratetype(combination) == default


end

@testset "type coercion" begin
    @test Yields.__coerce_rate(0.05, Periodic(1)) == Periodic(0.05,1)
    @test Yields.__coerce_rate(Periodic(0.05,12), Periodic(1)) == Periodic(0.05,12)
end

#Issue #117
@testset "DecFP" begin
    import DecFP

    @test Yields.Periodic(1/DecFP.Dec64(1/6)) == Yields.Periodic(6)
    mats = convert.(DecFP.Dec64,[1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30])
    rates = convert.(DecFP.Dec64,[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100)
    y = Yields.CMT(rates,mats)
    @test y isa Yields.AbstractYieldCurve
end