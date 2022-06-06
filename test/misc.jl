@testset "Rate type" begin

    default = Yields.Rate{Float64,typeof(Yields.DEFAULT_COMPOUNDING)}
    
    y = Yields.Constant(0.05)
    @test Yields.__ratetype(y) == Yields.Rate{Float64, Periodic}

    y = Yields.Constant(Continuous(0.05))
    @test Yields.__ratetype(y) == Yields.Rate{Float64, Continuous}
    
    y = Yields.Step([0.02,0.05], [1,2])
    @test Yields.__ratetype(y) == Yields.Rate{Float64, Periodic}
    
    y = Yields.Forward( [0.01,0.02,0.03] )
    @test Yields.__ratetype(y) == default
    
    rates =[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
    mats = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]
	
    y = Yields.CMT(rates,mats)
    @test Yields.__ratetype(y) == default
    
    
    combination = y + y
    @test Yields.__ratetype(combination) == default


end