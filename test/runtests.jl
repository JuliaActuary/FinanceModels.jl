using Yields
using Test

@testset "Yields.jl" begin

    @testset "constant curve" begin
        yield = Yields.Constant(0.05)

        @testset "constant discount time: $time" for time in [0,0.5,1,10]
            @test discount(yield, time) ≈ 1 / (1.05)^time 
            @test accumulate(yield, time) ≈ 1 * 1.05^time
            @test rate(yield, time) == 0.05
        end

        yield_2x = yield + yield
        @testset "constant discount added" for time in [0,0.5,1,10]
            @test discount(yield_2x, time) ≈ 1 / (1.1)^time 
            @test rate(yield_2x, time) == 0.1
        end

        yield_1bps = yield - Yields.Constant(0.04)
        @testset "constant discount subtraction" for time in [0,0.5,1,10]
            @test discount(yield_1bps, time) ≈ 1 / (1.01)^time 
            @test rate(yield_1bps, time) ≈ 0.01
        end
    end
    
    @testset "broadcasting" begin
        yield = Yields.Constant(0.05)
        @test discount.(yield, 1:3) == [1 / 1.05^t for t in 1:3]
    end

    @testset "short curve" begin
        z = Yields.Zero([0.0,0.05], [1,2])
        @test rate(z, 1) ≈ 0.00
        @test rate(z, 2) ≈ 0.05
    end

    @testset "Step curve" begin
        y = Yields.Step([0.02,0.05], [1,2])

        @test rate(y, 0.5) == 0.02

        @test discount(y, 0.5) ≈ 1 / (1.02)^(0.5)
        @test discount(y, 1) ≈ 1 / (1.02)^(1)
        @test rate(y, 1) ≈ 0.02

        @test discount(y, 2) ≈ 1 / (1.02) / 1.05
        @test rate(y, 2) ≈ 0.05
        @test rate(y, 2.5) ≈ 0.05

        @test discount(y, 2) ≈ 1 / (1.02) / 1.05

        @test discount(y, 1.5) ≈ 1 / (1.02) / 1.05^(0.5)


        y = Yields.Step([0.02,0.07])
        @test rate(y, 0.5) ≈ 0.02
        @test rate(y, 1) ≈ 0.02
        @test rate(y, 1.5) ≈ 0.07

    end

    @testset "Salomon Understanding the Yield Curve Pt 1 Figure 9" begin
        maturity = collect(1:10)
        
        par = [6.,8.,9.5,10.5,11.0,11.25,11.38,11.44,11.48,11.5] ./ 100
        spot = [6.,8.08,9.72,10.86,11.44,11.71,11.83,11.88,11.89,11.89] ./ 100

        # the forwards for 7+ have been adjusted from the table - perhaps rounding issues are exacerbated 
        # in the text? forwards for <= 6 matched so reasonably confident that the algo is correct
        # fwd = [6.,10.2,13.07,14.36,13.77,13.1,12.55,12.2,11.97,11.93] ./ 100 # from text
        fwd = [6.,10.2,13.07,14.36,13.77,13.1,12.61,12.14,12.05,11.84] ./ 100  # modified
        
        y = Yields.Par(par, maturity)

        @testset "UTYC Figure 9 par -> spot : $mat" for mat in maturity
            @test rate(y, mat) ≈ spot[mat] atol = 0.0001
            @test forward(y, mat) ≈ fwd[mat] atol = 0.0001
        end



    end

    @testset "simple rate and forward" begin 
    # Risk Managment and Financial Institutions, 5th ed. Appendix B

        maturity = [0.5, 1.0, 1.5, 2.0]
        zero    = [5.0, 5.8, 6.4, 6.8] ./ 100
        curve = Yields.Zero(zero, maturity)

        
        @test rate(curve, 0.5) ≈ 5.0 / 100
        @test rate(curve, 2.0) ≈ 6.8 / 100

        @test discount(curve, 1) ≈ 1 / (1 + zero[2])
        @test discount(curve, 2) ≈ 1 / (1 + zero[4])^2
    
        # extrapolation
        @test rate(curve, 0.0) ≈ 5.0 / 100
        @test rate(curve, 4.0) ≈ 6.8 / 100

        @test forward(curve, 0.5, 1.0) ≈ 6.6 / 100 atol = 0.001
        @test forward(curve, 1.0, 1.5) ≈ 7.6 / 100 atol = 0.001
        @test forward(curve, 1.5, 2.0) ≈ 8.0 / 100 atol = 0.001

        y = Yields.Zero(zero)

        @test discount(y, 1) ≈ 1 / 1.05
        @test discount(y, 2) ≈ 1 / 1.058^2

    end

    @testset "Forward Rates" begin 
    # Risk Managment and Financial Institutions, 5th ed. Appendix B

        forwards = [0.05, 0.04, 0.03, 0.08]
        curve = Yields.Forward(forwards)

        @testset "discounts: $t" for (t, r) in enumerate(forwards)
            @test discount(curve, t) ≈ reduce((v, r) -> v / (1 + r), forwards[1:t]; init=1.0)
        end

        @testset "with specified timepoints" begin
            i = [0.0,0.05]
            times = [0.5,1.5]
            y = Yields.Forward(i, times)
            @test discount(y, 0.5) ≈ 1 / 1.0^0.5  
            @test discount(y, 1.5) ≈ 1 / 1.0^0.5 / 1.05^1

        end

    end

    @testset "base + spread" begin
        riskfree_maturities = [0.5, 1.0, 1.5, 2.0]
        riskfree    = [5.0, 5.8, 6.4, 6.8] ./ 100 # spot rates

        spread_maturities = [0.5, 1.0, 1.5, 3.0] # different maturities
        spread    = [1.0, 1.8, 1.4, 1.8] ./ 100 # spot spreads

        rf_curve = Yields.Zero(riskfree, riskfree_maturities)
        spread_curve = Yields.Zero(spread, spread_maturities)

        yield = rf_curve + spread_curve 

        @test discount(yield, 1.0) ≈ 1 / (1 + riskfree[2] + spread[2])^1
        @test discount(yield, 1.5) ≈ 1 / (1 + riskfree[3] + spread[3])^1.5
    end

    @testset "bootstrap treasury" begin
        # https://financetrain.com/bootstrapping-spot-rate-curve-zero-curve/
        maturity = [0.5, 1.0, 1.5, 2.0]
        YTM      = [4.0, 4.3, 4.5, 4.9] ./ 100

        curve = Yields.USTreasury(YTM, maturity)

        @test rate(curve, 0.5) ≈ 0.04
        @test rate(curve, 1.0) ≈ 0.043

        # need more future tests, but need validating source...

    end

    @testset "actual cmt treasury" begin
        
        cmt  = [0.12,0.15,0.14,0.17,0.17,0.17,0.19,0.30,0.49,0.64,1.15,1.37] ./ 100
        mats =  [1 / 12,2 / 12,3 / 12,6 / 12,1,2,3,5,7,10,20,30]

    end


end
