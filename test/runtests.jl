using Yields
using Test

@testset "Yields.jl" begin
    @testset "simple rate and forward" begin 
    # Risk Managment and Financial Institutions, 5th ed. Appendix B

        maturity = [0.5, 1.0, 1.5, 2.0]
        zero    = [5.0, 5.8, 6.4, 6.8] ./ 100
        curve = ZeroCurve(zero,maturity)

        
        @test rate(curve,0.5) ≈ 5.0 / 100
        @test rate(curve,2.0) ≈ 6.8 / 100

        @test disc(curve, 1) ≈ 1 / (1 + zero[2])
        @test disc(curve, 2) ≈ 1 / (1 + zero[4]) ^2
    
        # extrapolation
        @test rate(curve,0.0) ≈ 5.0 / 100
        @test rate(curve,4.0) ≈ 6.8 / 100

        @test forward(curve,0.5,1.0) ≈ 6.6 / 100
        @test forward(curve,1.0,1.5) ≈ 7.6 / 100
        @test forward(curve,1.5,2.0) ≈ 8.0 / 100
    end

    @testset "bootstrap treasury" begin
        # https://financetrain.com/bootstrapping-spot-rate-curve-zero-curve/
        maturity = [0.5, 1.0, 1.5, 2.0]
        YTM      = [4.0, 4.3, 4.5, 4.9] ./ 100

        curve = TreasuryYieldCurve(YTM,maturity)

        @test rate(curve, 0.5) ≈ 0.04
        @test rate(curve, 1.0) ≈ 0.043

        # need more future tests, but need validating source...

    end

    @testset "actual cmt treasury" begin
        
        cmt  = [0.12,0.15,0.14,0.17,0.17,0.17,0.19,0.30,0.49,0.64,1.15,1.37] ./ 100
        mats =  [1/12,2/12,3/12,6/12,1,2,3,5,7,10,20,30]

    end

end
