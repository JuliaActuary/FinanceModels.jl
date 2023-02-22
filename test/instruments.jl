# a Zero Cooupon Bond Quote
@testset "instruments" begin
    @test ZCBPrice(0.9,1.) == Quote(0.9,Cashflow(1.,1.))
    @test ZCBYield(0.1,1.) == Quote(1/1.1,Cashflow(1.,1.))

    @test Cashflow(1.,1.) + Cashflow(1.,1.) == Cashflow(2.,1.)
    @test Cashflow(1.,1.) + Cashflow(1.,2.) == Yields.Composite(Cashflow(1.,1.), Cashflow(1.,2.))

    @test ParYield(0.1,1.) == Quote(1.,Bond(0.1,Periodic(2),1.))

    @test CMTYield(0.05,0.5) ≈ Quote(1/(1.05)^0.5,Bond(0.00,Periodic(1),0.5))
    @test CMTYield(0.05,1.) ≈ Quote(1/(1.05)^1,Bond(0.00,Periodic(1),1.))
    @test CMTYield(0.05,2.) ≈ Quote(1.,Bond(0.05,Periodic(2),2.))
    

    @testset "Forward" begin
        fy = ForwardYield.([0.01,0.02],[1.,2.])
        @test first(fy) == Quote(1/1.01,Forward(0.0,Cashflow(1.,1.)))
        @test last(fy) == Quote(1/1.02,Forward(1.0,Cashflow(1.,1.)))
    end

    @testset "Iteration" begin
        #TODO
    end

    @testset "CashflowMatrix & TimeMatrix" begin
        @testset "basic cf" begin
            times = 1:4
            obs = Cashflow.(1.0,times)
            @test Yields.cashflow_matrix(obs) == [i == j ? 1.0 : 0.0 for i in times, j in times]
            @test Yields.timesteps(obs) == times
        end

        @testset "bonds" begin
            obs = Bond.(0.1,Periodic(1),2.:3)
            @test Yields.cashflow_matrix(obs) == [0.1 1.1 0.0; 0.1 0.1 1.1]'
            @test Yields.timesteps(obs) == 1.:3

            obs = Bond.(0.1,Periodic(2),2.:3)
            @test Yields.timesteps(obs) == 0.5:0.5:3
            @test Yields.cashflow_matrix(obs) == [0.05 0.05 0.05 1.05 0.0 0.0; 0.05 0.05 0.05 0.05 0.05 1.05]'
        end
    end
end