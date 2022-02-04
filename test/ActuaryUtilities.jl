using ActuaryUtilities

@testset "ActuaryUtilities.jl integration tests" begin
    cfs = [5, 5, 105]
    times    = [1, 2, 3]

    discount_rates  = [0.03,Yields.Periodic(0.03,1), Yields.Constant(0.03)]

    for d in discount_rates
        @test present_value(d, cfs, times)           ≈ 105.65722270978935
        @test duration(Macaulay(), d, cfs, times)    ≈ 2.86350467067113
        @test duration(d, cfs, times)                ≈ 2.7801016220108057
        @test convexity(d, cfs, times)               ≈ 10.625805482685939
    end
end
    