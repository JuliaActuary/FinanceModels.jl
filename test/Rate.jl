@testset "Rates" begin
    @testset "rate types" begin
        rs = Rate.([0.1, 0.02], Yields.Continuous())
        @test rs[1] == Rate(0.1, Yields.Continuous())
        @test rate(rs[1]) == 0.1
    end

    @testset "constructor" begin
        @test Yields.Continuous(0.05) == Rate(0.05, Yields.Continuous())
        @test Yields.Periodic(0.02, 2) == Rate(0.02, Yields.Periodic(2))
        @test Rate(0.02, 2) == Rate(0.02, Yields.Periodic(2))
        @test Rate(0.02, Inf) == Rate(0.02, Yields.Continuous())
    end

    @testset "rate conversions" begin
        m = Rate(0.1, Yields.Periodic(2))
        @test rate(convert(Yields.Continuous(), m)) ≈ rate(Rate(0.09758, Yields.Continuous())) atol = 1e-5
        c = Rate(0.09758, Yields.Continuous())
        @test convert(Yields.Continuous(), c) == c
        @test rate(convert(Yields.Periodic(2), c)) ≈ rate(Rate(0.1, Yields.Periodic(2))) atol = 1e-5
        @test rate(convert(Yields.Periodic(4), m)) ≈ rate(Rate(0.09878030638383972, Yields.Periodic(4))) atol = 1e-5

    end

    @testset "rate equality" begin
        a = Yields.Periodic(0.02, 2)
        b = Yields.Periodic(0.03, 2)
        c = Yields.Continuous(0.02)

        @test a == a
        @test a != b
        @test ~(a ≈ b)
        @test (a ≈ a)
        @test ~(a ≈ c)

    end
end