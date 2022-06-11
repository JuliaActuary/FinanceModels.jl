@testset "Rates" begin
    @testset "rate types" begin
        rs = Yields.Rate.([0.1, 0.02], Yields.Continuous())
        @test rs[1] == Yields.Rate(0.1, Yields.Continuous())
        @test rs[1] == Yields.Continuous(0.1)
        @test rate(rs[1]) == 0.1
    end

    @testset "constructor" begin
        @test Yields.Continuous(0.05) == Yields.Rate(0.05, Yields.Continuous())
        @test Yields.Periodic(0.02, 2) == Yields.Rate(0.02, Yields.Periodic(2))
        @test Yields.Rate(0.02, 2) == Yields.Rate(0.02, Yields.Periodic(2))
        @test Yields.Rate(0.02, Inf) == Yields.Rate(0.02, Yields.Continuous())
    end

    @testset "rate conversions" begin
        m = Yields.Rate(0.1, Yields.Periodic(2))
        @test rate(convert(Yields.Continuous(), m)) ≈ rate(Yields.Rate(0.09758, Yields.Continuous())) atol = 1e-5
        c = Yields.Rate(0.09758, Yields.Continuous())
        @test convert(Yields.Continuous(), c) == c
        @test rate(convert(Yields.Periodic(2), c)) ≈ rate(Yields.Rate(0.1, Yields.Periodic(2))) atol = 1e-5
        @test rate(convert(Yields.Periodic(4), m)) ≈ rate(Yields.Rate(0.09878030638383972, Yields.Periodic(4))) atol = 1e-5

    end

    @testset "rate equality" begin
        a = Yields.Periodic(0.02, 2)
        a_eq = Yields.Periodic((1+.02/2)^2-1, 1)
        b = Yields.Periodic(0.03, 2)
        c = Yields.Continuous(0.02)

        @test a == a
        @test !(a == a_eq) # not equal due to floating point error
        @test a ≈ a_eq
        @test a != b
        @test ~(a ≈ b)
        @test (a ≈ a)
        @test ~(a ≈ c)

    end

    @testset "discounting and accumulation" for t in [-1.3, 2.46, 6.7]
        
        unspecified_rate = 0.035
        periodic_rate = Yields.Periodic(0.02, 2)
        continuous_rate = Yields.Continuous(0.03)

        @test discount(unspecified_rate, t) ≈ (1 + 0.035)^(-t)
        @test discount(periodic_rate, t) ≈ (1 + 0.02 / 2)^(-t * 2)
        @test discount(continuous_rate, t) ≈ exp(-0.03 * t)

        @test accumulation(unspecified_rate, t) ≈ (1 + 0.035)^t
        @test accumulation(periodic_rate, t) ≈ (1 + 0.02 / 2)^(t * 2)
        @test accumulation(continuous_rate, t) ≈ exp(0.03 * t)

    end

    @testset "rate over interval" begin
        
        from = -0.45
        to = 3.4
        rate = 0.15

        @test discount(rate, from, to) ≈ discount(rate, to - from)
        @test accumulation(rate, from, to) ≈ accumulation(rate, to - from)
        
    end

    @testset "AbstractYield Interface" begin
        c = Continuous(0.03)
        p = Periodic(0.04,2)

        @test zero(c,2) ≈ c
        @test zero(p,2) ≈ p
        @test forward(c,2) ≈ c
        @test forward(p,2) ≈ p

        @test discount(c,2) ≈ exp(-2*0.03)
        @test discount(p,2) ≈ 1 / (1 + .04/2)^(2*2)

        @test discount(c,2) ≈ 1 / accumulation(c,2)
        @test discount(p,2) ≈ 1 / accumulation(p,2)


    end

    @testset "rate algebra" begin

        a = 0.03
        b = 0.02
        
        @testset "addition" begin
            c(x) = Yields.Continuous(x)
            p(x) = Yields.Periodic(x, 1)

            @test c(a) + b ≈ Yields.Continuous(0.05)
            @test a + c(b) ≈ Yields.Continuous(0.05)
            
            @test p(a) + b ≈ Yields.Periodic(0.05,1)
            @test a + p(b) ≈ Yields.Periodic(0.05,1)
        end

        @testset "multiplication" begin
            c(x) = Yields.Continuous(x)
            p(x) = Yields.Periodic(x, 1)

            @test c(a) * b ≈ Yields.Continuous(a * b)
            @test a * c(b) ≈ Yields.Continuous(a * b)
            
            @test p(a) * b ≈ Yields.Periodic(a * b,1)
            @test a * p(b) ≈ Yields.Periodic(a * b,1)
        end

        @testset "division" begin
            c(x) = Yields.Continuous(x)
            p(x) = Yields.Periodic(x, 1)

            @test c(a) / b ≈ Yields.Continuous(a / b)
            @test_throws MethodError a / c(b) ≈ Yields.Continuous(a / b)
            
            @test p(a) / b ≈ Yields.Periodic(a / b,1)
            @test_throws MethodError a / p(b) ≈ Yields.Periodic(a / b,1)
        end

        @testset "subtraction" begin
            c(x) = Yields.Continuous(x)
            p(x) = Yields.Periodic(x, 1)

            @test c(a) - b ≈ Yields.Continuous(0.01)
            @test a - c(b) ≈ Yields.Continuous(0.01)
            
            @test p(a) - b ≈ Yields.Periodic(0.01,1)
            @test a - p(b) ≈ Yields.Periodic(0.01,1)
        end

        @testset "Rate and Rate" begin
            r = Yields.Periodic(0.04,2) - Yields.Periodic(0.01,2) 
            @test r ≈ Yields.Periodic(0.03,2)
            r = Yields.Periodic(0.04,2) + Yields.Periodic(0.01,2) 
            @test r ≈ Yields.Periodic(0.05,2)

            @test Yields.Periodic(0.04,1) > Yields.Periodic(0.03,2)
            @test Yields.Periodic(0.03,1) < Yields.Periodic(0.04,2)
            @test ~(Yields.Periodic(0.04,1) < Yields.Periodic(0.03,2))
            @test ~(Yields.Periodic(0.03,1) > Yields.Periodic(0.04,2))

            @test Yields.Periodic(0.03,1) < Yields.Periodic(0.03,2)
            @test Yields.Periodic(0.03,100) < Yields.Continuous(0.03)
            @test Yields.Periodic(0.03,2) > Yields.Periodic(0.03,1)
            @test Yields.Continuous(0.03) > Yields.Periodic(0.03,100) 
        end
    end

end