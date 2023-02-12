@testset "constant curve and rate -> Constant" begin
    yield = Yields.Constant(0.05)
    rate = Yields.Yields.Rate(0.05, Yields.Periodic(1))

    @test Yields.zero(yield, 1) == Yields.Rate(0.05, Yields.Periodic(1))
    @test Yields.zero(Yields.Constant(Yields.Periodic(0.05,2)), 10) == Yields.Rate(0.05, Yields.Periodic(2))
    @test Yields.zero(yield, 5, Yields.Periodic(2)) == convert(Yields.Periodic(2), Yields.Rate(0.05, Yields.Periodic(1)))

    @testset "constant discount time: $time" for time in [0, 0.5, 1, 10]
        @test discount(yield, time) ≈ 1 / (1.05)^time
        @test discount(yield, 0, time) ≈ 1 / (1.05)^time
        @test discount(yield, 3, time + 3) ≈ 1 / (1.05)^time
    end
    @testset "constant discount scalar time: $time" for time in [0, 0.5, 1, 10]
        @test discount(0.05, time) ≈ 1 / (1.05)^time
    end

    @testset "constant accumulation scalar time: $time" for time in [0, 0.5, 1, 10]
        @test accumulation(0.05, time) ≈ 1 * (1.05)^time
    end

    @testset "broadcasting" begin
        @test all(discount.(yield, [1, 2, 3]) .≈ 1 ./ 1.05 .^ (1:3))
        @test all(accumulation.(yield, [1, 2, 3]) .≈ 1.05 .^ (1:3))
    end

    @testset "constant accumulation time: $time" for time in [0, 0.5, 1, 10]
        @test accumulation(yield, time) ≈ 1 * 1.05^time
        @test accumulation(yield, 0, time) ≈ 1 * 1.05^time
        @test accumulation(yield, 3, time + 3) ≈ 1 * 1.05^time
    end


    yield_2x = yield + yield
    yield_add = yield + 0.05
    add_yield = 0.05 + yield
    @testset "constant discount added" for time in [0, 0.5, 1, 10]
        @test discount(yield_2x, time) ≈ 1 / (1.1)^time
        @test discount(yield_add, time) ≈ 1 / (1.1)^time
        @test discount(add_yield, time) ≈ 1 / (1.1)^time
    end

    yield_1bps = yield - Yields.Constant(0.04)
    yield_minus = yield - 0.01
    minus_yield = 0.05 - Yields.Constant(0.01)
    @testset "constant discount subtraction" for time in [0, 0.5, 1, 10]
        @test discount(yield_1bps, time) ≈ 1 / (1.01)^time
        @test discount(yield_minus, time) ≈ 1 / (1.04)^time
        @test discount(minus_yield, time) ≈ 1 / (1.04)^time
    end
end