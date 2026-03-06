@testset "CairnsPritchard" begin

    @testset "Direct construction" begin
        cp = Yield.CairnsPritchard(0.5, 1.5, 0.04, -0.02, -0.01)

        # discount(0) == 1
        @test discount(cp, 0) == 1.0

        # discount * accumulation == 1
        @test discount(cp, 10) ≈ 1 / accumulation(cp, 10)

        # zero rate at known params: r(t) = b₀ + b₁*exp(-c₁*t) + b₂*exp(-c₂*t)
        t = 5.0
        expected_r = 0.04 + (-0.02) * exp(-0.5 * t) + (-0.01) * exp(-1.5 * t)
        @test FinanceModels.zero(cp, t) ≈ Continuous(expected_r)

        # asymptotic: as t→∞, zero rate → b₀
        @test FinanceModels.rate(FinanceModels.zero(cp, 1000.0)) ≈ 0.04 atol = 1e-10
    end

    @testset "DomainError for negative c" begin
        @test_throws DomainError Yield.CairnsPritchard(-1.0, 1.0, 0.0, 0.0, 0.0)
        @test_throws DomainError Yield.CairnsPritchard(1.0, -1.0, 0.0, 0.0, 0.0)
        @test_throws DomainError Yield.CairnsPritchardExtended(-1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0)
        @test_throws DomainError Yield.CairnsPritchardExtended(1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0)
        @test_throws DomainError Yield.CairnsPritchardExtended(1.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0)
    end

    @testset "Fit to Cairns (1998) Figure 1 data" begin
        # Source: Cairns (1998), Figure 1 — digitized via WebPlotDigitizer
        # Cairns, A.J.G. (1998). "Descriptive Bond-Yield and Forward-Rate Models for
        # the British Government Securities Market". British Actuarial Journal, 4(2), 265-321.
        target = [4.106762688183153, 4.264064261935675, 4.403887883049028, 4.550725327846389,
                  4.739607977033906, 4.925342222661497, 5.017616324256814, 5.01874845567391,
                  5.009529671277556, 4.990299610492881, 4.85329014461631, 4.732249122822781,
                  4.625801814105818] ./ 100
        mats = [0.4511004586929275, 0.850176783219295, 1.2049112939094009, 1.6052006163826604,
                2.2324553324737497, 3.277385674943506, 5.2141390634757325, 6.605043375908062,
                8.136251117530515, 10.225033582072788, 17.61340407657162, 23.191577301663614,
                29.55631230653806]

        zqs = ZCBYield.(Continuous.(target), mats)
        c = fit(Yield.CairnsPritchard(), zqs)

        @test discount(c, 0) ≈ 1.0

        @testset "zero rates: $t" for (t, r) in zip(mats, target)
            @test FinanceModels.rate(FinanceModels.zero(c, t)) ≈ r atol = 0.005
        end
    end

    @testset "Fit to Bank of England gilt data (2 components)" begin
        # Source: Bank of England yield curve data (exact date unknown), via PR #128 comment
        gilt = [4.38, 4.83, 5.09, 5.21, 5.26, 5.26, 5.24, 5.20, 5.16, 5.11,
                5.07, 5.03, 4.98, 4.95, 4.91, 4.88, 4.86, 4.83, 4.82, 4.80,
                4.79, 4.78, 4.77, 4.75, 4.74, 4.73, 4.72, 4.70, 4.69, 4.67,
                4.66, 4.64, 4.62, 4.61, 4.59, 4.56, 4.54, 4.52, 4.50, 4.47,
                4.45, 4.42, 4.40, 4.37, 4.34, 4.31, 4.28, 4.26, 4.23, 4.20] ./ 100
        gmats = collect(0.5:0.5:25.0)

        zqs = ZCBYield.(Continuous.(gilt), gmats)
        c = fit(Yield.CairnsPritchard(), zqs)

        @test discount(c, 0) ≈ 1.0

        @testset "zero rates: $t" for (t, r) in zip(gmats, gilt)
            @test FinanceModels.rate(FinanceModels.zero(c, t)) ≈ r atol = 0.01
        end
    end

    @testset "Fit to Bank of England gilt data (3 components)" begin
        # Source: Bank of England yield curve data (exact date unknown), via PR #128 comment
        gilt = [4.38, 4.83, 5.09, 5.21, 5.26, 5.26, 5.24, 5.20, 5.16, 5.11,
                5.07, 5.03, 4.98, 4.95, 4.91, 4.88, 4.86, 4.83, 4.82, 4.80,
                4.79, 4.78, 4.77, 4.75, 4.74, 4.73, 4.72, 4.70, 4.69, 4.67,
                4.66, 4.64, 4.62, 4.61, 4.59, 4.56, 4.54, 4.52, 4.50, 4.47,
                4.45, 4.42, 4.40, 4.37, 4.34, 4.31, 4.28, 4.26, 4.23, 4.20] ./ 100
        gmats = collect(0.5:0.5:25.0)

        zqs = ZCBYield.(Continuous.(gilt), gmats)
        c = fit(Yield.CairnsPritchardExtended(), zqs)

        @test discount(c, 0) ≈ 1.0

        @testset "zero rates: $t" for (t, r) in zip(gmats, gilt)
            @test FinanceModels.rate(FinanceModels.zero(c, t)) ≈ r atol = 0.005
        end
    end

    @testset "Broadcasting" begin
        cp = Yield.CairnsPritchard(0.5, 1.5, 0.04, -0.02, -0.01)
        dfs = discount.(cp, [1, 2, 3])
        @test length(dfs) == 3
        @test all(0 .< dfs .< 1)
    end

    @testset "CairnsPritchardExtended direct construction" begin
        cp = Yield.CairnsPritchardExtended(0.5, 1.5, 3.0, 0.04, -0.02, -0.01, 0.005)

        @test discount(cp, 0) == 1.0
        @test discount(cp, 10) ≈ 1 / accumulation(cp, 10)

        t = 5.0
        expected_r = 0.04 + (-0.02) * exp(-0.5 * t) + (-0.01) * exp(-1.5 * t) + 0.005 * exp(-3.0 * t)
        @test FinanceModels.zero(cp, t) ≈ Continuous(expected_r)

        # asymptotic: as t→∞, zero rate → b₀
        @test FinanceModels.rate(FinanceModels.zero(cp, 1000.0)) ≈ 0.04 atol = 1e-10
    end

end
