@testset "Monotone Convex" begin

    # Interpolation of the yield curve, Gustaf Dehlbom
    # http://uu.diva-portal.org/smash/get/diva2:1477828/FULLTEXT01.pdf

    prices = [0.98, 0.955, 0.92, 0.88, 0.830]
    times = [1, 2, 3, 4, 5]
    quotes = ZCBPrice.(prices, times)
    rates = @. -log(prices) / times
    c = Yield.MonotoneConvex(rates, times)
    # fit(Yield.MonotoneConvex(), quotes)

    f, fᵈ = Yield.__monotone_convex_fs(rates, times)

    @test all(f .== c.f)
    @test all(fᵈ .== c.fᵈ)

    @test fᵈ[1] ≈ 0.0202 atol = 0.0001
    @test fᵈ[2] ≈ 0.0258 atol = 0.0001
    @test fᵈ[3] ≈ 0.0373 atol = 0.0001
    @test fᵈ[4] ≈ 0.0445 atol = 0.0001
    @test fᵈ[5] ≈ 0.0585 atol = 0.0001

    @test f[1] ≈ 0.0188 atol = 0.0001
    @test f[2] ≈ 0.023 atol = 0.0001
    @test f[3] ≈ 0.0316 atol = 0.0001
    @test f[4] ≈ 0.0409 atol = 0.0001
    @test f[5] ≈ 0.0515 atol = 0.0001
    @test f[6] ≈ 0.0620 atol = 0.0001

    @test FinanceModels.Yield.g(0.5, 0.018793076350927487, 0.023021969250703423, 0.020202707317519466) ≈ 0.0042 * (0.5)^2 - 0.0014 atol = 0.0001
    @test FinanceModels.Yield.g(0, 0.023021969250703423, 0.03158945081076577, 0.02584123118388738) ≈ -0.0028 atol = 0.0001
    @test FinanceModels.Yield.g(0.5, 0.023021969250703423, 0.03158945081076577, 0.02584123118388738) ≈ 0.0087 * (0.5)^2 − 0.0002 * 0.5 - 0.0028 atol = 0.0001
    @test FinanceModels.Yield.g(0.5, 0.03158945081076577, 0.04089471650423902, 0.03733767043764417) ≈ -0.0063 * (0.5)^2 + 0.0156 * 0.5 - 0.0057 atol = 0.0001
    @test FinanceModels.Yield.g(0.5, 0.04089471650423902, 0.051473984626221235, 0.04445176257083387) ≈ 0.0102 * (0.5)^2 + 0.0004 * 0.5 - 0.0036 atol = 0.0001


    function r(t)
        if 0 <= t <= 1
            return 0.0014t^2 + 0.0188
        elseif 1 <= t <= 1.0233
            return -0.0028 / t + 0.0230
        elseif 1.0233 <= t <= 2
            return 0.0029t^2 - 0.0088t - 0.0058 / t + 0.0319
        elseif 2 <= t <= 3
            return -0.0022t^2 + 0.0212t + 0.0324 / t - 0.0268
        elseif 3 <= t <= 4
            return 0.0031t^2 - 0.0274t - 0.1188 / t + 0.1217
        elseif 4 <= t <= 5
            return -0.0035t^2 + 0.0525t + 0.314 / t - 0.2005
        else
            error("t is out of the defined range")
        end
    end


    # Test zero rates against reference function (tolerance needed due to rounded coefficients in r(t))
    @testset for t in range(0, 5, 30)
        @test rate(zero(c, t)) ≈ r(t) atol = 0.0001
    end

    # Test that zero rates match exactly at the knot points
    @testset "knot points" begin
        for (i, t) in enumerate(times)
            @test rate(zero(c, t)) ≈ rates[i]
        end
    end

    # Test that discount factors match original prices at knot points
    @testset "discount factors" begin
        for (i, t) in enumerate(times)
            @test discount(c, t) ≈ prices[i]
        end
    end


    # TODO: More tests from
    # https://repository.up.ac.za/bitstream/handle/2263/25882/dissertation.pdf?sequence=1&isAllowed=y
end