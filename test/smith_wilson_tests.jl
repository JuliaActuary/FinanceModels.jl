@testset "SmithWilson" begin
    ufr = 0.03
    alpha = 0.1

    # A trivial Qb vector (=0) should result in a flat yield curve
    ufr_curve = SmithWilsonYield(ufr, alpha, [5.0, 7.0], [0.0, 0.0])
    @test discount(ufr_curve, 10.0) == exp(-ufr * 10.0)

    # A single payment at time 4, zero interest
    ci_single = CalibrationInstruments([4.0], reshape([1.0], 1, 1), [1.0])
    curve_with_zero_yield = SmithWilsonYield(ufr, alpha, ci_single)
    @test discount(curve_with_zero_yield, 4.0) == 1.0

    # Still, in the long end it's still just UFR
    @test discount(curve_with_zero_yield, 1000.0) ≈ exp(-ufr * 1000.0) rtol=0.5

    # Three maturities have known discount factors
    # times = [1.0, 2.5, 6.7]
    # prices = [0.9, 0.7, 0.5]
    # cfs = [1 0 0
    #        0 1 0
    #        0 0 1]
    times = [1.0, 2.5, 5.6]
    prices = [0.9, 0.7, 0.5]
    cfs = [1 0 0
           0 1 0
           0 0 1]

    ci_three = CalibrationInstruments(times, cfs, prices)
    curve_three = SmithWilsonYield(ufr, alpha, ci_three)
    @testset "three discounts known" for i in 1:2
        @test discount(curve_three, times[i]) ≈ prices[i] atol = 1e-14
    end

    prices = [1.0, 0.9]
    cfs = [0.1 0.1
           1.0 0.1
           0.0 1.0]
    ci_nondiag = CalibrationInstruments(times, cfs, prices)
    curve_nondiag = SmithWilsonYield(ufr, alpha, ci_nondiag)
    @testset "nondiagonal" for i in 1:2
        @test sum([discount(curve_nondiag, times[tidx]) * cfs[tidx, i] for tidx in 1:3]) ≈ prices[i] atol = 1e-14
    end
       
end
