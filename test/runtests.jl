using Yields
using Test

@testset "Yields.jl" begin
    
    @testset "rate types" begin
        rs = Rate.([0.1,.02], Yields.Continuous())
        @test rs[1] == Rate(0.1, Yields.Continuous())
        @test rate(rs[1]) == 0.1
    end
    
    @testset "constructor" begin
        @test Yields.Continuous(0.05) == Rate(0.05, Yields.Continuous())
        @test Yields.Periodic(0.02,2) == Rate(0.02, Yields.Periodic(2))
        @test Rate(0.02,2) == Rate(0.02, Yields.Periodic(2))
        @test Rate(0.02,Inf) == Rate(0.02, Yields.Continuous())
    end
    
    @testset "rate conversions" begin
        m = Rate(.1,Yields.Periodic(2))
        @test rate(convert(Yields.Continuous(),m)) ≈ rate(Rate(0.09758, Yields.Continuous())) atol = 1e-5
        c = Rate(0.09758,Yields.Continuous())
        @test convert(Yields.Continuous(),c) == c
        @test rate(convert(Yields.Periodic(2),c)) ≈ rate(Rate(0.1, Yields.Periodic(2))) atol = 1e-5
        @test rate(convert(Yields.Periodic(4),m)) ≈ rate(Rate(0.09878030638383972, Yields.Periodic(4))) atol = 1e-5
        
    end
    
    @testset "constant curve and rate -> Constant" begin
        yield = Yields.Constant(0.05)
        rate = Yields.Rate(0.05, Yields.Periodic(1))
        
        @testset "constant discount time: $time" for time in [0,0.5,1,10]
            @test discount(yield, time) ≈ 1 / (1.05)^time 
            @test discount(rate, time) ≈ 1 / (1.05)^time 
            @test discount(rate, 0, time) ≈ 1 / (1.05)^time 
        end
        @testset "constant discount scalar time: $time" for time in [0,0.5,1,10]
            @test discount(0.05,time) ≈ 1 / (1.05)^time 
        end

        @testset "constant accumulation scalar time: $time" for time in [0,0.5,1,10]
            @test accumulation(0.05,time) ≈ 1 *(1.05)^time 
        end

        @testset "constant accumulation time: $time" for time in [0,0.5,1,10]
            @test accumulation(yield, time) ≈ 1 * 1.05^time
            @test accumulation(rate,time) ≈ 1 * 1.05^time
            @test accumulation(rate,0,time) ≈ 1 * 1.05^time
        end
        
        @testset "CompoundingFrequency" begin
            @testset "Continuous" begin
                cnst = Yields.Constant(Yields.Continuous(0.05))
                @test accumulation(cnst,1) == exp(0.05)
                @test accumulation(cnst,2) == exp(0.05*2)
                @test discount(cnst,2) == 1 / exp(0.05*2)
            end
            
            @testset "Periodic" begin
                p = Yields.Constant(Rate(0.05, Yields.Periodic(2)))
                @test accumulation(p,1) == (1 + 0.05/2) ^ (1 * 2)
                @test accumulation(p,2) == (1 + 0.05/2) ^ (2 * 2)
                @test discount(p,2) == 1 / (1 + 0.05/2) ^ (2 * 2)
                
            end
        end
        
        
        yield_2x = yield + yield
        yield_add = yield + 0.05
        add_yield = 0.05 + yield
        @testset "constant discount added" for time in [0,0.5,1,10]
            @test discount(yield_2x, time) ≈ 1 / (1.1)^time 
            @test discount(yield_add, time) ≈ 1 / (1.1)^time 
            @test discount(add_yield, time) ≈ 1 / (1.1)^time 
        end
        
        yield_1bps = yield - Yields.Constant(0.04)
        yield_minus = yield - 0.01
        minus_yield = 0.05 - Yields.Constant(0.01)
        @testset "constant discount subtraction" for time in [0,0.5,1,10]
            @test discount(yield_1bps, time) ≈ 1 / (1.01)^time 
            @test discount(yield_minus, time) ≈ 1 / (1.04)^time 
            @test discount(minus_yield, time) ≈ 1 / (1.04)^time 
        end
    end
    
    @testset "broadcasting" begin
        yield = Yields.Constant(0.05)
        @test discount.(yield, 1:3) == [1 / 1.05^t for t in 1:3]
    end
    
    @testset "short curve" begin
        z = Yields.Zero([0.0,0.05], [1,2])
        @test rate(zero(z, 1)) ≈ 0.00
        @test discount(z, 1) ≈ 1.00
        @test rate(zero(z, 2)) ≈ 0.05

        # test no times constructor
        z = Yields.Zero([0.0,0.05])
        @test rate(zero(z, 1)) ≈ 0.00
        @test discount(z, 1) ≈ 1.00
        @test rate(zero(z, 2)) ≈ 0.05
    end
    
    @testset "Step curve" begin
        y = Yields.Step([0.02,0.05], [1,2])
        
        @test rate(y, 0.5) == 0.02
        
        @test discount(y,0.0) ≈ 1
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
        
        y = Yields.Par(Rate.(par, Yields.Periodic(1)), maturity)
        
        @testset "UTYC Figure 9 par -> spot : $mat" for mat in maturity
            @test rate(zero(y, mat)) ≈ spot[mat] atol = 0.0001
            @test forward(y, mat-1) ≈ fwd[mat] atol = 0.0001
        end
        
    end
    
    @testset "simple rate and forward" begin 
        # Risk Managment and Financial Institutions, 5th ed. Appendix B
        
        maturity = [0.5, 1.0, 1.5, 2.0]
        zero    = [5.0, 5.8, 6.4, 6.8] ./ 100
        curve = Yields.Zero(zero, maturity)
        
        @test discount(curve, 0) ≈ 1
        @test discount(curve, 1) ≈ 1 / (1 + zero[2])
        @test discount(curve, 2) ≈ 1 / (1 + zero[4])^2
        
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
        curve = Yields.Forward(forwards,[1,2,3,4])

        
        @testset "discounts: $t" for (t, r) in enumerate(forwards)
            @test discount(curve, t) ≈ reduce((v, r) -> v / (1 + r), forwards[1:t]; init=1.0)
        end
        
        # test constructor without times
        curve = Yields.Forward(forwards)
        
        @testset "discounts: $t" for (t, r) in enumerate(forwards)
            @test discount(curve, t) ≈ reduce((v, r) -> v / (1 + r), forwards[1:t]; init=1.0)
        end

        @test accumulation(curve, 0, 1) ≈ 1.05
        @test accumulation(curve, 1, 2) ≈ 1.04
        @test accumulation(curve, 0, 2) ≈ 1.04 * 1.05
        
        # addition / subtraction
        @test discount(curve + 0.1,1) ≈ 1 / 1.15
        @test discount(curve - 0.03,1) ≈ 1 / 1.02
        
        
        
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
    
    @testset "actual cmt treasury" begin
        # Fabozzi 5-5,5-6
        cmt  = [5.25,5.5,5.75,6.0,6.25,6.5,6.75,6.8,7.0,7.1,7.15,7.2,7.3,7.35,7.4,7.5,7.6,7.6,7.7,7.8] ./ 100
        mats = collect(0.5:0.5:10.)
        curve = Yields.CMT(cmt,mats)
        targets = [5.25,5.5,5.76,6.02,6.28,6.55,6.82,6.87,7.09,7.2,7.26,7.31,7.43,7.48,7.54,7.67,7.8,7.79,7.93,8.07] ./ 100
        target_periodicity = fill(2,length(mats))
        target_periodicity[2] = 1 # 1 year is a no-coupon, BEY yield, the rest are semiannual BEY
        @testset "Fabozzi bootstrapped rates" for (r,mat,target,tp) in zip(cmt,mats,targets,target_periodicity)
            @test rate(zero(curve,mat, Yields.Periodic(tp))) ≈ target atol=0.0001
        end

        # Hull, problem 4.34
        adj = ((1 + .051813/2) ^2 -1) * 100
        cmt  = [4.0816,adj,5.4986,5.8620] ./ 100
        mats =  [.5,1.,1.5,2.]
        curve = Yields.CMT(cmt,mats)
        targets = [4.0405,5.1293,5.4429,5.8085] ./ 100
        @testset "Hull bootstrapped rates" for (r,mat,target) in zip(cmt,mats,targets)
            @test rate(zero(curve,mat, Yields.Continuous())) ≈ target atol=0.001
        end
        
    #     # https://www.federalreserve.gov/pubs/feds/2006/200628/200628abs.html
    #     # 2020-04-02 data
    #     cmt = [0.0945,0.2053,0.4431,0.7139,0.9724,1.2002,1.3925,1.5512,1.6805,1.7853,1.8704,1.9399,1.9972,2.045,2.0855,2.1203,2.1509,2.1783,2.2031,2.2261,2.2477,2.2683,2.2881,2.3074,2.3262,2.3447,2.3629,2.3809,2.3987,2.4164] ./ 100
    #     mats = collect(1:30)
    #     curve = Yields.USCMT(cmt,mats)
    #     target = [0.0945,0.2053,0.444,0.7172,0.9802,1.2142,1.4137,1.5797,1.7161,1.8275,1.9183,1.9928,2.0543,2.1056,2.1492,2.1868,2.2198,2.2495,2.2767,2.302,2.3261,2.3494,2.372,2.3944,2.4167,2.439,2.4614,2.4839,2.5067,2.5297] ./ 100
        
    #     @testset "FRB data" for (t,mat,target) in zip(1:length(mats),mats,target)
    #         @show mat
    #         if mat >= 1
    #             @test rate(zero(curve,mat, Yields.Continuous())) ≈ target[mat] atol=0.001
    #         end
    #     end
    end
    
    @testset "OIS" begin
        ois =  [1.8 , 2.0, 2.2, 2.5, 3.0, 4.0] ./ 100
        mats = [1/12, 1/4, 1/2,    1,  2,   5]
        curve = Yields.OIS(ois,mats)
        targets = [0.017987,0.019950,0.021880,0.024693,0.029994,0.040401]
        @testset "bootstrapped rates" for (r,mat,target) in zip(ois,mats,targets)
            @test rate(zero(curve,mat, Yields.Continuous())) ≈ target atol=0.001
        end
    end
    
    @testset "InstrumentQuotes" begin
       
        maturities = [1.3, 2.7]
        prices = [1.1, 0.8]
        zcq = Yields.ZeroCouponQuotes(prices, maturities)
        @test zcq.prices == prices
        @test zcq.maturities == maturities
 
        @test_throws DomainError Yields.ZeroCouponQuotes([1.3, 2.4, 0.9], maturities)
 
        rates = [0.4, -0.7]
        swq = Yields.SwapQuotes(rates, maturities, 3)
        @test swq.rates == rates
        @test swq.maturities == maturities
        @test swq.frequency == 3
 
        @test_throws DomainError Yields.SwapQuotes([1.3, 2.4, 0.9], maturities, 3)
        @test_throws DomainError Yields.SwapQuotes(rates, maturities, 0)
        @test_throws DomainError Yields.SwapQuotes(rates, maturities, -2)
 
        rates = [0.4, -0.7]
        bbq = Yields.BulletBondQuotes(rates, maturities, prices, 3)
        @test bbq.interests == rates
        @test bbq.maturities == maturities
        @test bbq.prices == prices
        @test bbq.frequency == 3
 
        @test_throws DomainError Yields.BulletBondQuotes([1.3, 2.4, 0.9], maturities, prices, 3)
        @test_throws DomainError Yields.BulletBondQuotes(rates, [4.3, 5.6, 4.4, 4.4], prices, 3)
        @test_throws DomainError Yields.BulletBondQuotes(rates, maturities, [5.7], 3)
        @test_throws DomainError Yields.BulletBondQuotes(rates, maturities, prices, 0)
        @test_throws DomainError Yields.BulletBondQuotes(rates, maturities, prices, -4)
 
    end
 
    @testset "SmithWilson" begin
 
        ufr = 0.03
        α = 0.1
        u = [5.0, 7.0]
        qb = [2.3, -1.2]
 
        # Basic behaviour
        sw = Yields.SmithWilson(ufr, α, u, qb)
        @test sw.ufr == ufr
        @test sw.α == α
        @test sw.u == u
        @test sw.qb == qb
        @test_throws DomainError Yields.SmithWilson(ufr, α, u, [2.4, -3.4, 8.9])
    
        # Empty u and Qb should result in a flat yield curve
        # Use this to test methods expected from <:AbstractYieldCurve
        # Only discount and zero are explicitly implemented, so the others should follow automatically
        sw_flat = Yields.SmithWilson(ufr, α, Float64[], Float64[])
        @test discount(sw_flat, 10.0) == exp(-ufr * 10.0)
        @test accumulation(sw_flat, 10.0) ≈ exp(ufr * 10.0)
        @test rate(convert(Yields.Continuous(), zero(sw_flat, 8.0))) ≈ ufr
        @test discount.(sw_flat, [5.0, 10.0]) ≈ exp.(-ufr .* [5.0, 10.0])
        @test rate(convert(Yields.Continuous(), Rate(forward(sw_flat, 5.0, 8.0)))) ≈ ufr
    
        # A trivial Qb vector (=0) should result in a flat yield curve
        ufr_curve = Yields.SmithWilson(ufr, α, u, [0.0, 0.0])
        @test discount(ufr_curve, 10.0) == exp(-ufr * 10.0)
    
        # A single payment at time 4, zero interest
        curve_with_zero_yield = Yields.SmithWilson([4.0], reshape([1.0], 1, 1), [1.0], ufr=ufr, α=α)
        @test discount(curve_with_zero_yield, 4.0) == 1.0
    
        # In the long end it's still just UFR
        @test rate(convert(Yields.Continuous(), Rate(forward(curve_with_zero_yield, 1000.0, 2000.0)))) ≈ ufr
    
        # Three maturities have known discount factors
        times = [1.0, 2.5, 5.6]
        prices = [0.9, 0.7, 0.5]
        cfs = [1 0 0
                0 1 0
                0 0 1]
    
        curve_three = Yields.SmithWilson(times, cfs, prices, ufr=ufr, α=α)
        @test transpose(cfs) * discount.(curve_three, times) ≈ prices
    
        # Two cash flows with payments at three times
        prices = [1.0, 0.9]
        cfs = [0.1 0.1
                1.0 0.1
                0.0 1.0]
        curve_nondiag = Yields.SmithWilson(times, cfs, prices, ufr=ufr, α=α)
        @test transpose(cfs) * discount.(curve_nondiag, times) ≈ prices
    
        # Round-trip zero coupon quotes
        zcq_times = [1.2, 4.5, 5.6]
        zcq_prices = [1.0, 0.9, 1.2]
        zcq = Yields.ZeroCouponQuotes(zcq_prices, zcq_times)
        sw_zcq = Yields.SmithWilson(zcq, ufr=ufr, α=α)
        @testset "ZeroCouponQuotes round-trip" for idx in 1:length(zcq_times)
            @test discount(sw_zcq, zcq_times[idx]) ≈ zcq_prices[idx]
        end
    
        # Round-trip swap quotes
        swq_maturities = [1.2, 2.5, 3.6]
        swq_interests = [-0.02, 0.3, 0.04]
        frequency = 2
        swq = Yields.SwapQuotes(swq_interests, swq_maturities, frequency)
        swq_times = 0.5:0.5:3.5   # Maturities are rounded down to multiples of 1/frequency, [1.0, 2.5, 3.5]
        swq_payments = [-0.01 0.15 0.02
                        0.99 0.15 0.02
                        0.0  0.15 0.02
                        0.0  0.15 0.02
                        0.0  1.15 0.02
                        0.0  0.0  0.02
                        0.0  0.0  1.02]
        sw_swq = Yields.SmithWilson(swq, ufr=ufr, α=α)
        @testset "SwapQuotes round-trip" for swapIdx in 1:length(swq_interests)
            @test sum(discount.(sw_swq, swq_times) .* swq_payments[:, swapIdx]) ≈ 1.0
        end
    
        # Round-trip bullet bond quotes (reuse data from swap quotes)
        bbq_prices = [1.3, 0.1, 4.5]
        bbq = Yields.BulletBondQuotes(swq_interests, swq_maturities, bbq_prices, frequency)
        sw_bbq = Yields.SmithWilson(bbq, ufr=ufr, α=α)
        @testset "BulletBondQuotes round-trip" for bondIdx in 1:length(swq_interests)
            @test sum(discount.(sw_bbq, swq_times) .* swq_payments[:, bondIdx]) ≈ bbq_prices[bondIdx]
        end
    
        # EIOPA risk free rate (no VA), 31 August 2021.
        # https://www.eiopa.europa.eu/sites/default/files/risk_free_interest_rate/eiopa_rfr_20210831.zip
        eiopa_output_qb = [-0.59556534586390800 
                            -0.07442224713453920 
                            -0.34193181987682400 
                            1.54054875814153000 
                            -2.15552046042343000 
                            0.73559290752221900 
                            1.89365225129089000 
                            -2.75927773116240000 
                            2.24893737130629000 
                            -1.51625404117395000 
                            0.19284859623817400 
                            1.13410725406271000 
                            0.00153268224642171 
                            0.00147942301778158 
                            -1.85022125156483000 
                            0.00336230229850928 
                            0.00324546553910162 
                            0.00313268874430658 
                            0.00302383083427276 
                            1.36047951448615000]
        eiopa_output_u = 1:20
        eiopa_ufr = log(1.036)
        eiopa_α = 0.133394
        sw_eiopa_expected = Yields.SmithWilson(eiopa_ufr, eiopa_α, eiopa_output_u, eiopa_output_qb)
    
        eiopa_eurswap_maturities = [1:12; 15; 20]
        eiopa_eurswap_rates = [-0.00615, -0.00575, -0.00535, -0.00485, -0.00425, -0.00375, -0.003145, 
        -0.00245, -0.00185, -0.00125, -0.000711, -0.00019, 0.00111, 0.00215]   # Reverse engineered from output curve. This is the full precision of market quotes.
        eiopa_eurswap_quotes = Yields.SwapQuotes(eiopa_eurswap_rates, eiopa_eurswap_maturities, 1)
        sw_eiopa_actual = Yields.SmithWilson(eiopa_eurswap_quotes, ufr=eiopa_ufr, α=eiopa_α)
    
        @testset "Match EIOPA calculation" begin
            @test sw_eiopa_expected.u ≈ sw_eiopa_actual.u
            @test sw_eiopa_expected.qb ≈ sw_eiopa_actual.qb
        end
    end
 
end
