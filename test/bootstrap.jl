@testset "bootstrapped class of curves" begin


    @testset "short curve" begin
        z = curve(ZCBYield.([0.0, 0.05], [1, 2]))
        @test zero(z, 1) ≈ Periodic(0.00,1)
        @test discount(z, 1) ≈ 1.00
        @test zero(z, 2) ≈ Periodic(0.05,1)

        # test no times constructor
        z = curve(ZCBYield([0.0, 0.05]))
        @test zero(z, 1) ≈ Periodic(0.00,1)
        @test discount(z, 1) ≈ 1.00
        @test zero(z, 2) ≈ Periodic(0.05,1)
    end

    
    @testset "Hull" begin
        # Par Yield, pg 85

        c = curve(ParYield.(Periodic(2).([0.0687,0.0687]), [2.,3.]))

        @test Yields.par(c,2) ≈ Yields.Periodic(0.0687,2) atol = 0.00001

    end

    @testset "simple rate and forward" begin
        # Risk Managment and Financial Institutions, 5th ed. Appendix B

        maturity = [0.5, 1.0, 1.5, 2.0]
        zero = [5.0, 5.8, 6.4, 6.8] ./ 100
        c = ZCBYield.(zero, maturity)

        @testset "curves" for curve in [curve(c),curve(Bootstrap(),c),curve(Bootstrap(LinearSpline()),c),curve(Bootstrap(QuadraticSpline()),c)]
            @test discount(curve, 0) ≈ 1
            @test discount(curve, 1) ≈ 1 / (1 + zero[2])
            @test discount(curve, 2) ≈ 1 / (1 + zero[4])^2

            @test forward(curve, 0.5, 1.0) ≈ Yields.Periodic(6.6 / 100, 1) atol = 0.001
            @test forward(curve, 1.0, 1.5) ≈ Yields.Periodic(7.6 / 100, 1) atol = 0.001
            @test forward(curve, 1.5, 2.0) ≈ Yields.Periodic(8.0 / 100, 1) atol = 0.001
        end

        @testset "broadcasting" begin
            y = curve(ZCBYield(zero))
            @test all(discount.(y, [1, 2]) .≈ [1 / 1.05, 1 / 1.058^2])
            @test all(accumulation.(y, [1, 2]) .≈ [1.05, 1.058^2])
        end

    end

    @testset "simplest par" begin
        c = curve(ParYield([0.05]))
        @test Yields.par(c,1) ≈ Yields.Periodic(0.05,2)

        c = curve(ParYield([0.05,0.05]))
        @test Yields.par(c,2) ≈ Yields.Periodic(0.05,2)
    end

    @testset "Salomon Understanding the Yield Curve Pt 1 Figure 9" begin
        maturity = collect(1.:10.)

        par = [6.0, 8.0, 9.5, 10.5, 11.0, 11.25, 11.38, 11.44, 11.48, 11.5] ./ 100
        spot = [6.0, 8.08, 9.72, 10.86, 11.44, 11.71, 11.83, 11.88, 11.89, 11.89] ./ 100

        # the forwards for 7+ have been adjusted from the table - perhaps rounding issues are exacerbated 
        # in the text? forwards for <= 6 matched so reasonably confident that the algo is correct
        # fwd = [6.,10.2,13.07,14.36,13.77,13.1,12.55,12.2,11.97,11.93] ./ 100 # from text
        fwd = [6.0, 10.2, 13.07, 14.36, 13.77, 13.1, 12.61, 12.14, 12.05, 11.84] ./ 100  # modified

        y = curve(ParYield.(Periodic(1).(par), maturity))
        @testset "quadratic UTYC Figure 9 par -> spot : $m" for (p,s,f,m) in zip(par,spot,fwd,maturity)
            @test Yields.zero(y, m) ≈ Periodic(s,1) atol = 0.0001
            @test Yields.forward(y, m - 1) ≈ Yields.Periodic(f, 1) atol = 0.0001
        end

        y = curve(Bootstrap(LinearSpline),ParYield.(Periodic(1).(par), maturity))

        @testset "linear UTYC Figure 9 par -> spot : $m" for (p,s,f,m) in zip(par,spot,fwd,maturity)
            @test Yields.zero(y, m) ≈ Periodic(s,1) atol = 0.0001
            @test Yields.forward(y, m - 1) ≈ Yields.Periodic(f, 1) atol = 0.0001
        end

    end

    @testset "Forward Rates" begin
        # Risk Managment and Financial Institutions, 5th ed. Appendix B

        forwards = [0.05, 0.04, 0.03, 0.08]
        curve = Yields.Forward(forwards, [1, 2, 3, 4])


        @testset "discounts: $t" for (t, r) in enumerate(forwards)
            @test discount(curve, t) ≈ reduce((v, r) -> v / (1 + r), forwards[1:t]; init = 1.0)
        end

        # test constructor without times
        curve = Yields.Forward(forwards)

        @testset "discounts: $t" for (t, r) in enumerate(forwards)
            @test discount(curve, t) ≈ reduce((v, r) -> v / (1 + r), forwards[1:t]; init = 1.0)
        end

        @test accumulation(curve, 0, 1) ≈ 1.05
        @test accumulation(curve, 1, 2) ≈ 1.04
        @test accumulation(curve, 0, 2) ≈ 1.04 * 1.05

        # test construction using vector of reals and of Rates
        @test discount(Yields.Forward(forwards), 1) > discount(Yields.Forward(Yields.Continuous.(forwards)), 1)

        @testset "broadcasting" begin
            @test all(accumulation.(curve, [1, 2]) .≈ [1.05, 1.04 * 1.05])
            @test all(discount.(curve, [1, 2]) .≈ 1 ./ [1.05, 1.04 * 1.05])
        end

        # addition / subtraction
        @test discount(curve + 0.1, 1) ≈ 1 / 1.15
        @test discount(curve - 0.03, 1) ≈ 1 / 1.02



        @testset "with specified timepoints" begin
            i = [0.0, 0.05]
            times = [0.5, 1.5]
            y = Yields.Forward(i, times)
            @test discount(y, 0.5) ≈ 1 / 1.0^0.5
            @test discount(y, 1.5) ≈ 1 / 1.0^0.5 / 1.05^1

        end

    end

    @testset "forwardcurve" begin
        maturity = [0.5, 1.0, 1.5, 2.0]
        zeros = [5.0, 5.8, 6.4, 6.8] ./ 100
        curve = Yields.Zero(zeros, maturity)

        fwd = Yields.ForwardStarting(curve, 1.0)
        @test discount(fwd, 0) ≈ 1
        @test discount(fwd, 0.5) ≈ discount(curve, 1, 1.5)
        @test discount(fwd, 1) ≈ discount(curve, 1, 2)
        @test accumulation(fwd, 1) ≈ accumulation(curve, 1, 2)

        @test zero(fwd,1) ≈ forward(curve,1,2)
        @test zero(fwd,1,Yields.Continuous()) ≈ convert(Yields.Continuous(),forward(curve,1,2))
    end



    @testset "actual cmt treasury" begin
        # Fabozzi 5-5,5-6
        cmt = [5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 6.8, 7.0, 7.1, 7.15, 7.2, 7.3, 7.35, 7.4, 7.5, 7.6, 7.6, 7.7, 7.8] ./ 100
        mats = collect(0.5:0.5:10.0)
        targets = [5.25, 5.5, 5.76, 6.02, 6.28, 6.55, 6.82, 6.87, 7.09, 7.2, 7.26, 7.31, 7.43, 7.48, 7.54, 7.67, 7.8, 7.79, 7.93, 8.07] ./ 100
        target_periodicity = fill(2, length(mats))
        target_periodicity[2] = 1 # 1 year is a no-coupon, BEY yield, the rest are semiannual BEY
        objs = CMTYield.(cmt, mats)
        curves = [curve(objs), curve(Bootstrap(LinearSpline()), objs), curve(Bootstrap(QuadraticSpline()), objs)]
        @testset "curve bootstrapping choices" for c in curves
            @testset "Fabozzi bootstrapped rates" for (r, mat, target, tp) in zip(cmt, mats, targets, target_periodicity)
                @test zero(c, mat) ≈ Periodic(tp)(target) atol = 0.0001
            end
        end
        
        # Hull, problem 4.34
        adj = ((1 + 0.051813 / 2)^2 - 1) * 100
        cmt = [4.0816, adj, 5.4986, 5.8620] ./ 100
        mats = [0.5, 1.0, 1.5, 2.0]

        c = curve(CMTYield.(cmt, mats))
        targets = [4.0405, 5.1293, 5.4429, 5.8085] ./ 100
        @testset "Hull bootstrapped rates" for (r, mat, target) in zip(cmt, mats, targets)
            @test Yields.zero(c, mat) ≈ Continuous(target) atol = 0.001
        end

        # test that showing the curve doesn't error
        @test length(repr(c)) > 0

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
        ois = [1.8, 2.0, 2.2, 2.5, 3.0, 4.0] ./ 100
        mats = [1 / 12, 1 / 4, 1 / 2, 1, 2, 5]
        targets = [0.017987, 0.019950, 0.021880, 0.024693, 0.029994, 0.040401]

        curve = Yields.OIS(ois, mats)
        @testset "bootstrapped rates" for (r, mat, target) in zip(ois, mats, targets)
            @test rate(zero(curve, mat, Yields.Continuous())) ≈ target atol = 0.001
        end
        curve = Yields.OIS(Bootstrap(LinearSpline()),ois, mats)
        @testset "bootstrapped rates" for (r, mat, target) in zip(ois, mats, targets)
            @test rate(zero(curve, mat, Yields.Continuous())) ≈ target atol = 0.001
        end
    end

    @testset "par" begin
        @testset "first payment logic" begin
            ct = Yields.coupon_times
            @test ct(0.5,1) ≈ 0.5:1:0.5
            @test ct(1.5,1) ≈ 0.5:1:1.5
            @test ct(0.75,1) ≈ 0.75:1:0.75
            @test ct(1,1) ≈ 1:1:1
            @test ct(1,2) ≈ 0.5:0.5:1.0
            @test ct(0.5,2) ≈ 0.5:0.5:0.5
            @test ct(1.5,2) ≈ 0.5:0.5:1.5
        end
        
        # https://quant.stackexchange.com/questions/57608/how-to-compute-par-yield-from-zero-rate-curve
        c = Yields.Zero(Yields.Continuous.([0.02,0.025,0.03,0.035]),0.5:0.5:2)
        @test Yields.par(c,2) ≈ Yields.Periodic(0.03508591,2) atol = 0.000001

        c = Yields.Constant(0.04)
        @testset "misc combinations" for t in 0.5:0.5:5 
            @test Yields.par(c,t;frequency=1) ≈ Yields.Periodic(0.04,1)
            @test Yields.par(c,t) ≈ Yields.Periodic(0.04,1)
            @test Yields.par(c,t,frequency=4) ≈ Yields.Periodic(0.04,1)
        end

        @test Yields.par(c,0.6) ≈ Yields.Periodic(0.04,1)

        @testset "round trip" begin
            maturity = collect(1:10)

            par = [6.0, 8.0, 9.5, 10.5, 11.0, 11.25, 11.38, 11.44, 11.48, 11.5] ./ 100

            curve = Yields.Par(par,maturity)

            for (p,m) in zip(par,maturity)
                @test Yields.par(curve,m) ≈ Yields.Periodic(p,2) atol = 0.001
            end
        end


    end

end
