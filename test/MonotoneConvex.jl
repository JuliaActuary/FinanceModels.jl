@testset "Monotone Convex" begin

    # Interpolation of the yield curve, Gustaf Dehlbom
    # http://uu.diva-portal.org/smash/get/diva2:1477828/FULLTEXT01.pdf

    prices = [0.98, 0.955, 0.92, 0.88, 0.83]
    times = [1, 2, 3, 4, 5]
    quotes = ZCBPrice.(prices, times)
    rates = @. -log(prices) / times

    curves = [
        ("direct constructor", Yield.MonotoneConvex(rates, times)),
        ("fit", fit(Yield.MonotoneConvex(), quotes)),
    ]

    f, fᵈ = Yield.__monotone_convex_fs(rates, times)

    @testset "$name" for (name, c) in curves
        @test all(isapprox.(f, c.f; atol=1e-8))
        @test all(isapprox.(fᵈ, c.fᵈ; atol=1e-8))

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
        @test f[6] ≈ 0.062 atol = 0.0001

        @test FinanceModels.Yield.g(0.5, 0.018793076350927487, 0.023021969250703423, 0.020202707317519466) ≈ 0.0042 * (0.5)^2 - 0.0014 atol = 0.0001
        @test FinanceModels.Yield.g(0, 0.023021969250703423, 0.03158945081076577, 0.02584123118388738) ≈ -0.0028 atol = 0.0001
        @test FinanceModels.Yield.g(0.5, 0.023021969250703423, 0.03158945081076577, 0.02584123118388738) ≈ 0.0087 * (0.5)^2 − 0.0002 * 0.5 - 0.0028 atol = 0.0001
        @test FinanceModels.Yield.g(0.5, 0.03158945081076577, 0.04089471650423902, 0.03733767043764417) ≈ -0.0063 * (0.5)^2 + 0.0156 * 0.5 - 0.0057 atol = 0.0001
        @test FinanceModels.Yield.g(0.5, 0.04089471650423902, 0.051473984626221235, 0.04445176257083387) ≈ 0.0102 * (0.5)^2 + 0.0004 * 0.5 - 0.0036 atol = 0.0001
    end

    # Reference function for zero rates (tolerance needed due to rounded coefficients)
    function r(t)
        if 0 <= t <= 1
            return 0.0014t^2 + 0.0188
        elseif 1 <= t <= 1.0233
            return -0.0028 / t + 0.023
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

    @testset "zero rates, knots, discounts - $name" for (name, c) in curves
        # Test zero rates against reference function
        @testset "zero rate at t=$t" for t in range(0, 5, 30)
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
    end


    # Tests from Google TF Quant Finance
    # https://github.com/google/tf-quant-finance/blob/master/tf_quant_finance/rates/hagan_west/monotone_convex_test.py

    @testset "TF Quant Finance - interpolated yields" begin
        # test_interpolated_yields_with_yields Example1
        # Reference times and yields (zero rates in percent)
        ref_times = [1.0, 2.0, 3.0, 4.0]
        yields = [5.0, 4.75, 4.533333, 4.775] ./ 100  # convert to decimal

        c = Yield.MonotoneConvex(yields, ref_times)

        # Interpolation times and expected yields
        interp_times = [0.25, 0.5, 1.0, 2.0, 3.0, 1.1, 2.5, 2.9, 3.6, 4.0]
        expected = [5.1171875, 5.09375, 5.0, 4.75, 4.533333, 4.9746, 4.624082, 4.535422, 4.661777, 4.775] ./ 100

        @testset "t=$t" for (i, t) in enumerate(interp_times)
            @test rate(zero(c, t)) ≈ expected[i] atol = 0.0001
        end
    end

    @testset "TF Quant Finance - yields at knot points" begin
        # test_interpolated_yields_with_yields Example2
        # Verifies exact interpolation at knot points
        ref_times = [0.1, 0.2, 0.21]
        yields = [0.1, 0.2, 0.3]

        c = Yield.MonotoneConvex(yields, ref_times)

        @testset "t=$t" for (i, t) in enumerate(ref_times)
            @test rate(zero(c, t)) ≈ yields[i] atol = 1.0e-10
        end
    end

    @testset "fit with ZCBPrice" begin
        prices = [0.98, 0.955, 0.92, 0.88, 0.83]
        times = [1, 2, 3, 4, 5]
        quotes = ZCBPrice.(prices, times)

        c = fit(Yield.MonotoneConvex(), quotes)

        # Verify discount factors match prices
        @testset "discount at t=$t" for (i, t) in enumerate(times)
            @test discount(c, t) ≈ prices[i] atol = 1.0e-6
        end
    end

    @testset "fit with ParYield" begin
        par_rates = [0.02, 0.025, 0.03, 0.035, 0.04]
        times = [1, 2, 3, 4, 5]
        quotes = ParYield.(par_rates, times)

        c = fit(Yield.MonotoneConvex(), quotes)

        # Verify par rates match
        @testset "par at t=$t" for (i, t) in enumerate(times)
            @test rate(par(c, t)) ≈ par_rates[i] atol = 0.0001
        end
    end

    @testset "forward rates" begin
        prices = [0.98, 0.955, 0.92, 0.88, 0.83]
        times = [1, 2, 3, 4, 5]
        rates = @. -log(prices) / times
        quotes = ZCBPrice.(prices, times)

        curves = [
            ("direct constructor", Yield.MonotoneConvex(rates, times)),
            ("fit", fit(Yield.MonotoneConvex(), quotes)),
        ]

        @testset "$name" for (name, c) in curves
            # Forward at t=0 should equal instantaneous forward f[1]
            @test Yield.instantaneous_forward(c, 0.0) ≈ c.f[1] atol = 1.0e-10

            # Forward at internal knot points equals the instantaneous forward f[i+1]
            # (at x=1 of the interval ending at that knot)
            @testset "forward at knot t=$t" for (i, t) in enumerate(times[1:(end - 1)])
                @test Yield.instantaneous_forward(c, t) ≈ c.f[i + 1] atol = 1.0e-10
            end

            # Forward at/beyond the last knot equals the boundary instantaneous
            # forward f[end], keeping the forward curve continuous at t_n
            @test Yield.instantaneous_forward(c, times[end]) ≈ c.f[end] atol = 1.0e-10

            # Forward should be continuous and positive everywhere, including across the last knot
            @testset "forward positive and continuous" begin
                ts = range(0.01, 6.0, 120)
                fwds = [Yield.instantaneous_forward(c, t) for t in ts]
                @test all(fwds .> 0)

                # Check continuity: adjacent forward values shouldn't jump excessively
                for i in 2:length(fwds)
                    @test abs(fwds[i] - fwds[i - 1]) < 0.01
                end
            end

            # Extrapolation: flat at the boundary instantaneous forward
            @test Yield.instantaneous_forward(c, 6.0) ≈ c.f[end] atol = 1.0e-10
            @test Yield.instantaneous_forward(c, 10.0) ≈ c.f[end] atol = 1.0e-10

            # Continuity exactly at the last knot (this used to jump from f[end] to fᵈ[end])
            @testset "continuity at the last knot" begin
                ε = 1.0e-9
                @test Yield.instantaneous_forward(c, times[end] - ε) ≈
                    Yield.instantaneous_forward(c, times[end] + ε) atol = 1.0e-6
                @test rate(zero(c, times[end] + 1.0e-6)) ≈ rate(zero(c, times[end])) atol = 1.0e-6
            end

            # Mid-interval forwards should be bounded by fᵈ ± max deviation
            @testset "mid-interval forward bounds" begin
                # First interval: t ∈ (0, 1)
                f_mid = Yield.instantaneous_forward(c, 0.5)
                @test f_mid > 0  # Must be positive
                @test abs(f_mid - c.fᵈ[1]) < max(abs(c.f[1] - c.fᵈ[1]), abs(c.f[2] - c.fᵈ[1])) + 0.001

                # Second interval: t ∈ (1, 2)
                f_mid = Yield.instantaneous_forward(c, 1.5)
                @test f_mid > 0  # Must be positive
                @test abs(f_mid - c.fᵈ[2]) < max(abs(c.f[2] - c.fᵈ[2]), abs(c.f[3] - c.fᵈ[2])) + 0.001
            end
        end
    end

    @testset "positivity collar" begin
        times = [1.0, 2.0, 3.0, 4.0]
        # discrete forwards with sharp spikes/dips so the Hagan-West collar must bind
        cases = [
            [0.001, 0.20, 0.05, 0.04],
            [0.10, 0.001, 0.10, 0.002],
            [0.05, 0.001, 0.15, 0.001],
            [0.2, 0.02, 0.001, 0.1],
        ]
        @testset "fᵈ = $fd" for fd in cases
            # zero rates whose discrete forwards equal `fd`: r_i = (Σ_{j≤i} fdⱼ·Δtⱼ)/t_i
            rates = cumsum(fd) ./ times
            f, fᵈ = Yield.__monotone_convex_fs(rates, times)
            @test fᵈ ≈ fd

            # node forwards are collared against their *adjacent* discrete forwards:
            # boundary nodes have one neighbor; interior node t_j adjoins intervals j and j+1
            @test 0 <= f[1] <= 2 * fᵈ[1] + eps()
            @test 0 <= f[end] <= 2 * fᵈ[end] + eps()
            for j in 1:(length(times) - 1)
                @test 0 <= f[j + 1] <= 2 * min(fᵈ[j], fᵈ[j + 1]) + eps()
            end

            c = Yield.MonotoneConvex(rates, times)
            # the collared interpolant keeps instantaneous forwards nonnegative...
            @test minimum(Yield.instantaneous_forward(c, t) for t in range(1.0e-6, 4.0, 1001)) >= -1.0e-12
            # ...without giving up exact knot repricing (the g construction always
            # integrates to the discrete forwards regardless of node values)
            for (i, t) in enumerate(times)
                @test rate(zero(c, t)) ≈ rates[i] atol = 1.0e-12
            end
        end
    end

    @testset "negative rates" begin
        rates = [-0.005, -0.003, 0.001, 0.004]
        times = [1.0, 2.0, 3.0, 5.0]
        c = Yield.MonotoneConvex(rates, times)
        for (i, t) in enumerate(times)
            @test rate(zero(c, t)) ≈ rates[i] atol = 1.0e-12
        end
        # negative rates imply discount factors above 1 and finite forwards everywhere
        @test discount(c, 1.0) > 1.0
        @test all(isfinite(Yield.instantaneous_forward(c, t)) for t in range(0.01, 6.0, 200))
    end

    @testset "single knot" begin
        c = Yield.MonotoneConvex([0.03], [2.0])
        @test rate(zero(c, 2.0)) ≈ 0.03
        @test rate(zero(c, 1.0)) ≈ 0.03  # flat forward within the single interval
        @test Yield.instantaneous_forward(c, 0.5) ≈ 0.03
        @test rate(zero(c, 4.0)) ≈ 0.03  # constant-forward extrapolation
    end

    @testset "equal adjacent discrete forwards (g1 == 0 / η == 0)" begin
        # When two adjacent discrete forwards coincide (any flat-forward segment),
        # the interpolated node forward equals the discrete forward, so the
        # sector-(iii) deviation g1 is exactly 0 and η == 0. The zero-rate integral
        # G is evaluated at x == 0 exactly at each knot; the unguarded 0/0 there
        # made both the discount factor AND its ForwardDiff gradient NaN — the
        # latter silently stalled `fit` at its initial guess for every optimizer.

        # a flat zero-rate curve makes every adjacent pair of discrete forwards equal
        times = [1 / 12, 2 / 12, 3 / 12, 0.5, 1.0, 2.0, 3.0, 5.0]
        cflat = Yield.MonotoneConvex(fill(0.03, length(times)), times)
        @test all(isfinite(discount(cflat, t)) for t in range(1.0e-6, last(times); length = 200))
        for (i, t) in enumerate(times)
            @test rate(zero(cflat, t)) ≈ 0.03 atol = 1.0e-12
        end

        # the reported scenario: par quotes including short (≤6m) tenors. Before the
        # fix `fit` returned a NaN-repricing curve for every optimizer because the
        # loss gradient was NaN at the (degenerate) fixed ramp seed.
        tenors = [1 / 12, 2 / 12, 3 / 12, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0]
        pars = [0.0375, 0.0372, 0.0372, 0.0372, 0.0364, 0.0368,
            0.0369, 0.0380, 0.0400, 0.0423, 0.0483, 0.0486]
        qs = sort(CMTYield.(pars, tenors); by = maturity)
        reprice(c) = maximum(abs, present_value(c, q.instrument) - q.price for q in qs)
        for opt in (
                FinanceModels.OptimizationOptimJL.LBFGS(),
                FinanceModels.OptimizationOptimJL.BFGS(),
                FinanceModels.OptimizationOptimJL.Newton(),
            )
            c = fit(Yield.MonotoneConvex(), qs; optimizer = opt)
            @test reprice(c) < 1.0e-6   # finite (not NaN) and actually fits
        end

        # the Spline.MonotoneConvex tag routes to the native curve on every path
        @test reprice(fit(Spline.MonotoneConvex(), qs)) < 1.0e-6
        @test reprice(fit(Spline.MonotoneConvex(), qs, Fit.Loss(x -> x^2))) < 1.0e-6
        @test reprice(fit(Spline.MonotoneConvex(), qs, Fit.Bootstrap())) < 1.0e-6
    end

    @testset "g0 == 0 (flat-forward segment): continuous zero rates" begin
        # The mirror of the g1 == 0 case: a flat-forward segment makes the node
        # forward equal the discrete forward at the segment's right boundary, so
        # g0 == 0 there. That used to route to sector (iii) with η == 3 > 1, where
        # ∫₀¹ g ≠ 0, so the interval failed to reprice its right knot — a small
        # (~1e-3) zero-rate discontinuity (visible approaching the knot from the
        # left, since the knot itself is repriced via the next interval).
        rates = [0.02, 0.03, 0.035, 0.035, 0.035, 0.04, 0.05]
        times = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0]
        c = Yield.MonotoneConvex(rates, times)
        for (i, t) in enumerate(times)
            @test rate(zero(c, t)) ≈ rates[i] atol = 1.0e-12          # reprices at the knot
            # z is continuous, so the left limit must also approach rate[i]; the
            # sector-(iii) η>1 misroute made this off by ~9e-4 at the g0==0 knot
            i > 1 && @test rate(zero(c, t - 1.0e-6)) ≈ rates[i] atol = 1.0e-4
        end

        # g0 == 0 in the LAST interval forces the x == 1 boundary evaluation at the
        # last knot (i_time is capped at lastindex), exercising the x == 1 guard
        # that keeps sector (ii)'s η == 1 from dividing by zero.
        clast = Yield.MonotoneConvex([0.03, 0.035, 0.035], [1.0, 2.0, 3.0])
        @test isfinite(rate(zero(clast, 3.0)))
        @test rate(zero(clast, 3.0)) ≈ 0.035 atol = 1.0e-12
    end

    @testset "Accessors / generic fit compatibility" begin
        times = [1.0, 2.0, 3.0, 5.0, 7.0]
        rates = [0.03, 0.032, 0.034, 0.037, 0.039]
        c = Yield.MonotoneConvex(rates, times)

        # one bounded optic per knot (previously a broadcast over a ClosedInterval
        # that threw `iterate(::ClosedInterval)` the moment it was evaluated)
        optic = FinanceModels.__default_optic(c)
        @test length(optic) == length(rates)

        # reconstructable via Accessors, and the cached f/fᵈ stay consistent with
        # the updated rates (the struct caches derived fields but had no 4-arg ctor)
        c2 = Accessors.@set c.rates[2] = 0.05
        @test c2.rates[2] == 0.05
        f, fᵈ = Yield.__monotone_convex_fs(c2.rates, c2.times)
        @test c2.f == f
        @test c2.fᵈ == fᵈ

        # the generic, variables-driven `fit(model, quotes)` path now runs (it used
        # the broken optic and the non-reconstructable struct) and reprices
        qs = CMTYield.(rates, times)
        fitted = fit(c, qs)
        @test maximum(abs, present_value(fitted, q.instrument) - q.price for q in qs) < 1.0e-6
    end

    @testset "second-order optimizer fit is quiet and accurate" begin
        # A second-order optimizer (e.g. `Newton`) needs a Hessian. The loss function
        # declares `SecondOrder` AD so OptimizationBase finds it directly instead of
        # warning and auto-promoting `AutoForwardDiff` → `SecondOrder(...)` once per
        # fit. The bare `@test_logs` (no patterns, default `min_level = Warn`) asserts
        # the fit emits no warnings. The default `LBFGS` path is unaffected: it only
        # requests gradients, still computed via the inner `AutoForwardDiff`.
        qs = CMTYield.([0.04, 0.045, 0.05, 0.052], [1.0, 5.0, 10.0, 30.0])
        reprice(c) = maximum(abs, present_value(c, q.instrument) - q.price for q in qs)
        N = FinanceModels.OptimizationOptimJL.Newton()

        c = @test_logs fit(Yield.MonotoneConvex(), qs; optimizer = N)
        @test reprice(c) < 1.0e-10

        # the `Spline.MonotoneConvex` tag routes through the same loss function
        cs = @test_logs fit(Spline.MonotoneConvex(), qs, Fit.Loss(x -> x^2); optimizer = N)
        @test reprice(cs) < 1.0e-10
    end
end
