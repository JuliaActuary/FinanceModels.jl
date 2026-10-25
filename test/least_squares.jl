using ForwardDiff, Random

# A bond index: semiannual bonds on a monthly maturity grid, priced off a smooth curve plus a spread per
# issuer, so that no curve reprices them all.
function lsq_index(nb; seed = 11)
    rng = Random.Xoshiro(seed)
    truth = Yield.NelsonSiegel(2.0, 0.058, -0.02, 0.01)
    mats = round.(12 .* (0.5 .+ 19.5 .* rand(rng, nb))) ./ 12
    cpns = round.(800 .* (0.02 .+ 0.05 .* rand(rng, nb))) ./ 800
    spreads = 0.003 .* randn(rng, nb)
    bonds = [Bond.Fixed(c, Periodic(2), m) for (c, m) in zip(cpns, mats)]
    return [Quote(pv(truth + Yield.Constant(Continuous(s)), b), b) for (b, s) in zip(bonds, spreads)]
end

@testset "least-squares fits of knot curves" begin
    qs = lsq_index(200)
    knots = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
    norm2(v) = sqrt(sum(abs2, v))
    objective(spline, w = ones(length(qs))) = z -> sum(w[i] * (pv(ZeroRateCurve(z, knots, spline), q.instrument) - q.price)^2 for (i, q) in enumerate(qs))

    # (PCHIP and Akima fits also converge, but their objective has kinks where their formula switches, so
    # its gradient is not checked here.)
    @testset "$spline minimizes the squared pricing residuals" for spline in (Spline.Linear(), Spline.Quadratic(), Spline.Cubic(), Spline.MonotoneConvex())
        c0 = ZeroRateCurve(fill(0.04, length(knots)), knots, spline)
        c = fit(c0, qs)
        @test c == fit(c0, qs, Fit.LeastSquares())
        @test knot_tenors(c) == knots
        z = rate.(knot_rates(c))
        ϕ = objective(spline)
        @test norm2(ForwardDiff.gradient(ϕ, z)) <= 1.0e-9 * norm2(ForwardDiff.gradient(ϕ, fill(0.04, length(knots))))
        # the minimizer of the squared loss (which Optim's default for MonotoneConvex, LBFGS, does not reach
        # to this tolerance)
        if !(spline isa Spline.MonotoneConvex)
            zl = rate.(knot_rates(fit(c0, qs, Fit.Loss(abs2); solve_kwargs = (; g_tol = 1.0e-12))))
            @test z ≈ zl atol = 1.0e-9
        end
    end

    c0 = ZeroRateCurve(fill(0.04, length(knots)), knots, Spline.Cubic())
    z = rate.(knot_rates(fit(c0, qs)))

    @testset "weights" begin
        w = [1 / (q.price * maturity(q))^2 for q in qs]
        zw = rate.(knot_rates(fit(c0, qs, Fit.LeastSquares(weights = w))))
        ϕ = objective(Spline.Cubic(), w)
        @test norm2(ForwardDiff.gradient(ϕ, zw)) <= 1.0e-9 * norm2(ForwardDiff.gradient(ϕ, fill(0.04, length(knots))))
        @test maximum(abs, zw .- z) > 1.0e-4                       # the weights matter
        @test rate.(knot_rates(fit(c0, qs, Fit.LeastSquares(weights = 10 .* w)))) ≈ zw atol = 1.0e-9
        # a zero weight leaves its quote out
        keep = isodd.(eachindex(qs))
        @test rate.(knot_rates(fit(c0, qs, Fit.LeastSquares(weights = float.(keep))))) ≈ rate.(knot_rates(fit(c0, qs[keep]))) atol = 1.0e-10
        @test_throws DimensionMismatch fit(c0, qs, Fit.LeastSquares(weights = w[1:10]))
        @test_throws ArgumentError fit(c0, qs, Fit.LeastSquares(weights = -w))
        @test_throws ArgumentError fit(Spline.Cubic(), qs[1:3], Fit.LeastSquares(weights = [1.0, NaN, 1.0]))
    end

    @testset "repeated maturities" begin
        cmt = CMTYield.([0.04, 0.041, 0.045, 0.046], [2, 2, 5, 5])
        @test_throws "Fit.LeastSquares" fit(Spline.Cubic(), cmt)
        @test_throws "Fit.LeastSquares" fit(Spline.Cubic(), cmt, Fit.Loss(abs2))
        c = fit(ZeroRateCurve(fill(0.04, 2), [2.0, 5.0], Spline.Linear()), cmt)
        # the curve passes between each pair: one bond of a pair is cheap on it, the other rich
        r = [pv(c, q.instrument) - q.price for q in cmt]
        @test r[1] * r[2] < 0 && r[3] * r[4] < 0
    end

    @testset "knots at the quote maturities: the exact fit" begin
        cmt = CMTYield.([0.04, 0.042, 0.045, 0.047], [1, 2, 5, 10])
        for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            @test isequal(fit(spline, cmt, Fit.LeastSquares()), fit(spline, cmt))
            @test isequal(fit(spline, cmt, Fit.LeastSquares(weights = [1.0, 2.0, 3.0, 4.0])), fit(spline, cmt))
            @test isequal(fit(ZeroRateCurve(fill(0.03, 4), [1.0, 2.0, 5.0, 10.0], spline), cmt), fit(spline, cmt))
        end
    end

    @testset "quotes that do not determine the knots throw" begin
        @test_throws "at least as many quotes" fit(c0, qs[1:5])
        cl = ZeroRateCurve(fill(0.05, 4), [1.0, 2.0, 3.0, 40.0], Spline.Linear())
        @test_throws "the knot at 40.0 does not move any quote's price" fit(cl, ZCBPrice.([0.97, 0.94, 0.9, 0.88], [1.0, 2.0, 2.5, 3.0]))
        # nor a zero weight on the only quote that moves a knot
        cz = ZeroRateCurve(fill(0.05, 3), [1.0, 2.0, 3.0], Spline.Linear())
        zqs = ZCBPrice.([0.97, 0.94, 0.9, 0.91], [1.0, 2.0, 3.0, 3.0])
        @test_throws "the knot at 3.0" fit(cz, zqs, Fit.LeastSquares(weights = [1.0, 1.0, 0.0, 0.0]))
    end

    @testset "quotes whose cashflows the fit cannot read are priced one by one" begin
        sub = qs[1:40]
        c5 = ZeroRateCurve(fill(0.04, 5), [1.0, 3.0, 5.0, 10.0, 20.0], Spline.Cubic())
        wrapped = [Quote(q.price, WrappedContract(q.instrument)) for q in sub]
        @test rate.(knot_rates(fit(c5, wrapped))) ≈ rate.(knot_rates(fit(c5, sub))) atol = 1.0e-10
        @test_throws "cashflows depend on the curve" ForwardDiff.derivative(
            a -> discount(fit(c5, [Quote(q.price + a, q.instrument) for q in wrapped]), 4.0), 0.0
        )
    end

    @testset "derivatives" begin
        sub = qs[1:60]
        c5 = ZeroRateCurve(fill(0.04, 5), [1.0, 3.0, 5.0, 10.0, 20.0], Spline.Cubic())
        bonds = [q.instrument for q in sub]
        h = 1.0e-4
        central(f, x, i) = (e = zeros(length(x)); e[i] = h; (f(x .+ e) .- f(x .- e)) ./ 2h)
        # quote prices
        kp(p) = rate.(knot_rates(fit(c5, Quote.(p, bonds))))
        prices = [q.price for q in sub]
        Jp = ForwardDiff.jacobian(kp, prices)
        @test Jp[:, 1] ≈ central(kp, prices, 1) rtol = 1.0e-6
        @test Jp[:, 37] ≈ central(kp, prices, 37) rtol = 1.0e-6
        # coupons (cashflow amounts)
        kc(x) = rate.(knot_rates(fit(c5, [Quote(q.price, Bond.Fixed(x[i], Periodic(2), maturity(q))) for (i, q) in enumerate(sub)])))
        coupons = [b.coupon_rate for b in bonds]
        @test ForwardDiff.jacobian(kc, coupons)[:, 5] ≈ central(kc, coupons, 5) rtol = 1.0e-6
        # weights
        w = [1 / (q.price * maturity(q))^2 for q in sub]
        kw(x) = rate.(knot_rates(fit(c5, sub, Fit.LeastSquares(weights = x))))
        @test ForwardDiff.jacobian(kw, w)[:, 11] ≈ central(kw, w, 11) rtol = 1.0e-5
        # quoted yields, through a valuation
        tenors = [1, 2, 3, 5, 7, 10]
        rates = [0.03, 0.032, 0.034, 0.036, 0.037, 0.038]
        cfs = Cashflow.([3.0, 3.0, 103.0], [2.0, 5.0, 10.0])
        c3 = ZeroRateCurve(fill(0.03, 3), [1.0, 4.0, 10.0], Spline.Linear())
        value(x) = pv(fit(c3, CMTYield.(x, tenors)), cfs)
        g = ForwardDiff.gradient(value, rates)
        @test g ≈ [(value(rates .+ h .* (1:6 .== i)) - value(rates .- h .* (1:6 .== i))) / 2h for i in 1:6] rtol = 1.0e-6
        # a loss fit with another number of knots is not differentiated
        @test_throws "Fit.LeastSquares()" ForwardDiff.gradient(x -> pv(fit(c3, CMTYield.(x, tenors), Fit.Loss(abs2)), cfs), rates)
        # nor is a least-squares fit's extrapolation forward
        ce(f) = ZeroRateCurve(fill(0.03, 3), [1.0, 4.0, 10.0], Spline.Linear(); extrapolation = Yield.FlatForwardAt(Continuous(f)))
        @test_throws "extrapolation forward" ForwardDiff.derivative(f -> discount(fit(ce(f), CMTYield.(rates, tenors)), 12.0), 0.04)
    end

    @testset "FX forwards fit their foreign knot curve by least squares" begin
        eurusd = FX.Pair(:EUR, :USD)
        usd = Yield.Constant(Continuous(0.05))
        truth = FX.Forwards(eurusd, 1.1, usd, Yield.Constant(Continuous(0.03)))
        outrights = [FX.Outright(eurusd, forward(truth, t), t) for t in 0.5:0.5:6.0]
        foreign = ZeroRateCurve(fill(0.02, 3), [1.0, 3.0, 6.0], Spline.Linear())
        m = fit(FX.Forwards(eurusd, 1.1, usd, foreign), outrights)
        @test knot_tenors(m.foreign) == [1.0, 3.0, 6.0]
        @test all(abs(pv(m, q.instrument)) < 1.0e-6 for q in outrights)   # a flat curve is representable
        @test knot_rates(m.foreign) == knot_rates(fit(FX.Forwards(eurusd, 1.1, usd, foreign), outrights, Fit.LeastSquares()).foreign)
    end
end
