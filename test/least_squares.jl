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

    # (Monotone convex, PCHIP and Akima curves switch formulas, so their objective has kinks: see the
    # robustness testset below.)
    @testset "$spline minimizes the squared pricing residuals" for spline in (Spline.Linear(), Spline.Quadratic(), Spline.Cubic(), Spline.MonotoneConvex())
        c0 = ZeroRateCurve(fill(0.04, length(knots)), knots, spline)
        c = fit(c0, qs, Fit.LeastSquares())
        @test isequal(fit(c0, qs), fit(c0, qs, Fit.Loss()))
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
    z = rate.(knot_rates(fit(c0, qs, Fit.LeastSquares())))

    @testset "weights" begin
        w = [1 / (q.price * maturity(q))^2 for q in qs]
        zw = rate.(knot_rates(fit(c0, qs, Fit.LeastSquares(weights = w))))
        ϕ = objective(Spline.Cubic(), w)
        @test norm2(ForwardDiff.gradient(ϕ, zw)) <= 1.0e-9 * norm2(ForwardDiff.gradient(ϕ, fill(0.04, length(knots))))
        @test maximum(abs, zw .- z) > 1.0e-4                       # the weights matter
        @test rate.(knot_rates(fit(c0, qs, Fit.LeastSquares(weights = 10 .* w)))) ≈ zw atol = 1.0e-9
        # a zero weight leaves its quote out
        keep = isodd.(eachindex(qs))
        @test rate.(knot_rates(fit(c0, qs, Fit.LeastSquares(weights = float.(keep))))) ≈ rate.(knot_rates(fit(c0, qs[keep], Fit.LeastSquares()))) atol = 1.0e-10
        @test Fit.LeastSquares(weights = w).weights === w
        @test_throws MethodError Fit.LeastSquares(w)      # keyword only
        @test Fit.LeastSquares(maxiters = 5).maxiters == 5
        @test_throws "within 1 iterations" fit(c0, qs, Fit.LeastSquares(maxiters = 1))
        @test_throws DimensionMismatch fit(c0, qs, Fit.LeastSquares(weights = w[1:10]))
        @test_throws ArgumentError fit(c0, qs, Fit.LeastSquares(weights = -w))
        @test_throws ArgumentError fit(Spline.Cubic(), qs[1:3], Fit.LeastSquares(weights = [1.0, NaN, 1.0]))
    end

    @testset "repeated maturities" begin
        cmt = CMTYield.([0.04, 0.041, 0.045, 0.046], [2, 2, 5, 5])
        @test_throws "Fit.LeastSquares" fit(Spline.Cubic(), cmt)
        @test_throws "Fit.LeastSquares" fit(Spline.Cubic(), cmt, Fit.Loss(abs2))
        c = fit(ZeroRateCurve(fill(0.04, 2), [2.0, 5.0], Spline.Linear()), cmt, Fit.LeastSquares())
        # the curve passes between each pair: one bond of a pair is cheap on it, the other rich
        r = [pv(c, q.instrument) - q.price for q in cmt]
        @test r[1] * r[2] < 0 && r[3] * r[4] < 0
    end

    @testset "knots at the quote maturities" begin
        cmt = CMTYield.([0.04, 0.042, 0.045, 0.047], [1, 2, 5, 10])
        reprice(c) = maximum(abs(pv(c, q.instrument) - q.price) / q.price for q in cmt)
        for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            exact = fit(spline, cmt, Fit.Exact())
            # an explicit exact refit matches the spline fit
            @test isequal(fit(ZeroRateCurve(fill(0.03, 4), [1.0, 2.0, 5.0, 10.0], spline), cmt, Fit.Exact()), exact)
            # least squares minimizes, and reaches the same curve, whatever the weights
            for method in (Fit.LeastSquares(), Fit.LeastSquares(weights = [1.0, 2.0, 3.0, 4.0]))
                c = fit(spline, cmt, method)
                @test reprice(c) < 1.0e-12
                @test rate.(knot_rates(c)) ≈ rate.(knot_rates(exact)) atol = 1.0e-12
            end
        end
        # Quotes that no curve with the interpolation reprices, here a quadratic spline through alternating
        # zero-coupon and par yields at 15%: the exact fit throws, and least squares returns a local minimum.
        t = [0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30]
        y = 0.15 .+ 0.01 .* (1 .- exp.(-t ./ 5))
        qq = [isodd(i) ? ZCBYield(y[i], t[i]) : CMTYield(y[i], t[i]) for i in eachindex(t)]
        @test_throws FitConvergenceError fit(Spline.Quadratic(), qq, Fit.Exact())
        c = fit(Spline.Quadratic(), qq, Fit.LeastSquares())
        ϕ(z) = sum(abs2, [pv(ZeroRateCurve(z, float.(t), Spline.Quadratic()), q.instrument) - q.price for q in qq])
        @test norm2(ForwardDiff.gradient(ϕ, rate.(knot_rates(c)))) <= 1.0e-9 * norm2(ForwardDiff.gradient(ϕ, fill(0.05, length(t))))
        @test maximum(abs(pv(c, q.instrument) - q.price) / q.price for q in qq) > 0.01
        # all-zero weights leave no quote to fit
        c4 = ZeroRateCurve(fill(0.03, 4), [1.0, 2.0, 5.0, 10.0], Spline.Linear())
        @test_throws "positive weight" fit(c4, cmt, Fit.LeastSquares(weights = zeros(4)))
        @test_throws "positive weight" fit(Spline.Linear(), cmt, Fit.LeastSquares(weights = zeros(4)))
    end

    @testset "robustness: random indices" begin
        # (nine knots each, so that the gradient checks compile once per interpolation)
        knotsets = (knots, [1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 19.0], [0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 9.0, 13.0, 20.0])
        stationary(c, q) = (
            ϕ = z -> sum(abs2, [pv(ZeroRateCurve(z, knot_tenors(c), c.spline), x.instrument) - x.price for x in q]);
            norm2(ForwardDiff.gradient(ϕ, rate.(knot_rates(c)))) <= 1.0e-9 * norm2(ForwardDiff.gradient(ϕ, fill(0.04, length(knot_tenors(c)))))
        )
        for seed in 1:20
            kn = knotsets[mod1(seed, 3)]
            q = lsq_index(80; seed)
            # every fourth index also through quotes priced one by one
            routes = seed % 4 == 0 ? (q, [Quote(x.price, WrappedContract(x.instrument)) for x in q[1:40]]) : (q,)
            for spline in (Spline.Linear(), Spline.Quadratic(), Spline.Cubic(), Spline.MonotoneConvex(), Spline.PCHIP()), qx in routes
                c = try
                    fit(ZeroRateCurve(fill(0.04, length(kn)), kn, spline), qx, Fit.LeastSquares())
                catch e
                    e
                end
                if spline isa Union{Spline.MonotoneConvex, Spline.PCHIP}
                    # They switch formulas, so a minimum can lie on a switch: such a fit throws, and says so.
                    # (A fit that converges may still end next to one, where ForwardDiff cannot check its gradient.)
                    @test !(c isa Exception) || (c isa FitConvergenceError && occursin("switches formulas", sprint(showerror, c)))
                else
                    @test !(c isa Exception) && stationary(c, qx)
                end
            end
        end
    end

    @testset "a fit is accepted on the full Hessian, not on a Gauss-Newton step" begin
        # One knot and two zero-coupon prices, chosen so that at 5% the residuals' curvature cancels JᵀJ but
        # for 1e-5, with a gradient of 1e-9: a Gauss-Newton step there is 5.6e-11, while the minimum is
        # 1.1e-6 away, where the curvature is 1.8e-3. (From a review of the fit.)
        ts = [1.0, 50.0]
        ps = [0.5704183799730305, 0.17034440723747699]
        c1 = ZeroRateCurve([0.04], [1.5], Spline.Linear())
        grad(x) = sum(ts .* exp.(-ts .* x) .* (ps .- exp.(-ts .* x)))   # of half the sum of squares
        curv(x) = sum(ts .^ 2 .* exp.(-ts .* x) .* (2 .* exp.(-ts .* x) .- ps))
        lo, hi = 0.04999, 0.05   # brackets the minimum
        for _ in 1:60
            mid = (lo + hi) / 2
            grad(mid) > 0 ? (hi = mid) : (lo = mid)
        end
        zmin = (lo + hi) / 2
        kn(p) = rate(only(knot_rates(fit(c1, ZCBPrice.(p, ts), Fit.LeastSquares()))))
        z = kn(ps)
        @test z ≈ zmin atol = 1.0e-12
        # stationary by the analytic Hessian: a minimum, whose Newton step is at rounding
        @test curv(z) > 0 && abs(grad(z) / curv(z)) <= 1.0e-12
        # the derivatives of the minimum, dz/dp = -(∂grad/∂p)/curv, and central differences of the fit. A
        # price lower by 3e-9 (the first) or 3e-10 (the second) removes this minimum, so the difference step
        # is 3e-12: its truncation and the fitted rates' rounding (1e-13) each move a difference by less
        # than 1e-4 of it.
        dz = ForwardDiff.gradient(kn, ps)
        @test dz ≈ -ts .* exp.(-ts .* zmin) ./ curv(zmin) rtol = 1.0e-6
        ε = 3.0e-12
        @test dz ≈ [(kn(ps .+ ε .* (1:2 .== i)) - kn(ps .- ε .* (1:2 .== i))) / 2ε for i in 1:2] rtol = 1.0e-3
    end

    @testset "quotes that do not determine the knots throw" begin
        @test_throws "at least as many quotes" fit(c0, qs[1:5], Fit.LeastSquares())
        cl = ZeroRateCurve(fill(0.05, 4), [1.0, 2.0, 3.0, 40.0], Spline.Linear())
        @test_throws "the knot at 40.0 does not move any quote's price" fit(cl, ZCBPrice.([0.97, 0.94, 0.9, 0.88], [1.0, 2.0, 2.5, 3.0]), Fit.LeastSquares())
        # nor a zero weight on the only quote that moves a knot
        cz = ZeroRateCurve(fill(0.05, 3), [1.0, 2.0, 3.0], Spline.Linear())
        zqs = ZCBPrice.([0.97, 0.94, 0.941, 0.9, 0.91], [1.0, 2.0, 2.0, 3.0, 3.0])
        @test_throws "the knot at 3.0" fit(cz, zqs, Fit.LeastSquares(weights = [1.0, 1.0, 1.0, 0.0, 0.0]))
    end

    @testset "quotes whose cashflows the fit cannot read are priced one by one" begin
        sub = qs[1:40]
        c5 = ZeroRateCurve(fill(0.04, 5), [1.0, 3.0, 5.0, 10.0, 20.0], Spline.Cubic())
        wrapped = [Quote(q.price, WrappedContract(q.instrument)) for q in sub]
        @test rate.(knot_rates(fit(c5, wrapped, Fit.LeastSquares()))) ≈ rate.(knot_rates(fit(c5, sub, Fit.LeastSquares()))) atol = 1.0e-10
        @test_throws "cashflows depend on the curve" ForwardDiff.derivative(
            a -> discount(fit(c5, [Quote(q.price + a, q.instrument) for q in wrapped], Fit.LeastSquares()), 4.0), 0.0
        )
    end

    @testset "derivatives" begin
        sub = qs[1:60]
        c5 = ZeroRateCurve(fill(0.04, 5), [1.0, 3.0, 5.0, 10.0, 20.0], Spline.Cubic())
        bonds = [q.instrument for q in sub]
        h = 1.0e-4
        central(f, x, i) = (e = zeros(length(x)); e[i] = h; (f(x .+ e) .- f(x .- e)) ./ 2h)
        # quote prices
        kp(p) = rate.(knot_rates(fit(c5, Quote.(p, bonds), Fit.LeastSquares())))
        prices = [q.price for q in sub]
        Jp = ForwardDiff.jacobian(kp, prices)
        @test Jp[:, 1] ≈ central(kp, prices, 1) rtol = 1.0e-6
        @test Jp[:, 37] ≈ central(kp, prices, 37) rtol = 1.0e-6
        # coupons (cashflow amounts)
        kc(x) = rate.(knot_rates(fit(c5, [Quote(q.price, Bond.Fixed(x[i], Periodic(2), maturity(q))) for (i, q) in enumerate(sub)], Fit.LeastSquares())))
        coupons = [b.coupon_rate for b in bonds]
        @test ForwardDiff.jacobian(kc, coupons)[:, 5] ≈ central(kc, coupons, 5) rtol = 1.0e-6
        # weights
        w = [1 / (q.price * maturity(q))^2 for q in sub]
        kw(x) = rate.(knot_rates(fit(c5, sub, Fit.LeastSquares(weights = x))))
        Jw = ForwardDiff.jacobian(kw, w)
        @test Jw[:, 11] ≈ central(kw, w, 11) rtol = 1.0e-5
        # weights in a vector of abstract element type carry their derivatives too
        @test ForwardDiff.jacobian(x -> kw(Real[x...]), w) ≈ Jw rtol = 1.0e-14
        # quoted yields, through a valuation
        tenors = [1, 2, 3, 5, 7, 10]
        rates = [0.03, 0.032, 0.034, 0.036, 0.037, 0.038]
        cfs = Cashflow.([3.0, 3.0, 103.0], [2.0, 5.0, 10.0])
        c3 = ZeroRateCurve(fill(0.03, 3), [1.0, 4.0, 10.0], Spline.Linear())
        value(x) = pv(fit(c3, CMTYield.(x, tenors), Fit.LeastSquares()), cfs)
        g = ForwardDiff.gradient(value, rates)
        @test g ≈ [(value(rates .+ h .* (1:6 .== i)) - value(rates .- h .* (1:6 .== i))) / 2h for i in 1:6] rtol = 1.0e-6
        # a loss fit with another number of knots is not differentiated
        @test_throws "Fit.LeastSquares()" ForwardDiff.gradient(x -> pv(fit(c3, CMTYield.(x, tenors), Fit.Loss(abs2)), cfs), rates)
        # nor is a least-squares fit's extrapolation forward
        ce(f) = ZeroRateCurve(fill(0.03, 3), [1.0, 4.0, 10.0], Spline.Linear(); extrapolation = Yield.FlatForwardAt(Continuous(f)))
        @test_throws "extrapolation forward" ForwardDiff.derivative(f -> discount(fit(ce(f), CMTYield.(rates, tenors), Fit.LeastSquares()), 12.0), 0.04)
    end

    @testset "FX forwards fit their foreign knot curve by least squares" begin
        eurusd = FX.Pair(:EUR, :USD)
        usd = Yield.Constant(Continuous(0.05))
        truth = FX.Forwards(eurusd, 1.1, usd, Yield.Constant(Continuous(0.03)))
        outrights = [FX.Outright(eurusd, forward(truth, t), t) for t in 0.5:0.5:6.0]
        foreign = ZeroRateCurve(fill(0.02, 3), [1.0, 3.0, 6.0], Spline.Linear())
        m = fit(FX.Forwards(eurusd, 1.1, usd, foreign), outrights, Fit.LeastSquares())
        @test knot_tenors(m.foreign) == [1.0, 3.0, 6.0]
        @test all(abs(pv(m, q.instrument)) < 1.0e-6 for q in outrights)   # a flat curve is representable
        @test isequal(fit(FX.Forwards(eurusd, 1.1, usd, foreign), outrights).foreign, fit(FX.Forwards(eurusd, 1.1, usd, foreign), outrights, Fit.Loss()).foreign)
    end
end
