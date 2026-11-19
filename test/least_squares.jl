using ForwardDiff, Random, LinearAlgebra
using Accessors, IntervalSets

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

# A contract type that `fit` knows nothing about.
struct LsqOpaque{C} <: FinanceCore.AbstractContract
    inner::C
end
FinanceCore.maturity(c::LsqOpaque) = FinanceCore.maturity(c.inner)
FinanceCore.present_value(m, c::LsqOpaque) = FinanceCore.present_value(m, c.inner)

# A contract priced through a least-squares fit of its own: the 5-year discount factor of a two-knot curve
# fitted to `quotes` with their prices scaled by the model's 1-year discount factor.
struct LsqRefit{Q} <: FinanceCore.AbstractContract
    quotes::Q
end
function FinanceCore.present_value(m, c::LsqRefit)
    scaled = [Quote(q.price * discount(m, 1.0), q.instrument) for q in c.quotes]
    return discount(fit(ZeroRateCurve(fill(0.04, 2), [2.0, 10.0], Spline.Linear()), scaled, Fit.LeastSquares()), 5.0)
end

@testset "least-squares fits" begin
    qs = lsq_index(200)
    knots = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
    zr(c) = rate.(knot_rates(c))
    # half the weighted sum of squares, written independently of the fit
    objective(spline, q = qs, kn = knots, w = ones(length(q))) =
        z -> sum(w[i] * (pv(ZeroRateCurve(z, kn, spline), x.instrument) - x.price)^2 for (i, x) in enumerate(q)) / 2
    # a minimum by ForwardDiff's own Hessian: positive definite, and a Newton correction at rounding
    function certified(ϕ, z; atol = 1.0e-10)
        H = ForwardDiff.hessian(ϕ, z)
        return isposdef(Symmetric(H)) && maximum(abs, H \ ForwardDiff.gradient(ϕ, z)) <= atol
    end

    @testset "$spline: a certified minimum of the squared pricing residuals" for spline in (Spline.Linear(), Spline.Quadratic(), Spline.Cubic())
        c0 = ZeroRateCurve(fill(0.04, length(knots)), knots, spline)
        c = fit(c0, qs, Fit.LeastSquares())
        @test knot_tenors(c) == knots
        @test c.extrapolation == c0.extrapolation
        @test certified(objective(spline), zr(c))
        # the minimizer that a tightened loss fit reaches
        @test zr(c) ≈ zr(fit(c0, qs, Fit.Loss(abs2); solve_kwargs = (; g_tol = 1.0e-12))) atol = 1.0e-9
    end

    c0 = ZeroRateCurve(fill(0.04, length(knots)), knots, Spline.Cubic())
    z = zr(fit(c0, qs, Fit.LeastSquares()))

    @testset "the default is unchanged" begin
        @test isequal(fit(c0, qs), fit(c0, qs, Fit.Loss(x -> x^2)))
        @test !isequal(zr(fit(c0, qs)), z)   # (a loss fit stops at its gradient tolerance)
    end

    @testset "weights" begin
        w = [1 / (q.price * maturity(q))^2 for q in qs]
        zw = zr(fit(c0, qs, Fit.LeastSquares(weights = w)))
        @test certified(objective(Spline.Cubic(), qs, knots, w ./ sum(w)), zw)
        @test maximum(abs, zw .- z) > 1.0e-4                       # the weights matter
        @test zr(fit(c0, qs, Fit.LeastSquares(weights = 10 .* w))) ≈ zw atol = 1.0e-10
        # a zero weight leaves its quote out
        keep = isodd.(eachindex(qs))
        @test zr(fit(c0, qs, Fit.LeastSquares(weights = float.(keep)))) ≈ zr(fit(c0, qs[keep], Fit.LeastSquares())) atol = 1.0e-10
        @test Fit.LeastSquares(weights = w).weights === w
        @test_throws MethodError Fit.LeastSquares(w)      # keyword only
        @test Fit.LeastSquares(maxiters = 5).maxiters == 5
        @test_throws "within 1 iterations" fit(c0, qs, Fit.LeastSquares(maxiters = 1))
        @test_throws DimensionMismatch fit(c0, qs, Fit.LeastSquares(weights = w[1:10]))
        @test_throws DomainError fit(c0, qs, Fit.LeastSquares(weights = -w))
        @test_throws "not finite" fit(c0, qs, Fit.LeastSquares(weights = fill(NaN, length(qs))))
        # the loss fit's solver options are not least-squares options
        @test_throws MethodError fit(c0, qs, Fit.LeastSquares(); solve_kwargs = (; maxiters = 10))
    end

    @testset "repeated maturities" begin
        cmt = CMTYield.([0.04, 0.041, 0.045, 0.046], [2, 2, 5, 5])
        @test_throws "strictly increasing" fit(Spline.Cubic(), cmt, Fit.LeastSquares())   # one knot per quote maturity
        c = fit(ZeroRateCurve(fill(0.04, 2), [2.0, 5.0], Spline.Linear()), cmt, Fit.LeastSquares())
        # the curve passes between each pair: one bond of a pair is cheap on it, the other rich
        r = [pv(c, q.instrument) - q.price for q in cmt]
        @test r[1] * r[2] < 0 && r[3] * r[4] < 0
    end

    @testset "knots at the quote maturities: $spline reprices" for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex(), Spline.PCHIP(), Spline.Akima())
        cmt = CMTYield.([0.04, 0.042, 0.045, 0.047], [1, 2, 5, 10])
        loss = fit(spline, cmt, Fit.Loss(abs2); solve_kwargs = (; g_tol = 1.0e-14))
        # the same curve whatever the weights, and the one an existing curve with those knots refits to
        for method in (Fit.LeastSquares(), Fit.LeastSquares(weights = [1.0, 2.0, 3.0, 4.0]))
            c = fit(spline, cmt, method)
            @test maximum(abs(pv(c, q.instrument) - q.price) for q in cmt) < 1.0e-14
            @test zr(c) ≈ zr(loss) atol = 1.0e-7
            @test typeof(c) == typeof(loss)
            @test zr(fit(ZeroRateCurve(fill(0.03, 4), [1.0, 2.0, 5.0, 10.0], spline), cmt, method)) ≈ zr(c) atol = 1.0e-13
        end
        spline == Spline.Linear() && @test zr(fit(spline, cmt, Fit.LeastSquares())) ≈ zr(fit(spline, cmt, Fit.Bootstrap())) atol = 1.0e-14
        @test fit(spline, cmt, Fit.LeastSquares(); extrapolation = :flat_zero).extrapolation == :flat_zero
    end

    @testset "quotes that no curve reprices" begin
        # A quadratic spline through alternating zero-coupon and par yields at 15%: least squares returns a
        # local minimum, 9% from repricing.
        t = [0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30]
        y = 0.15 .+ 0.01 .* (1 .- exp.(-t ./ 5))
        qq = [isodd(i) ? ZCBYield(y[i], t[i]) : CMTYield(y[i], t[i]) for i in eachindex(t)]
        c = fit(Spline.Quadratic(), qq, Fit.LeastSquares())
        @test certified(objective(Spline.Quadratic(), qq, float.(t)), zr(c))
        @test maximum(abs(pv(c, q.instrument) - q.price) / q.price for q in qq) > 0.01
    end

    @testset "robustness: random indices" begin
        # (nine knots each, so that the Hessians of the checks compile once per interpolation)
        knotsets = (knots, [1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 19.0], [0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 9.0, 13.0, 20.0])
        converged = Dict{Any, Int}()
        unexplained = String[]
        for seed in 1:20
            kn = knotsets[mod1(seed, 3)]
            q = lsq_index(80; seed)
            # every fourth index also fits a smaller subset
            routes = seed % 4 == 0 ? (q, q[1:40]) : (q,)
            for spline in (Spline.Linear(), Spline.Quadratic(), Spline.Cubic(), Spline.MonotoneConvex(), Spline.PCHIP()), qx in routes
                c = try
                    fit(ZeroRateCurve(fill(0.04, length(kn)), kn, spline), qx, Fit.LeastSquares())
                catch e
                    e
                end
                if spline isa Union{Spline.MonotoneConvex, Spline.PCHIP}
                    # They switch formulas, so a minimum can lie on a switch: such a fit throws, and says so.
                    converged[spline] = get(converged, spline, 0) + !(c isa Exception)
                    (c isa Exception && !(c isa FitConvergenceError && occursin("switches formulas", sprint(showerror, c)))) &&
                        push!(unexplained, "index $seed, $(length(qx)) bonds, $spline: " * sprint(showerror, c))
                else
                    # Some of these have directions in which the residuals' curvature offsets most of JᵀJ: the
                    # first-order search alone ends up to 2e-7 from the minimum there.
                    @test !(c isa Exception) && certified(objective(spline, qx, kn), zr(c); atol = 1.0e-11)
                end
            end
        end
        @test isempty(unexplained)   # (shown when not)
        @test converged[Spline.MonotoneConvex()] >= 8 && converged[Spline.PCHIP()] >= 8   # of 25 each
    end

    @testset "a fit is accepted on the full Hessian, not on a Gauss-Newton step" begin
        # One knot and two zero-coupon prices, chosen so that at 5% the residuals' curvature cancels JᵀJ but
        # for 1e-5, with a gradient of 1e-9: a Gauss-Newton step there is 5.6e-11, while the minimum is
        # 1.1e-6 away, where the curvature is 1.8e-3. (From a review of an earlier least-squares fit.) A
        # one-knot fit starts at 5%.
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
        z1 = kn(ps)
        @test z1 ≈ zmin atol = 1.0e-12
        @test curv(z1) > 0 && abs(grad(z1) / curv(z1)) <= 1.0e-12
        # the derivatives of the minimum, dz/dp = -(∂grad/∂p)/curv, and central differences of the fit. A
        # price lower by 3e-9 (the first) or 3e-10 (the second) removes this minimum, so the difference step
        # is 3e-12: its truncation and the fitted rates' rounding (1e-13) each move a difference by less
        # than 1e-4 of it.
        dz = ForwardDiff.gradient(kn, ps)
        @test dz ≈ -ts .* exp.(-ts .* zmin) ./ curv(zmin) rtol = 1.0e-6
        ε = 3.0e-12
        @test dz ≈ [(kn(ps .+ ε .* (1:2 .== i)) - kn(ps .- ε .* (1:2 .== i))) / 2ε for i in 1:2] rtol = 1.0e-3
    end

    @testset "quotes that do not determine the parameters throw" begin
        # fewer quotes than knots: the Hessian is singular
        @test_throws "not positive definite" fit(c0, qs[1:5], Fit.LeastSquares())
        @test_throws FitConvergenceError fit(c0, qs[1:5], Fit.LeastSquares())
        cl = ZeroRateCurve(fill(0.05, 4), [1.0, 2.0, 3.0, 40.0], Spline.Linear())
        @test_throws "the knot rate at 40.0 does not move any quote's price" fit(cl, ZCBPrice.([0.97, 0.94, 0.9, 0.88], [1.0, 2.0, 2.5, 3.0]), Fit.LeastSquares())
        # nor does a zero weight on the only quote that moves a knot
        cz = ZeroRateCurve(fill(0.05, 3), [1.0, 2.0, 3.0], Spline.Linear())
        zqs = ZCBPrice.([0.97, 0.94, 0.941, 0.9, 0.91], [1.0, 2.0, 2.0, 3.0, 3.0])
        @test_throws "the knot rate at 3.0" fit(cz, zqs, Fit.LeastSquares(weights = [1.0, 1.0, 1.0, 0.0, 0.0]))
        @test_throws ArgumentError fit(cz, zqs, Fit.LeastSquares(weights = zeros(5)))
    end

    @testset "any contract" begin
        sub = qs[1:60]
        c5 = ZeroRateCurve(fill(0.04, 5), [1.0, 3.0, 5.0, 10.0, 20.0], Spline.Cubic())
        ref = zr(fit(c5, sub, Fit.LeastSquares()))
        # contracts that `fit` cannot look inside
        @test zr(fit(c5, [Quote(q.price, LsqOpaque(q.instrument)) for q in sub], Fit.LeastSquares())) == ref
        # a quote vector of abstract element type, mixing contract types
        mixed = Any[sub..., ZCBPrice(0.97, 0.5), ParYield(0.05, 12.0), Quote(1.9, Composite(Cashflow(1.0, 1.0), Cashflow(1.0, 2.5)))]
        @test certified(objective(Spline.Cubic(), mixed, [1.0, 3.0, 5.0, 10.0, 20.0]), zr(fit(c5, mixed, Fit.LeastSquares())))
        # option quotes: one volatility for five strikes that imply different ones
        vol(k) = 0.2 + 0.02 * (k - 1.0)
        calls = [Option.EuroCall(CommonEquity(), k, 1.0) for k in 0.8:0.1:1.2]
        oq = [Quote(pv(Equity.BlackScholesMerton(0.03, 0.01, Volatility.Constant(vol(c.strike))), c), c) for c in calls]
        σ = fit(Equity.BlackScholesMerton(0.03, 0.01, Volatility.Constant()), oq, Fit.LeastSquares()).σ.σ
        sse(s) = sum(abs2, pv(Equity.BlackScholesMerton(0.03, 0.01, Volatility.Constant(s)), q.instrument) - q.price for q in oq) / 2
        @test abs(ForwardDiff.derivative(sse, σ)) <= 1.0e-12 * ForwardDiff.derivative(s -> ForwardDiff.derivative(sse, s), σ)
        @test σ ≈ fit(Equity.BlackScholesMerton(0.03, 0.01, Volatility.Constant()), oq).σ.σ atol = 1.0e-6
    end

    @testset "any model" begin
        zq(m) = [Quote(pv(m, Cashflow(1.0, t)), Cashflow(1.0, t)) for t in 0.5:0.5:30]
        # one parameter, four prices
        q4 = ZCBPrice([0.9, 0.8, 0.7, 0.6])
        m = fit(Yield.Constant(), q4, Fit.LeastSquares())
        half(r) = sum(abs2, exp(-r * t) - q.price for (t, q) in zip(1:4, q4)) / 2
        r = m.rate.continuous_value
        @test abs(ForwardDiff.derivative(half, r)) <= 1.0e-12
        @test m isa Yield.Constant
        @test r ≈ fit(Yield.Constant(), q4).rate.continuous_value atol = 1.0e-8
        # Nelson–Siegel recovers the curve that priced the quotes, from the package's default start
        truth = Yield.NelsonSiegel(2.0, 0.058, -0.02, 0.03)
        ns = fit(Yield.NelsonSiegel(), zq(truth), Fit.LeastSquares())
        @test [ns.τ₁, ns.β₀, ns.β₁, ns.β₂] ≈ [2.0, 0.058, -0.02, 0.03] atol = 1.0e-10
        # chosen variables: the decay time stays where the model has it
        betas = (@optic(_.β₀) => -10.0 .. 10.0, @optic(_.β₁) => -10.0 .. 10.0, @optic(_.β₂) => -10.0 .. 10.0)
        nb = fit(Yield.NelsonSiegel(2.0, 0.0, 0.0, 0.0), zq(truth), Fit.LeastSquares(); variables = betas)
        @test nb.τ₁ == 2.0
        @test [nb.β₀, nb.β₁, nb.β₂] ≈ [0.058, -0.02, 0.03] atol = 1.0e-10
        # a short-rate model: at least as good a fit as the loss fit's
        vq = zq(ShortRate.Vasicek(0.3, 0.045, 0.012, Continuous(0.02)))
        v0 = ShortRate.Vasicek(0.1, 0.03, 0.01, Continuous(0.03))
        vsse(v) = sum(abs2, pv(v, q.instrument) - q.price for q in vq)
        @test vsse(fit(v0, vq, Fit.LeastSquares())) <= vsse(fit(v0, vq)) * (1 + 1.0e-8)
        # a minimum on or beyond a variable's bound throws
        tight = (@optic(_.rate.continuous_value) => 0.0 .. 0.05,)
        @test_throws "upper end of its interval" fit(Yield.Constant(), q4, Fit.LeastSquares(); variables = tight)
        # more variables than the quotes determine: these prices came from a Nelson–Siegel curve
        @test_throws FitConvergenceError fit(Yield.NelsonSiegelSvensson(), zq(truth), Fit.LeastSquares())
    end

    @testset "derivatives" begin
        sub = qs[1:60]
        c5 = ZeroRateCurve(fill(0.04, 5), [1.0, 3.0, 5.0, 10.0, 20.0], Spline.Cubic())
        bonds = [q.instrument for q in sub]
        h = 1.0e-4
        central(f, x, i) = (e = zeros(length(x)); e[i] = h; (f(x .+ e) .- f(x .- e)) ./ 2h)
        # quote prices
        kp(p) = zr(fit(c5, Quote.(p, bonds), Fit.LeastSquares()))
        prices = [q.price for q in sub]
        Jp = ForwardDiff.jacobian(kp, prices)
        @test Jp[:, 1] ≈ central(kp, prices, 1) rtol = 1.0e-6
        @test Jp[:, 37] ≈ central(kp, prices, 37) rtol = 1.0e-6
        # the values are the primal fit's exactly
        duals = ForwardDiff.Dual{:lsq}.(prices, 1.0)
        zd = zr(fit(c5, Quote.(duals, bonds), Fit.LeastSquares()))
        @test ForwardDiff.value.(zd) == kp(prices)
        # A dual number built by hand, with a tag ForwardDiff never ordered: the fit does not mix it with
        # its own.
        @test ForwardDiff.partials.(zd, 1) ≈ Jp * ones(length(prices)) rtol = 1.0e-10
        # coupons, inside the contracts
        kc(x) = zr(fit(c5, [Quote(q.price, Bond.Fixed(x[i], Periodic(2), maturity(q))) for (i, q) in enumerate(sub)], Fit.LeastSquares()))
        coupons = [b.coupon_rate for b in bonds]
        @test ForwardDiff.jacobian(kc, coupons)[:, 5] ≈ central(kc, coupons, 5) rtol = 1.0e-6
        # weights
        w = [1 / (q.price * maturity(q))^2 for q in sub]
        kw(x) = zr(fit(c5, sub, Fit.LeastSquares(weights = x)))
        Jw = ForwardDiff.jacobian(kw, w)
        @test Jw[:, 11] ≈ central(kw, w, 11) rtol = 1.0e-5
        # weights in a vector of abstract element type carry their derivatives too, and a zero weight has one
        @test ForwardDiff.jacobian(x -> kw(Real[x...]), w) ≈ Jw rtol = 1.0e-12
        w0 = [iseven(i) ? 0.0 : w[i] for i in eachindex(w)]
        @test ForwardDiff.jacobian(kw, w0)[:, 2] ≈ (kw(w0 .+ 1.0e-6 .* (eachindex(w0) .== 2)) .- kw(w0)) ./ 1.0e-6 rtol = 1.0e-3
        # quoted yields, through a valuation, with more quotes than knots
        tenors = [1, 2, 3, 5, 7, 10]
        rates = [0.03, 0.032, 0.034, 0.036, 0.037, 0.038]
        cfs = Cashflow.([3.0, 3.0, 103.0], [2.0, 5.0, 10.0])
        c3 = ZeroRateCurve(fill(0.03, 3), [1.0, 4.0, 10.0], Spline.Linear())
        value(x) = pv(fit(c3, CMTYield.(x, tenors), Fit.LeastSquares()), cfs)
        g = ForwardDiff.gradient(value, rates)
        @test g ≈ [(value(rates .+ h .* (1:6 .== i)) - value(rates .- h .* (1:6 .== i))) / 2h for i in 1:6] rtol = 1.0e-6
        # the extrapolation forward, which prices the quotes beyond the last knot
        ce(f) = ZeroRateCurve(fill(0.03, 3), [1.0, 4.0, 7.0], Spline.Linear(); extrapolation = Yield.FlatForwardAt(Continuous(f)))
        tail(f) = discount(fit(ce(f), CMTYield.(rates, tenors), Fit.LeastSquares()), 12.0)
        @test ForwardDiff.derivative(tail, 0.04) ≈ (tail(0.04 + h) - tail(0.04 - h)) / 2h rtol = 1.0e-6
        # a loss fit with another number of knots is still not differentiated
        @test_throws DimensionMismatch ForwardDiff.gradient(x -> pv(fit(c3, CMTYield.(x, tenors), Fit.Loss(abs2)), cfs), rates)
        # any model: Nelson–Siegel parameters with respect to prices
        truth = Yield.NelsonSiegel(2.0, 0.058, -0.02, 0.03)
        cashflows = [Cashflow(1.0, t) for t in 1.0:1.0:30.0]
        nsp(p) = (m = fit(Yield.NelsonSiegel(), Quote.(p, cashflows), Fit.LeastSquares()); [m.τ₁, m.β₀, m.β₁, m.β₂])
        p0 = [pv(truth, c) for c in cashflows] .* (1 .+ 0.001 .* sinpi.(eachindex(cashflows) ./ 3))
        Jn = ForwardDiff.jacobian(nsp, p0)
        hn = 1.0e-6
        @test Jn[:, 7] ≈ (nsp(p0 .+ hn .* (eachindex(p0) .== 7)) .- nsp(p0 .- hn .* (eachindex(p0) .== 7))) ./ 2hn rtol = 1.0e-5

        @testset "what is not differentiated throws" begin
            # a dual number where the fit cannot see it
            opaque(a) = [Quote(0.97, LsqOpaque(Cashflow(a, 1.0))), Quote(0.94, LsqOpaque(Cashflow(a, 2.0)))]
            c1 = ZeroRateCurve([0.03], [1.0], Spline.Linear())
            @test_throws "Fit.LeastSquares differentiates" ForwardDiff.derivative(a -> zero(fit(c1, opaque(a), Fit.LeastSquares()), 1.0), 1.0)
            # in a field of the model that is not a variable
            oq = [Quote(0.08, Option.EuroCall(CommonEquity(), 1.0, 1.0)), Quote(0.05, Option.EuroCall(CommonEquity(), 1.1, 1.0))]
            @test_throws "Fit.LeastSquares differentiates" ForwardDiff.derivative(
                r -> fit(Equity.BlackScholesMerton(r, 0.01, Volatility.Constant()), oq, Fit.LeastSquares()).σ.σ, 0.03
            )
            # second derivatives
            @test_throws "first-order" ForwardDiff.hessian(x -> value(x), rates)
            # maturities
            @test_throws "maturities or cashflow times" ForwardDiff.derivative(t -> zero(fit(c1, [ZCBPrice(0.97, t), ZCBPrice(0.94, 2.0)], Fit.LeastSquares()), 1.0), 1.0)
            # a fit on a kink of its interpolation: flat quotes make adjacent forwards equal
            flat(x) = discount(fit(Spline.MonotoneConvex(), ZCBYield.(x, [1.0, 2.0, 5.0, 10.0]), Fit.LeastSquares()), 7.0)
            @test flat(fill(0.04, 4)) ≈ exp(-7 * log(1.04))
            @test_throws "switches shape" ForwardDiff.gradient(flat, fill(0.04, 4))
        end
    end

    @testset "a fit inside a fit" begin
        # A loss fit whose quote is priced through a least-squares fit: the inner fit is differentiated with
        # respect to the outer solve's dual numbers, and each solve has tags of its own.
        refit = LsqRefit(qs[1:30])
        target = pv(Yield.Constant(Continuous(0.03)), refit)
        m = fit(Yield.Constant(), [Quote(target, refit)], Fit.Loss(abs2); solve_kwargs = (; g_tol = 1.0e-14))
        @test m.rate.continuous_value ≈ 0.03 atol = 1.0e-8
        # A least-squares fit certifies its minimum with second derivatives, which a fit inside it does not have.
        @test_throws "first-order" fit(Yield.Constant(), [Quote(target, refit)], Fit.LeastSquares())
    end

    @testset "FX forwards fit their foreign knot curve by least squares" begin
        eurusd = FX.Pair(:EUR, :USD)
        usd = Yield.Constant(Continuous(0.05))
        truth = FX.Forwards(eurusd, 1.1, usd, Yield.Constant(Continuous(0.03)))
        outrights = [FX.Outright(eurusd, forward(truth, t), t) for t in 0.5:0.5:6.0]
        foreign = ZeroRateCurve(fill(0.02, 3), [1.0, 3.0, 6.0], Spline.Linear())
        m = fit(FX.Forwards(eurusd, 1.1, usd, foreign), outrights, Fit.LeastSquares())
        @test knot_tenors(m.foreign) == [1.0, 3.0, 6.0]
        @test all(abs(pv(m, q.instrument)) < 1.0e-12 for q in outrights)   # a flat curve is representable
        @test zr(m.foreign) ≈ fill(0.03, 3) atol = 1.0e-12
    end
end
