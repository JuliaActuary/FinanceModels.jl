# Derivatives of least-squares fits, a fit inside a fit, and FX forwards. Continues least_squares.jl.
using ForwardDiff

include("helpers/least_squares_setup.jl")   # `lsq_index`, `lsq_knots`, `LsqOpaque`, `zr`

# A contract priced through a least-squares fit of its own: the 5-year discount factor of a two-knot curve
# fitted to `quotes` with their prices scaled by the model's 1-year discount factor.
struct LsqRefit{Q} <: FinanceCore.AbstractContract
    quotes::Q
end
function FinanceCore.present_value(m, c::LsqRefit)
    scaled = [Quote(q.price * discount(m, 1.0), q.instrument) for q in c.quotes]
    return discount(fit(ZeroRateCurve(fill(0.04, 2), [2.0, 10.0], Spline.Linear()), scaled, Fit.LeastSquares()), 5.0)
end

@testset "least-squares fits: derivatives" begin
    qs = lsq_index(200)
    knots = lsq_knots

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

        # Outrights and basis swaps together: the residuals minimized are the quotes' own, whatever curve the
        # model holds. With a spot of 100, a one-year outright at 90 and a par basis leg under zero rates
        # are worth 100 DF - 90 and DF, so the weighted minimum is in closed form.
        zero = Yield.Constant(Continuous(0.0))
        mixed = [FX.Outright(eurusd, 90.0, 1.0), FX.ParBasisSwap(eurusd, 0.0, 1.0; reference = zero)]
        knot = ZeroRateCurve([0.02], [1.0], Spline.Linear())
        df(foreign, spot, method) = discount(fit(FX.Forwards(eurusd, spot, zero, foreign), mixed, method).foreign, 1.0)
        @test df(knot, 100.0, Fit.LeastSquares()) ≈ 9001 / 10001 atol = 1.0e-13
        @test df(Yield.Constant(Continuous(0.02)), 100.0, Fit.LeastSquares()) ≈ 9001 / 10001 atol = 1.0e-13
        closed(s) = (2 * 90 * s + 0.5) / (2 * s^2 + 0.5)   # weights of 2 and 0.5
        weighted = Fit.LeastSquares(weights = [2.0, 0.5])
        @test df(knot, 100.0, weighted) ≈ closed(100.0) atol = 1.0e-13
        @test ForwardDiff.derivative(s -> df(knot, s, weighted), 100.0) ≈ ForwardDiff.derivative(closed, 100.0) rtol = 1.0e-10
        @test_throws DimensionMismatch df(knot, 100.0, Fit.LeastSquares(weights = [1.0]))
        # Each outright must imply a positive discount factor, even when the average
        # of the implied factors would give a feasible least-squares optimum.
        invalid = [Quote(p, FX.Forward(eurusd, 1.0, 1.0)) for p in (-2.0, 4.0)]
        @test_throws "implies a non-positive" fit(FX.Forwards(eurusd, 1.0, zero, knot), invalid, Fit.LeastSquares())
        # A quote left out by its weight is not converted to foreign-curve cashflows either, with knots of the
        # curve's own or at the quotes' maturities.
        good, bad = FX.Outright(eurusd, 90.0, 1.0), FX.Outright(eurusd, NaN, 2.0)
        for foreign in (knot, Spline.Linear())
            fitted(quotes, method) = discount(fit(FX.Forwards(eurusd, 100.0, zero, foreign), quotes, method).foreign, 1.0)
            @test fitted([good, bad], Fit.LeastSquares(weights = [1.0, 0.0])) == fitted([good], Fit.LeastSquares())
            @test fitted([good], Fit.LeastSquares()) ≈ 0.9 atol = 1.0e-12
            # nor does it stop a derivative with respect to spot in the outright's foreign receipt
            other = [good, FX.Outright(eurusd, 85.0, 2.0)]
            @test ForwardDiff.derivative(s -> discount(fit(FX.Forwards(eurusd, s, zero, foreign), other, Fit.LeastSquares(weights = [1.0, 0.0])).foreign, 1.0), 100.0) ≈ -0.009 rtol = 1.0e-10
        end

        # A contract with no implied quote on the foreign curve is priced under the FX model itself: two
        # one-year forwards at 90 in one contract are worth 2 (spot DF - 90), so a price p is fitted by
        # DF = (180 + p) / (2 spot).
        pair_of = Composite(good.instrument, good.instrument)
        swapfit(foreign, spot, p; domestic = zero) = discount(fit(FX.Forwards(eurusd, spot, domestic, foreign), [Quote(p, pair_of)], Fit.LeastSquares()).foreign, 1.0)
        for foreign in (knot, Spline.Linear(), Yield.Constant(Continuous(0.02)))
            @test swapfit(foreign, 100.0, 0.0) ≈ 0.9 atol = 1.0e-12
        end
        for foreign in (knot, Spline.Linear())
            @test ForwardDiff.derivative(p -> swapfit(foreign, 100.0, p), 0.0) ≈ 1 / 200 rtol = 1.0e-10
            @test ForwardDiff.derivative(s -> swapfit(foreign, s, 0.0), 100.0) ≈ -0.009 rtol = 1.0e-10
            # the strikes inside the contract
            struck(k) = discount(fit(FX.Forwards(eurusd, 100.0, zero, foreign), [Quote(0.0, Composite(FX.Forward(eurusd, k, 1.0), good.instrument))], Fit.LeastSquares()).foreign, 1.0)
            @test ForwardDiff.derivative(struck, 90.0) ≈ 1 / 200 rtol = 1.0e-10
            # The discounted strike carries the domestic curve derivative: DF = 0.9 exp(-r) for outrights.
            @test ForwardDiff.derivative(r -> discount(fit(FX.Forwards(eurusd, 100.0, Yield.Constant(Continuous(r)), foreign), [good], Fit.LeastSquares()).foreign, 1.0), 0.0) ≈ -0.9 rtol = 1.0e-10
            @test_throws "Fit.LeastSquares differentiates" ForwardDiff.derivative(r -> swapfit(foreign, 100.0, 0.0; domestic = Yield.Constant(Continuous(r))), 0.0)
        end

        # Net PV can be many orders of magnitude smaller than either discounted leg.
        for spot in (1.0, 100.0, 1.0e6, 1.0e12), ratio in (0.5e-12, 0.8e-12, 1.0e-9, 1.0e-6)
            strike = spot * exp(-0.05)
            c = FX.Forward(eurusd, strike, 1.0)
            p = spot * ratio
            fx = fit(FX.Forwards(eurusd, spot, zero, knot), [Quote(p, c)], Fit.LeastSquares())
            @test pv(fx, c) ≈ p atol = 4eps(strike) rtol = 1.0e-10
            if spot == 1.0e12 && p == 0.5
                @test pv(fx, c) == p
            end
        end
        for r in (0.049, 0.05, 0.051)
            truth = FX.Forwards(eurusd, 1.1, usd, Yield.Constant(Continuous(r)))
            ts = 0.5:0.5:6.0
            qs = [FX.Outright(eurusd, forward(truth, t), t) for t in ts]
            curve = ZeroRateCurve(fill(0.05, 4), [1.0, 2.0, 4.0, 6.0], Spline.Linear())
            fitted = fit(FX.Forwards(eurusd, 1.1, usd, curve), qs, Fit.LeastSquares())
            @test maximum(t -> abs(discount(fitted.foreign, t) - exp(-r * t)), ts) < 1.0e-12
        end
        function kink(r)
            ts = [1.0, 2.0, 3.0, 5.0]
            curve = ZeroRateCurve(fill(0.03, 4), ts, Spline.MonotoneConvex())
            qs = map(ts) do t
                f = FX.Forward(eurusd, 0.0, t)
                Quote(2exp(-r * t), Composite(f, f))
            end
            return discount(fit(FX.Forwards(eurusd, 1.0, zero, curve), qs, Fit.LeastSquares()).foreign, 1.5)
        end
        @test kink(0.03) ≈ exp(-0.03 * 1.5)
        @test_throws ArgumentError ForwardDiff.derivative(kink, 0.03)
        @test_throws "switches shape" ForwardDiff.derivative(kink, 0.03)

        priced = FX.Forwards(eurusd, 100.0, zero, Yield.Constant(Continuous(0.03)))
        for q in mixed
            legs = FinanceModels.__foreign_lsq_quote(priced, q)
            @test FinanceModels.__quote_value(priced, q) - q.price ≈ FinanceModels.__quote_value(priced.foreign, legs) - legs.price rtol = 1.0e-12
        end
    end
end
