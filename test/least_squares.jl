# Least-squares fits of knot curves to quotes. least_squares_generic.jl (any contract, any model) and
# least_squares_derivatives.jl continue these tests: they are separate files so that they run in parallel.
using ForwardDiff
using Accessors, IntervalSets

include("helpers/least_squares_setup.jl")   # `lsq_index`, `lsq_knots`, `zr`, `errors`, `certified`

# A contract type that `fit` knows nothing about.
struct LsqNoMaturity <: FinanceCore.AbstractContract end

struct LsqLinearPrice <: FinanceCore.AbstractContract end
FinanceCore.present_value(m::Yield.Constant, ::LsqLinearPrice) = m.rate.continuous_value

# The price is flat both at an interior minimum (r=0) and at a different value at r=1.
struct LsqSmoothstep <: FinanceCore.AbstractContract end
function FinanceCore.present_value(m::Yield.Constant, ::LsqSmoothstep)
    r = m.rate.continuous_value
    return r^2 * (3 - 2r)
end

# Two parameters seen through either of two sets of coordinates: `a` and `b` themselves, or their sum and
# difference. A contract whose price is flat in the first where it is zero (-2 exp(-s) + exp(-2s)), and a
# short bill that depends on the second.
struct LsqPair{A, B, C}
    a::A
    b::B
    rotated::C
end
struct LsqFlatPrice <: FinanceCore.AbstractContract end
struct LsqShortBill <: FinanceCore.AbstractContract end
lsq_s(m::LsqPair) = m.rotated ? m.a + m.b : m.a
lsq_d(m::LsqPair) = m.rotated ? m.a - m.b : m.b
FinanceCore.present_value(m::LsqPair, ::LsqFlatPrice) = -2 * exp(-lsq_s(m)) + exp(-2 * lsq_s(m))
FinanceCore.present_value(m::LsqPair, ::LsqShortBill) = exp(-0.01 * lsq_d(m))

@testset "least-squares fits" begin
    qs = lsq_index(200)
    knots = lsq_knots

    @testset "$spline: a certified minimum of the squared pricing residuals" for spline in (Spline.Linear(), Spline.Quadratic(), Spline.Cubic())
        c0 = ZeroRateCurve(fill(0.04, length(knots)), knots, spline)
        c = fit(c0, qs, Fit.LeastSquares())
        @test knot_tenors(c) == knots
        @test c.extrapolation == c0.extrapolation
        @test certified(errors(spline, qs, knots), zr(c), qs)
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
        @test certified(errors(Spline.Cubic(), qs, knots), zw, qs, w)
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
        @test_throws "finite and non-negative" fit(c0, qs, Fit.LeastSquares(weights = -w))
        @test_throws "finite and non-negative" fit(c0, qs, Fit.LeastSquares(weights = [NaN; w[2:end]]))
        @test_throws "finite and non-negative" fit(c0, qs, Fit.LeastSquares(weights = [Inf; w[2:end]]))
        # the loss fit's solver options are not least-squares options
        @test_throws MethodError fit(c0, qs, Fit.LeastSquares(); solve_kwargs = (; maxiters = 10))
    end

    @testset "a weight of zero leaves its quote out" begin
        # The quote is not priced at all: one that cannot be priced does not disturb the fit of the others.
        good, bad = ZCBPrice(0.9, 1.0), ZCBPrice(NaN, 2.0)
        alone = fit(Yield.Constant(), [good], Fit.LeastSquares())
        @test fit(Yield.Constant(), [good, bad], Fit.LeastSquares(weights = [1.0, 0.0])) == alone
        @test_throws "not finite" fit(Yield.Constant(), [good, bad], Fit.LeastSquares(weights = [1.0, 1.0]))
        # nor does it get a knot in a spline fit
        @test isequal(fit(Spline.Linear(), [good, bad], Fit.LeastSquares(weights = [1.0, 0.0])), fit(Spline.Linear(), [good], Fit.LeastSquares()))
        for model in (Yield.Constant(), Spline.Linear(), ZeroRateCurve([0.05], [1.0], Spline.Linear()))
            method = Fit.LeastSquares(weights = [1.0, 0.0])
            @test discount(fit(model, [good, Quote(1.0, LsqNoMaturity())], method), 1.0) ≈ 0.9
            excluded_time(x) = discount(fit(model, [ZCBPrice(x, 1.0), Quote(1.0, Cashflow(1.0, zero(x) + 2))], method), 1.0)
            @test ForwardDiff.derivative(excluded_time, 0.9) ≈ 1
        end
        # A zero weight has a derivative on a curve with knots of its own (see "derivatives"), but not in a
        # spline fit, where a positive weight would add a knot.
        two = ZCBPrice.([0.9, 0.8], [1.0, 2.0])
        @test_throws "would add a knot" ForwardDiff.derivative(x -> discount(fit(Spline.Linear(), two, Fit.LeastSquares(weights = [1.0, x])), 1.5), 0.0)
        @test ForwardDiff.derivative(x -> discount(fit(Spline.Linear(), two, Fit.LeastSquares(weights = [1.0, x])), 1.5), 1.0) ≈ 0 atol = 1.0e-10   # an exact fit
        # A zero weight that only looks differentiated (a dual number whose partials are zero) is left out all
        # the same, and does not stand in the way of another derivative.
        blank(p) = discount(fit(Spline.Linear(), ZCBPrice.([p, 0.8], [1.0, 2.0]), Fit.LeastSquares(weights = [1.0, 0 * p])), 1.0)
        @test ForwardDiff.derivative(blank, 0.9) ≈ 1 rtol = 1.0e-10
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
        @test certified(errors(Spline.Quadratic(), qq, float.(t)), zr(c), qq)
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
                    # first-order search alone ends up to 2e-7 from the minimum there, in knot rates.
                    @test !(c isa Exception) && certified(errors(spline, qx, kn), zr(c), qx)
                end
            end
        end
        @test isempty(unexplained)   # (shown when not)
        @test converged[Spline.MonotoneConvex()] >= 8 && converged[Spline.PCHIP()] >= 8   # of 25 each
    end

    @testset "one quote: a one-knot MonotoneConvex curve has no formula switch" begin
        # One knot gives a flat curve for every rate, so the certificate's Hessian exists. A two-year par
        # bond with a 5% semiannual coupon has the continuous rate 2 log(1 + 5%/2); its interim coupons are
        # discounted inside the curve's only interval.
        q = CMTYield(0.05, 2)
        c = fit(Spline.MonotoneConvex(), [q], Fit.LeastSquares())
        @test rate(only(knot_rates(c))) ≈ 2 * log1p(0.05 / 2) atol = 1.0e-13
        @test pv(c, q.instrument) ≈ q.price atol = 1.0e-13
    end

    @testset "two quotes at one yield: a flat two-knot MonotoneConvex curve has no formula switch" begin
        # Both intervals of a two-knot curve are end intervals, whose deviations keep one ratio as the
        # knot rates move: the interim coupons of two par bonds at 4% are priced where the two discrete
        # forwards are equal, and the certificate's Hessian exists there.
        qq = CMTYield.(0.04, [2, 4])
        c = fit(Spline.MonotoneConvex(), qq, Fit.LeastSquares())
        @test zr(c) ≈ fill(2 * log1p(0.04 / 2), 2) atol = 1.0e-13
        @test maximum(abs(pv(c, q.instrument) - q.price) for q in qq) < 1.0e-13
        # and the fit has derivatives there: the sensitivities of the fitted rates to the prices invert
        # the sensitivities of the prices to the rates
        refit(p) = zr(fit(Spline.MonotoneConvex(), [Quote(p[i], qq[i].instrument) for i in eachindex(qq)], Fit.LeastSquares()))
        prices(z) = [pv(ZeroRateCurve(z, [2.0, 4.0], Spline.MonotoneConvex()), q.instrument) for q in qq]
        @test ForwardDiff.jacobian(prices, zr(c)) * ForwardDiff.jacobian(refit, [1.0, 1.0]) ≈ I atol = 1.0e-12
    end

    @testset "a search that stops on a formula switch throws a FitConvergenceError" begin
        # Bonds priced off a flat curve have their minimum at that curve, where all of MonotoneConvex's
        # adjacent forwards are equal. A search started there does not move, and the certificate asks for
        # a Hessian that the interpolation does not supply. (A fit reaches such a point by rounding: about
        # one in a hundred of the fits of index 19 above, with its prices moved in their eleventh digit.)
        flat = [Quote(pv(Yield.Constant(Continuous(0.05)), q.instrument), q.instrument) for q in lsq_index(40)]
        spline = Spline.MonotoneConvex()
        rates = FinanceModels.__KnotRates(spline, Yield.__interpolation(spline, length(knots)), knots, :flat_forward)
        err = try
            FinanceModels.__least_squares(rates, fill(0.05, length(knots)), flat, flat, Fit.LeastSquares(), "fit")
        catch e
            e
        end
        @test err isa FitConvergenceError
        @test occursin("could not supply the second derivatives needed to certify this fit", sprint(showerror, err))
        @test occursin("switches formulas", sprint(showerror, err))
    end

    @testset "a fit is accepted on the full Hessian, not on a Gauss-Newton step" begin
        # At the 5% start, residual curvature cancels all but 1e-5 of JᵀJ. The Gauss-Newton
        # step is only 5.6e-11, but the minimum is 1.1e-6 away and has curvature 1.8e-3.
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
        # The derivatives are those of the optimality condition at the returned rate, dz/dp = -(∂grad/∂p)/curv.
        # Compared with the minimum's own, and with central differences of the fit. A price lower by 3e-9
        # (the first) or 3e-10 (the second) removes this minimum, so the difference step is 3e-12: its
        # truncation and the fitted rates' rounding each move a difference by less than 1e-4 of it.
        dz = ForwardDiff.gradient(kn, ps)
        @test dz ≈ -ts .* exp.(-ts .* zmin) ./ curv(zmin) rtol = 1.0e-6
        ε = 3.0e-12
        @test dz ≈ [(kn(ps .+ ε .* (1:2 .== i)) - kn(ps .- ε .* (1:2 .== i))) / 2ε for i in 1:2] rtol = 1.0e-3
    end

    @testset "Hessian and curvature uncertainty across a partial chunk" begin
        # Each residual has Hessian 2aᵢaᵢ′. Mixed signs distinguish the uncertainty
        # bound (entrywise absolute values before summing) from the signed Hessian.
        A = [1.0 -2.0 3.0; -2.0 1.0 1.0]
        θ, w, noise = [1.0, 2.0, 3.0], [3.0, 2.0], [0.25, 0.5]
        residual(x) = (A * x) .^ 2 .- [2, 1]
        outer = ForwardDiff.Tag(Val(:outer_lsq_test), Float64)
        inner = ForwardDiff.Tag(Val(:inner_lsq_test), Float64)
        H, R = FinanceModels.__hessian(residual, w, θ, outer, inner, noise, ForwardDiff.Chunk{2}())
        @test H == A' * Diagonal(w .* (4 .* (A * θ) .^ 2 .+ 2 .* residual(θ))) * A
        @test R == 2 .* abs.(A)' * Diagonal(noise) * abs.(A)
    end

    @testset "tolerances are on prices, quote by quote" begin
        # Ordinary price spacing at the starting model must not widen the final
        # relative tolerance or make an interior solution look like a bound.
        linear(p, start) = fit(
            Yield.Constant(Continuous(start)), [Quote(p, LsqLinearPrice())], Fit.LeastSquares();
            variables = (@optic(_.rate.continuous_value) => 0.0 .. 2.0e15,)
        ).rate.continuous_value
        @test linear(0.123, 1.0e10) ≈ 0.123 rtol = 1.0e-10
        @test linear(0.03, 1.0e14) ≈ 0.03 rtol = 1.0e-10
        # A cashflow of 1e12 today has the same value on every curve. Adding it must preserve
        # both the one-year discount factor and its derivative with respect to the bill price.
        c1 = ZeroRateCurve([0.02], [1.0], Spline.Linear())
        today = Quote(1.0e12, Cashflow(1.0e12, 0.0))
        df(p, others...) = discount(fit(c1, [ZCBPrice(p, 1.0), others...], Fit.LeastSquares()), 1.0)
        @test df(0.9, today) == df(0.9)
        @test df(0.9) ≈ 0.9 rtol = 1.0e-10
        @test ForwardDiff.derivative(p -> df(p, today), 0.9) ≈ 1 rtol = 1.0e-10
        # One that they move a little does have a say, and does not loosen the others' tolerance: the same
        # cashflow a billionth of a year out. Its price has steps of 1e-4, which leave the rate uncertain by
        # 1e-7, so no rate is a minimum to within 1e-10 of the one-year zero's price, and the fit says so
        # rather than stop where the larger quote alone is satisfied.
        soon = Quote(1.0e12, Cashflow(1.0e12, 1.0e-9))
        @test_throws "still moves a price" fit(c1, [ZCBPrice(0.9, 1.0), soon], Fit.LeastSquares())
        # At zero rates PV equals spot minus strike. Even a single representable spot increment
        # must be fitted when the quote calls for it.
        eurusd = FX.Pair(:EUR, :USD)
        flat = Yield.Constant(Continuous(0.0))
        spot = 1.0e12
        model = FX.Forwards(eurusd, spot, flat, flat)
        outright = FX.Forward(eurusd, spot, 1.0)
        for price in (0.5, eps(spot))
            fx = fit(model, [Quote(price, outright)], Fit.LeastSquares(); variables = (@optic(_.spot) => 1.0 .. 1.0e13,))
            @test pv(fx, outright) == price
        end
    end

    @testset "quotes that do not determine the parameters throw" begin
        # A parameter with a zero pricing derivative is determined when the residuals' curvature says so: the
        # present value of -2 at year 1 and +1 at year 2 is (exp(-r) - 1)² - 1, flat at r = 0, so for a price
        # of -1.1 the sum of squares has its only minimum there, with a Hessian of 0.2.
        flat = [Quote(-1.1, Composite(Cashflow(-2.0, 1.0), Cashflow(1.0, 2.0)))]
        @test abs(fit(Yield.Constant(Continuous(0.0)), flat, Fit.LeastSquares()).rate.continuous_value) <= 1.0e-10
        # From elsewhere the price test alone would stop short: at a rate of 1e-8 a correction of the whole
        # distance moves the price by 1e-16. A parameter that the curvature determines is refined until its
        # correction settles, because the fit's derivatives depend on where it is. Every price below -1 has
        # its minimum at a rate of zero, so the rate's derivative with respect to the price is zero; at a
        # rate of 1e-8 the implicit formula gives 9 for a price of -1.000000001.
        for price in (-1.1, -1.000000001)
            curved(p) = fit(Yield.Constant(Continuous(0.01)), [Quote(p, flat[1].instrument)], Fit.LeastSquares())
            @test curved(price).rate.continuous_value ≈ 0 atol = 1.0e-12
            @test pv(curved(price), flat[1].instrument) ≈ -1 atol = 1.0e-14
            @test ForwardDiff.derivative(p -> curved(p).rate.continuous_value, price) ≈ 0 atol = 1.0e-9
        end
        flat_rate(p, start) = fit(Yield.Constant(Continuous(start)), [Quote(p, flat[1].instrument)], Fit.LeastSquares()).rate.continuous_value
        for (p, start, tol) in ((-1 - 1.0e-13, 0.01, 0.01), (-1 - 1.0e-9, 5.0e-7, 1.0e-6))
            @test abs(flat_rate(p, start)) < 1.0e-12
            @test abs(ForwardDiff.derivative(x -> flat_rate(x, start), p)) < tol
        end
        @test_throws FitConvergenceError ForwardDiff.derivative(p -> flat_rate(p, 0.01), -1.0)
        # One spacing below the turning price, r = 0 has a strictly positive Hessian.
        # Rounding the quote can remove that curvature, so only the plain fit is resolved.
        pflat = prevfloat(-1.0)
        @test flat_rate(pflat, 0.0) == 0.0
        @test_throws "two-sided sensitivity cannot be resolved" ForwardDiff.derivative(p -> flat_rate(p, 0.0), pflat)
        above = -1 + 1.0e-9
        root = sqrt(1 + above)
        @test ForwardDiff.derivative(p -> flat_rate(p, 0.01), above) ≈ inv(2root * (1 - root)) rtol = 1.0e-6
        vars = (@optic(_.rate.continuous_value) => -2.0 .. 1.0,)
        stationary(p) = fit(Yield.Constant(Continuous(0.0)), [Quote(p, LsqSmoothstep())], Fit.LeastSquares(); variables = vars).rate.continuous_value
        @test stationary(-0.1) == 0
        @test ForwardDiff.derivative(stationary, -0.1) == 0
        # A badly mispriced time-zero cashflow adds a constant to the objective. It must
        # preserve both an off-market minimum and the selected root when two roots exist.
        function add_constant(p, extra)
            q = Quote(p, flat[1].instrument)
            quotes = extra ? [q, Quote(1.0e12, Cashflow(1.0, 0.0))] : [q]
            return fit(Yield.Constant(Continuous(0.05)), quotes, Fit.LeastSquares()).rate.continuous_value
        end
        for p in (-1.23, -0.8)
            @test add_constant(p, true) ≈ add_constant(p, false) atol = 1.0e-10
            @test ForwardDiff.derivative(x -> add_constant(x, true), p) ≈ ForwardDiff.derivative(x -> add_constant(x, false), p) atol = 1.0e-9
        end
        # The direction in which the prices are flat need not be one parameter's. With a short bill that
        # the second parameter moves, the minimum is at s = d = 0 and ds/dprice is zero, whether the fit's
        # variables are s and d or a and b with s = a + b and d = a - b. In a and b each variable alone
        # moves the bill, while the common shift leaves the bill unchanged and is determined by curvature.
        ab = (@optic(_.a) => -1.0 .. 1.0, @optic(_.b) => -1.0 .. 1.0)
        for (rotated, a, b) in ((false, 0.01, 0.0), (true, 0.005, 0.005)), price in (-1.1, -1.000000001)
            s(p) = lsq_s(fit(LsqPair(a, b, rotated), [Quote(p, LsqFlatPrice()), Quote(1.0, LsqShortBill())], Fit.LeastSquares(); variables = ab))
            @test s(price) ≈ 0 atol = 1.0e-12
            @test ForwardDiff.derivative(s, price) ≈ 0 atol = 1.0e-9
        end
        # Seeding another, well-determined quote in the same AD lane cannot regularize a fold.
        joint(p) = lsq_s(fit(LsqPair(0.01, 0.0, false), [Quote(p, LsqFlatPrice()), Quote(1 + 1.0e8 * (p + 1), LsqShortBill())], Fit.LeastSquares(); variables = ab))
        @test_throws FitConvergenceError ForwardDiff.derivative(joint, -1.0)
        # (A price above -1 is repriced, at a rate where the price has a slope: -0.99 at a discount factor of 0.9.)
        sloped(p) = fit(Yield.Constant(Continuous(0.01)), [Quote(p, flat[1].instrument)], Fit.LeastSquares()).rate.continuous_value
        @test sloped(-0.99) ≈ -log(0.9) atol = 1.0e-10
        @test ForwardDiff.derivative(sloped, -0.99) ≈ 1 / (2 * 0.1 * 0.9) rtol = 1.0e-8
        # fewer quotes than knots: the Hessian is singular
        @test_throws "not positive definite" fit(c0, qs[1:5], Fit.LeastSquares())
        @test_throws FitConvergenceError fit(c0, qs[1:5], Fit.LeastSquares())
        cl = ZeroRateCurve(fill(0.05, 4), [1.0, 2.0, 3.0, 40.0], Spline.Linear())
        @test_throws "the knot rate at 40.0 does not move any quote's price" fit(cl, ZCBPrice.([0.97, 0.94, 0.9, 0.88], [1.0, 2.0, 2.5, 3.0]), Fit.LeastSquares())
        # nor does a zero weight on the only quote that moves a knot
        cz = ZeroRateCurve(fill(0.05, 3), [1.0, 2.0, 3.0], Spline.Linear())
        zqs = ZCBPrice.([0.97, 0.94, 0.941, 0.9, 0.91], [1.0, 2.0, 2.0, 3.0, 3.0])
        @test_throws "the knot rate at 3.0" fit(cz, zqs, Fit.LeastSquares(weights = [1.0, 1.0, 1.0, 0.0, 0.0]))
        @test_throws "the knot rate at 1.0 does not move any quote" fit(cz, zqs, Fit.LeastSquares(weights = zeros(5)))
    end
end
