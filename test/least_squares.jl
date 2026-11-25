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
struct LsqNoMaturity <: FinanceCore.AbstractContract end
struct LsqOpaque{C} <: FinanceCore.AbstractContract
    inner::C
end
FinanceCore.maturity(c::LsqOpaque) = FinanceCore.maturity(c.inner)
FinanceCore.present_value(m, c::LsqOpaque) = FinanceCore.present_value(m, c.inner)

struct LsqLinearPrice <: FinanceCore.AbstractContract end
FinanceCore.present_value(m::Yield.Constant, ::LsqLinearPrice) = m.rate.continuous_value

# The price is flat both at an interior minimum (r=0) and at a different value at r=1.
struct LsqSmoothstep <: FinanceCore.AbstractContract end
function FinanceCore.present_value(m::Yield.Constant, ::LsqSmoothstep)
    r = m.rate.continuous_value
    return r^2 * (3 - 2r)
end

# A contract priced through a least-squares fit of its own: the 5-year discount factor of a two-knot curve
# fitted to `quotes` with their prices scaled by the model's 1-year discount factor.
struct LsqRefit{Q} <: FinanceCore.AbstractContract
    quotes::Q
end
function FinanceCore.present_value(m, c::LsqRefit)
    scaled = [Quote(q.price * discount(m, 1.0), q.instrument) for q in c.quotes]
    return discount(fit(ZeroRateCurve(fill(0.04, 2), [2.0, 10.0], Spline.Linear()), scaled, Fit.LeastSquares()), 5.0)
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
    knots = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
    zr(c) = rate.(knot_rates(c))
    # the pricing errors, written independently of the fit
    errors(spline, q = qs, kn = knots) = z -> [pv(ZeroRateCurve(z, kn, spline), x.instrument) - x.price for x in q]
    # The fit's certificate, recomputed with ForwardDiff's own Hessian of half the weighted sum of squares:
    # positive definite, and a Newton correction of any one knot rate moves no quote's price by more than
    # 1e-10 of the quote's size.
    function certified(e, z, q = qs, w = ones(length(q)))
        ϕ(x) = sum(w .* abs2.(e(x))) / 2
        H, J, f = ForwardDiff.hessian(ϕ, z), ForwardDiff.jacobian(e, z), e(z)
        δ = H \ ForwardDiff.gradient(ϕ, z)
        size(i) = max(abs(q[i].price), abs(f[i] + q[i].price))
        return isposdef(Symmetric(H)) && all(abs(J[i, k] * δ[k]) <= 1.0e-10 * size(i) for i in axes(J, 1), k in axes(J, 2))
    end

    @testset "$spline: a certified minimum of the squared pricing residuals" for spline in (Spline.Linear(), Spline.Quadratic(), Spline.Cubic())
        c0 = ZeroRateCurve(fill(0.04, length(knots)), knots, spline)
        c = fit(c0, qs, Fit.LeastSquares())
        @test knot_tenors(c) == knots
        @test c.extrapolation == c0.extrapolation
        @test certified(errors(spline), zr(c))
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
        @test certified(errors(Spline.Cubic()), zw, qs, w)
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

    @testset "any contract" begin
        sub = qs[1:60]
        c5 = ZeroRateCurve(fill(0.04, 5), [1.0, 3.0, 5.0, 10.0, 20.0], Spline.Cubic())
        ref = zr(fit(c5, sub, Fit.LeastSquares()))
        # contracts that `fit` cannot look inside
        @test zr(fit(c5, [Quote(q.price, LsqOpaque(q.instrument)) for q in sub], Fit.LeastSquares())) == ref
        # a quote vector of abstract element type, mixing contract types
        mixed = Any[sub..., ZCBPrice(0.97, 0.5), ParYield(0.05, 12.0), Quote(1.9, Composite(Cashflow(1.0, 1.0), Cashflow(1.0, 2.5)))]
        @test certified(errors(Spline.Cubic(), mixed, [1.0, 3.0, 5.0, 10.0, 20.0]), zr(fit(c5, mixed, Fit.LeastSquares())), mixed)
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
        # a minimum beyond a variable's bound throws
        tight = (@optic(_.rate.continuous_value) => 0.0 .. 0.05,)
        @test_throws "upper end of its interval" fit(Yield.Constant(), q4, Fit.LeastSquares(); variables = tight)
        # So does one on the bound, which rounding would otherwise put just inside, with a two-sided
        # derivative: a one-year zero at a price of 1 is fitted by a rate of exactly 0.
        bounded(p) = fit(Yield.Constant(Continuous(0.02)), [ZCBPrice(p, 1.0)], Fit.LeastSquares(); variables = tight).rate.continuous_value
        @test_throws "lower end of its interval" bounded(1.0)
        @test_throws "lower end of its interval" ForwardDiff.derivative(bounded, 1.0)
        near = exp(-1.0e-5)
        @test bounded(near) ≈ 1.0e-5 atol = 1.0e-14
        @test ForwardDiff.derivative(bounded, near) ≈ -1 / near rtol = 1.0e-10
        # This forward requires a spot of exactly 100.001, the upper bound. A quote just
        # inside that bound has a two-sided spot sensitivity of one.
        eurusd = FX.Pair(:EUR, :USD)
        norates = Yield.Constant(Continuous(0.0))
        capped(p) = fit(FX.Forwards(eurusd, 100.0, norates, norates), [Quote(p, FX.Forward(eurusd, 100.0, 1.0))], Fit.LeastSquares(); variables = (@optic(_.spot) => 99.0 .. 100.001,)).spot
        @test_throws "upper end of its interval" capped(100.001 - 100)
        @test_throws "upper end of its interval" ForwardDiff.derivative(capped, 100.001 - 100)
        @test capped(0.00099) ≈ 100.00099 atol = 1.0e-10
        @test ForwardDiff.derivative(capped, 0.00099) ≈ 1 rtol = 1.0e-10
        # Widening the rate interval must preserve the fitted price of the 30-year zero.
        # A full Newton step from a zero rate would overshoot to 54.6.
        wide = (@optic(_.rate.continuous_value) => -1.0e9 .. 1.0e9,)
        @test fit(Yield.Constant(Continuous(0.0)), [ZCBPrice(1.8, 30.0)], Fit.LeastSquares(); variables = wide).rate.continuous_value ≈ -log(1.8) / 30 atol = 1.0e-14
        # Other number types: the solve works in Float64 at least, and keeps a wider type.
        function optimum(prices)
            r = big"0.1"
            for _ in 1:60
                d = [exp(-r * t) for t in 1:3]
                r -= sum(-t * d[t] * (d[t] - prices[t]) for t in 1:3) / sum(t^2 * d[t] * (2 * d[t] - prices[t]) for t in 1:3)
            end
            return r
        end
        for T in (Float32, BigFloat)
            prices = T.([0.91, 0.82, 0.75])
            vars = (@optic(_.rate.continuous_value) => T(-1) .. T(1),)
            fitted = fit(Yield.Constant(Continuous(T(0.05))), ZCBPrice.(prices, T.(1:3)), Fit.LeastSquares(); variables = vars).rate.continuous_value
            @test fitted isa (T == BigFloat ? BigFloat : Float64)
            @test abs(fitted - optimum(prices)) <= 1.0e-12
        end
        for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            ts = [1.0, 2.0, 4.0, 6.0]
            curve = ZeroRateCurve(fill(big"0.04", 4), ts, spline)
            prices = [exp(-t * (big"0.03" + big"0.001" * t)) for t in ts]
            first_rate(p) = first(zr(fit(curve, ZCBPrice.([p; prices[2:end]], ts), Fit.LeastSquares())))
            @test first_rate(prices[1]) isa BigFloat
            sensitivity = ForwardDiff.derivative(first_rate, prices[1])
            @test sensitivity isa BigFloat
            @test sensitivity ≈ -1 / prices[1] rtol = 1.0e-10
        end
        # More parameters than one ForwardDiff chunk, including the last partial chunk.
        ts17 = collect(1.0:17.0)
        prices17 = exp.(-0.03 .* ts17)
        long_curve = ZeroRateCurve(fill(0.04, 17), ts17, Spline.Linear())
        last_rate(p) = last(zr(fit(long_curve, ZCBPrice.([prices17[1:(end - 1)]; p], ts17), Fit.LeastSquares())))
        @test ForwardDiff.derivative(last_rate, last(prices17)) ≈ -1 / (17 * last(prices17)) rtol = 1.0e-10
        # A discrete set excludes the continuous optimum of 0.03 between its first two values.
        three = [ZCBPrice(exp(-0.03), 1.0)]
        @test_throws "variable 1 is not given an interval" fit(Yield.Constant(Continuous(0.05)), three, Fit.LeastSquares(); variables = (@optic(_.rate.continuous_value) => [0.0, 0.05, 0.1],))
        @test fit(Yield.Constant(Continuous(0.05)), three, Fit.LeastSquares(); variables = (@optic(_.rate.continuous_value) => 0.0 .. 0.1,)).rate.continuous_value ≈ 0.03 atol = 1.0e-12
        # The certificate measures a correction by the prices it moves, not in the variable's units: the
        # same contract in years and with cashflows 1e11 years out, where the rate that fits is 7e-12.
        # The two cashflows are worth d² - 3d for d = discount(t), and -1.23 at d = 0.49005.
        for t in (1.0, 1.0e11)
            far = [Quote(-1.23, Composite(Cashflow(-3.0, t), Cashflow(1.0, 2t)))]
            d(p) = discount(fit(Yield.Constant(Continuous(0.0)), [Quote(p, far[1].instrument)], Fit.LeastSquares()), t)
            @test d(-1.23) ≈ 0.4900495061637922 rtol = 1.0e-10
            @test ForwardDiff.derivative(d, -1.23) ≈ 1 / (2 * 0.4900495061637922 - 3) rtol = 1.0e-8
        end
        # What the method needs of a model: `variables`, and prices that ForwardDiff differentiates twice.
        @test_throws MethodError fit(Yield.SmithWilson(ufr = 0.03, α = 0.1), ZCBPrice([0.97, 0.94]), Fit.LeastSquares())
        hw = ShortRate.HullWhite(0.1, 0.015, Yield.Constant(Continuous(0.03)))
        swaptions = [Option.Swaption(e, e + 5.0, 0.03, 2) for e in 1.0:3.0]
        @test_throws "first-order" fit(ShortRate.HullWhite(0.05, 0.01, Yield.Constant(Continuous(0.03))), [Quote(pv(hw, s), s) for s in swaptions], Fit.LeastSquares())
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
