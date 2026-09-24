using ForwardDiff
const DI = FinanceModels.DifferentiationInterface

struct ImplicitFitTestTag end
struct OtherImplicitFitTestTag end

# A contract type the calibration walker does not know, carrying its cashflow inside.
struct WrappedContract{C} <: FinanceCore.AbstractContract
    inner::C
end
FinanceCore.maturity(c::WrappedContract) = FinanceCore.maturity(c.inner)
FinanceCore.present_value(m, c::WrappedContract, t = 0.0) = FinanceCore.present_value(m, c.inner, t)

@testset "differentiable fits" begin
    tenors = [1.0, 2.0, 3.0, 5.0, 10.0]
    rates = [0.03, 0.032, 0.035, 0.037, 0.04]
    # cash flows inside, at, and beyond the knot grid
    cfs = Cashflow.([4.0, 4.0, 4.0, 4.0, 104.0, 50.0], [0.5, 2.0, 3.5, 5.0, 10.0, 20.0])
    fitters = (
        "Linear bootstrap" => qs -> fit(Spline.Linear(), qs, Fit.Bootstrap()),
        "BSpline(1) bootstrap" => qs -> fit(Spline.BSpline(1), qs, Fit.Bootstrap()),
        "Linear loss" => qs -> fit(Spline.Linear(), qs),
        "Cubic loss" => qs -> fit(Spline.Cubic(), qs),
        "PCHIP loss" => qs -> fit(Spline.PCHIP(), qs),
        "MonotoneConvex loss" => qs -> fit(Spline.MonotoneConvex(), qs),
    )
    families = (
        "CMTYield" => CMTYield,
        "OISYield" => OISYield,
        "annual ParYield" => (r, t) -> ParYield(r, t; frequency = 1),
    )
    # derivative accuracy tracks each fit's repricing precision
    rtol(name) = occursin("bootstrap", name) ? 1.0e-7 : occursin("MonotoneConvex", name) ? 1.0e-4 : 1.0e-6
    # Reference refits for finite differences. Loss fits stop at optimizer tolerance
    # (residuals near 1e-11), which is noise of order 1e-4 in a difference quotient, so each
    # refit is polished to an exact solve of the repricing conditions with Newton's method.
    function exact(curve, qs)
        z = collect(knot_rates(curve))
        R(z) = [pv(reconstruct(curve; rates = z), q.instrument) - q.price for q in qs]
        for _ in 1:10
            r = R(z)
            maximum(abs, r) < 1.0e-15 && break
            z = z - ForwardDiff.jacobian(R, z) \ r
        end
        return reconstruct(curve; rates = z)
    end

    @testset "primal values are the primal fit, bitwise: $name" for (name, fitter) in fitters
        qs = CMTYield.(rates, tenors)
        dual_qs = CMTYield.(ForwardDiff.Dual{ImplicitFitTestTag}.(rates, 1.0), tenors)
        c, cd = fitter(qs), fitter(dual_qs)
        @test typeof(c) != typeof(cd)
        @test ForwardDiff.value.(knot_rates(cd)) == knot_rates(c)
        @test knot_tenors(cd) == knot_tenors(c) && cd.spline == c.spline
        # primal quotes take the primal path and return a primal curve
        @test eltype(knot_rates(c)) === Float64
    end

    @testset "zero-coupon quotes: ∂zᵢ/∂pⱼ = -δᵢⱼ/(tᵢpᵢ): $name" for (name, fitter) in fitters
        prices = exp.(-rates .* tenors)
        J = ForwardDiff.jacobian(p -> collect(knot_rates(fitter(ZCBPrice.(p, tenors)))), prices)
        @test J ≈ [i == j ? -1 / (tenors[i] * prices[i]) : 0.0 for i in eachindex(tenors), j in eachindex(tenors)] rtol = rtol(name)
    end

    @testset "coupon quotes against refit central differences: $fname, $name" for (fname, family) in families,
            (name, fitter) in fitters
        value(x) = pv(fitter(family.(x, tenors)), cfs)
        reference(x) = (qs = family.(x, tenors); pv(exact(fitter(qs), qs), cfs))
        g = ForwardDiff.gradient(value, rates)
        h = 1.0e-6
        fd = map(eachindex(rates)) do i
            e = h .* (eachindex(rates) .== i)
            (reference(rates .+ e) - reference(rates .- e)) / 2h
        end
        @test g ≈ fd rtol = 1.0e-6
    end

    @testset "method independence: linear loss equals linear bootstrap" begin
        value(fitter) = x -> pv(fitter(CMTYield.(x, tenors)), cfs)
        @test ForwardDiff.gradient(value(qs -> fit(Spline.Linear(), qs)), rates) ≈
            ForwardDiff.gradient(value(qs -> fit(Spline.Linear(), qs, Fit.Bootstrap())), rates) rtol = 1.0e-6
    end

    @testset "tag, chunk, and backend invariance" begin
        value(x) = pv(fit(Spline.Linear(), CMTYield.(x, tenors), Fit.Bootstrap()), cfs)
        g = ForwardDiff.gradient(value, rates)
        for chunk in (1, 2, 5)
            cfg = ForwardDiff.GradientConfig(value, rates, ForwardDiff.Chunk{chunk}())
            @test ForwardDiff.gradient(value, rates, cfg) ≈ g rtol = 1.0e-12
        end
        @test DI.gradient(value, DI.AutoForwardDiff(), rates) ≈ g rtol = 1.0e-12
        # quote prices and rates in the same differentiation
        mixed(x) = pv(fit(Spline.Linear(), [ZCBPrice(x[1], 1.0); CMTYield.(x[2:end], tenors[2:end])], Fit.Bootstrap()), cfs)
        x0 = [exp(-rates[1]); rates[2:end]]
        h = 1.0e-6
        fd = [(mixed(x0 .+ h .* (1:5 .== i)) - mixed(x0 .- h .* (1:5 .== i))) / 2h for i in 1:5]
        @test ForwardDiff.gradient(mixed, x0) ≈ fd rtol = 1.0e-6
    end

    @testset "a dual extrapolation forward" begin
        value(f) = pv(
            fit(
                Spline.Linear(), CMTYield.(rates, tenors), Fit.Bootstrap();
                extrapolation = Yield.FlatForwardAt(Continuous(f))
            ), cfs
        )
        h = 1.0e-6
        @test ForwardDiff.derivative(value, 0.04) ≈ (value(0.04 + h) - value(0.04 - h)) / 2h rtol = 1.0e-6
        # the knot rates do not depend on a tail beyond every quote's cash flows
        c = fit(
            Spline.Linear(), CMTYield.(rates, tenors), Fit.Bootstrap();
            extrapolation = Yield.FlatForwardAt(Continuous(ForwardDiff.Dual{ImplicitFitTestTag}(0.04, 1.0)))
        )
        @test all(z -> ForwardDiff.partials(z, 1) == 0, knot_rates(c))
        @test ForwardDiff.partials(rate(zero(c, 30.0)), 1) > 0
    end

    @testset "FX forwards: implied foreign quotes carry the derivatives" begin
        eurusd = FX.Pair(:EUR, :USD)
        usd = Yield.Constant(Continuous(0.05))
        ts = [0.5, 1.0, 2.0, 3.0]
        points = [0.0112, 0.022, 0.043, 0.062]
        value(spot) = begin
            qs = FX.Outright.(Ref(eurusd), spot .+ points, ts)
            m = fit(FX.Forwards(eurusd, spot, usd, Spline.Linear()), qs, Fit.Bootstrap())
            discount(m.foreign, 2.5)
        end
        h = 1.0e-6
        @test ForwardDiff.derivative(value, 1.1) ≈ (value(1.1 + h) - value(1.1 - h)) / 2h rtol = 1.0e-6
    end

    @testset "loud errors" begin
        value(x) = pv(fit(Spline.Linear(), CMTYield.(x, tenors), Fit.Bootstrap()), cfs)
        # nested differentiation
        @test_throws "nested dual numbers" ForwardDiff.derivative(
            s -> ForwardDiff.derivative(u -> value(rates .+ s .+ u), 0.0), 0.0
        )
        # maturities and cash flow times
        @test_throws "maturities" ForwardDiff.derivative(
            t -> discount(fit(Spline.Linear(), [ZCBPrice(0.95, 2.0 + t)], Fit.Bootstrap()), 1.0), 0.0
        )
        @test_throws "maturities" ForwardDiff.derivative(
            t -> discount(fit(Spline.Linear(), [Quote(0.95, Cashflow(1.0, 2.0 + t))]), 1.0), 0.0
        )
        # dual numbers from two differentiations
        mixed = [
            ZCBPrice(ForwardDiff.Dual{ImplicitFitTestTag}(0.97, 1.0), 1.0),
            ZCBPrice(ForwardDiff.Dual{OtherImplicitFitTestTag}(0.93, 1.0), 2.0),
        ]
        @test_throws "different ForwardDiff calls" fit(Spline.Linear(), mixed, Fit.Bootstrap())
        # a dual number inside a contract type the calibration does not know
        wrapped = [Quote(0.97, WrappedContract(Cashflow(ForwardDiff.Dual{ImplicitFitTestTag}(1.0, 1.0), 1.0)))]
        @test_throws "WrappedContract" fit(Spline.Linear(), wrapped, Fit.Bootstrap())
        @test_throws "WrappedContract" fit(Spline.Linear(), wrapped)
        # a loss fit that does not reprice its quotes has no implicit derivative
        offset = Fit.Loss(x -> (x - 0.01)^2)
        @test_throws "does not reprice" ForwardDiff.gradient(x -> pv(fit(Spline.Linear(), CMTYield.(x, tenors), offset), cfs), rates)
        @test fit(Spline.Linear(), CMTYield.(rates, tenors), offset) isa Yield.Spline    # the primal fit is fine
        # quote prices that do not determine a knot
        free = [ZCBPrice(ForwardDiff.Dual{ImplicitFitTestTag}(0.97, 1.0), 1.0), Quote(0.0, Cashflow(0.0, 2.0))]
        @test_throws "singular" fit(Spline.Linear(), free, Fit.Bootstrap())
        # model calibrations that are not differentiated
        @test_throws "reconstruct" ForwardDiff.derivative(
            x -> discount(fit(Yield.NelsonSiegel(), ZCBPrice.(exp.(-rates .* tenors) .+ x, tenors)), 3.0), 0.0
        )
        @test_throws "reconstruct" ForwardDiff.derivative(
            x -> discount(fit(ZeroRateCurve(rates, tenors), ZCBPrice.(exp.(-rates .* tenors) .+ x, tenors)), 3.0), 0.0
        )
    end
end
