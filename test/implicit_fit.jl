using ForwardDiff
const DI = FinanceModels.DifferentiationInterface

struct ImplicitFitTestTag end
struct OtherImplicitFitTestTag end

# A contract type the calibration walker does not know, carrying its cashflow inside.
struct WrappedContract{C} <: FinanceCore.AbstractContract
    inner::C
end
FinanceCore.maturity(c::WrappedContract) = FinanceCore.maturity(c.inner)
FinanceCore.present_value(m, c::WrappedContract) = FinanceCore.present_value(m, c.inner)

# Prices that depend on a dual number only in part of the parameter space (an external review's
# counterexample and its spline analogue). An optimizer fit's parameters take the dual number's type
# from the start (SciMLBase promotes them); a bootstrap meets it mid-solve, in its residual.
struct ConditionalPrice{T} <: FinanceCore.AbstractContract
    adjustment::T
end
FinanceCore.present_value(m::Yield.Constant, c::ConditionalPrice) = (r = rate(m.rate); r < 0.1 ? r : r + c.adjustment)
struct ConditionalZCB{T} <: FinanceCore.AbstractContract
    adjustment::T
    time::Float64
end
FinanceCore.maturity(c::ConditionalZCB) = c.time
FinanceCore.present_value(m::Yield.AbstractYieldModel, c::ConditionalZCB) =
    (d = discount(m, c.time); d < 0.96 ? d : d * (1 + c.adjustment))

# Contracts valued by an inner fit on quotes that depend on the outer candidate model.
struct NestedConstant <: FinanceCore.AbstractContract end
FinanceCore.present_value(m::Yield.Constant, ::NestedConstant) =
    discount(fit(Yield.Constant(0.0), [ZCBPrice(discount(m, 1.0), 1.0)]), 2.0)
struct NestedSpline <: FinanceCore.AbstractContract end
FinanceCore.present_value(m::Yield.Constant, ::NestedSpline) =
    discount(fit(Spline.Linear(), [ZCBPrice(discount(m, 1.0), 1.0), ZCBPrice(discount(m, 3.0), 3.0)]), 2.0)
# The same inner fit on quote vectors whose element type hides the dual number.
struct NestedAny <: FinanceCore.AbstractContract end
FinanceCore.present_value(m::Yield.Constant, ::NestedAny) =
    discount(fit(Yield.Constant(0.0), Any[ZCBPrice(discount(m, 1.0), 1.0)]), 2.0)
struct NestedQuotes <: FinanceCore.AbstractContract end
FinanceCore.present_value(m::Yield.Constant, ::NestedQuotes) =
    discount(fit(Yield.Constant(0.0), Quote[ZCBPrice(discount(m, 1.0), 1.0)]), 2.0)
# An inner fit on a `V[...]` quote vector whose price depends on the outer candidate only once the
# inner rate reaches 10%, after the inner fit's start.
struct NestedConditional{V, O} <: FinanceCore.AbstractContract
    optimizer::O
end
NestedConditional{V}(optimizer) where {V} = NestedConditional{V, typeof(optimizer)}(optimizer)
FinanceCore.present_value(m::Yield.Constant, c::NestedConditional{V}) where {V} = rate(
    fit(
        Yield.Constant(0.0), V[Quote(0.3, ConditionalPrice(rate(m.rate)))];
        optimizer = c.optimizer, solve_kwargs = (; g_tol = 1.0e-20, maxiters = 10_000)
    ).rate
)

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
        z = rate.(knot_rates(curve))
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
        @test ForwardDiff.value.(rate.(knot_rates(cd))) == rate.(knot_rates(c))
        @test knot_tenors(cd) == knot_tenors(c) && cd.spline == c.spline
        # primal quotes take the primal path and return a primal curve
        @test eltype(knot_rates(c)) === Rate{Float64, Continuous}
    end

    @testset "zero-coupon quotes: ∂zᵢ/∂pⱼ = -δᵢⱼ/(tᵢpᵢ): $name" for (name, fitter) in fitters
        prices = exp.(-rates .* tenors)
        J = ForwardDiff.jacobian(p -> rate.(knot_rates(fitter(ZCBPrice.(p, tenors)))), prices)
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

    @testset "notionals do not matter" begin
        # Mixed notionals (every other quote scaled by n) give the same knot rates and the same
        # derivatives with respect to the unit prices; conditioning and the repricing check are
        # in knot-rate units.
        p = exp.(-rates .* tenors)
        value(qs) = pv(fit(Spline.Linear(), qs, Fit.Bootstrap()), cfs)
        g = ForwardDiff.gradient(x -> value(ZCBPrice.(x, tenors)), p)
        for n in (1.0e-10, 1.0e8)
            scaled(x) = value([Quote(n^(i % 2) * x[i], Cashflow(n^(i % 2), tenors[i])) for i in eachindex(x)])
            @test ForwardDiff.gradient(scaled, p) ≈ g rtol = 1.0e-10
        end
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
        @test all(z -> ForwardDiff.partials(z, 1) == 0, rate.(knot_rates(c)))
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
        # a dual number inside a contract type the calibration does not know (made by ForwardDiff,
        # whose tags order against the optimizer's own)
        wrapped(a) = [Quote(0.97, WrappedContract(Cashflow(a, 1.0)))]
        @test_throws "reached the solve" ForwardDiff.derivative(a -> discount(fit(Spline.Linear(), wrapped(a), Fit.Bootstrap()), 1.0), 1.0)
        @test_throws "reached the solve" ForwardDiff.derivative(a -> discount(fit(Spline.Linear(), wrapped(a)), 1.0), 1.0)
        # a loss fit that does not reprice its quotes has no implicit derivative
        offset = Fit.Loss(x -> (x - 0.01)^2)
        @test_throws "does not reprice" ForwardDiff.gradient(x -> pv(fit(Spline.Linear(), CMTYield.(x, tenors), offset), cfs), rates)
        @test fit(Spline.Linear(), CMTYield.(rates, tenors), offset) isa Yield.Spline    # the primal fit is fine
        # The check bounds one Newton correction of the knots, not the scaled residual: this
        # ill-conditioned pair reprices to about 2e-9 with its second knot 1e-5 off.
        ε = 1.0e-4
        pair(p1, p2) = [ZCBPrice(p1, 1.0), Quote(p2, FinanceCore.Composite(Cashflow(1.0, 1.0), Cashflow(ε, 2.0)))]
        zx = [0.03, 0.04]
        px = [exp(-zx[1]), exp(-zx[1]) + ε * exp(-2 * zx[2])]
        dual_pair = pair(ForwardDiff.Dual{ImplicitFitTestTag}(px[1], 1.0, 0.0), ForwardDiff.Dual{ImplicitFitTestTag}(px[2], 0.0, 1.0))
        off_curve = ZeroRateCurve(zx .+ [0.0, 1.0e-5], [1.0, 2.0], Spline.Linear())
        @test maximum(q -> abs(pv(off_curve, q.instrument) - q.price), pair(px...)) < 1.0e-8
        @test_throws "Newton correction" FinanceModels.__implicit_knot_curve(off_curve, dual_pair, pair(px...), off_curve.extrapolation)
        # and it says how to tighten the fit
        @test_throws "solve_kwargs" FinanceModels.__implicit_knot_curve(off_curve, dual_pair, pair(px...), off_curve.extrapolation)
        exact_curve = ZeroRateCurve(zx, [1.0, 2.0], Spline.Linear())
        @test Yield.knot_rates(FinanceModels.__implicit_knot_curve(exact_curve, dual_pair, pair(px...), exact_curve.extrapolation)) isa AbstractVector{<:Rate{<:ForwardDiff.Dual, Continuous}}
        # quote prices that do not determine a knot
        free = [ZCBPrice(ForwardDiff.Dual{ImplicitFitTestTag}(0.97, 1.0), 1.0), Quote(0.0, Cashflow(0.0, 2.0))]
        @test_throws "singular" fit(Spline.Linear(), free, Fit.Bootstrap())
        # model calibrations that are not differentiated
        @test_throws "reconstruct" ForwardDiff.derivative(
            x -> discount(fit(Yield.NelsonSiegel(), ZCBPrice.(exp.(-rates .* tenors) .+ x, tenors)), 3.0), 0.0
        )
    end

    @testset "dual-number provenance: only the solve's own dual numbers" begin
        # A caller's dual number that reaches a primal solve through data it does not
        # differentiate throws where it appears, rather than returning the derivative of the
        # solver's iterations (0.0 under NelderMead).
        provenance = "reached the solve through data it does not differentiate"
        NM = FinanceModels.OptimizationOptimJL.NelderMead
        # its price depends on the dual number only once the rate reaches 10%
        conditional(a) = [Quote(0.3, ConditionalPrice(a))]
        @test_throws provenance ForwardDiff.derivative(a -> rate(fit(Yield.Constant(0.0), conditional(a); optimizer = NM()).rate), 0.1)
        @test_throws provenance ForwardDiff.derivative(a -> rate(fit(Yield.Constant(0.0), conditional(a)).rate), 0.1)
        # a dual strike: FinanceModels' own option contract, not stripped by the calibration
        bsm = Equity.BlackScholesMerton(0.01, 0.02, Volatility.Constant())
        price = pv(Equity.BlackScholesMerton(0.01, 0.02, 0.15), Option.EuroCall(CommonEquity(), 1.0, 1.0))
        vol(K; kw...) = fit(bsm, [Quote(price, Option.EuroCall(CommonEquity(), K, 1.0))]; kw...).σ.σ
        @test_throws provenance ForwardDiff.derivative(K -> vol(K; optimizer = NM()), 1.0)
        @test_throws provenance ForwardDiff.derivative(vol, 1.0)
        # a dual amount inside a wrapper the calibration does not strip
        ns(a) = [Quote(p, Forward(0.0, Cashflow(a, t))) for (p, t) in zip([0.97, 0.93, 0.86, 0.75], [1.0, 2.0, 4.0, 8.0])]
        @test_throws provenance ForwardDiff.derivative(a -> discount(fit(Yield.NelsonSiegel(), ns(a)), 3.0), 1.0)
        # spline fits: the loss fit through its parameters, the bootstrap mid-solve
        knot(a) = [Quote(0.9704, ConditionalZCB(a, 1.0)), ZCBPrice(0.94, 2.0)]
        @test_throws provenance ForwardDiff.derivative(a -> discount(fit(Spline.Linear(), knot(a)), 1.5), 0.01)
        @test_throws provenance ForwardDiff.derivative(a -> discount(fit(Spline.Linear(), knot(a), Fit.Bootstrap()), 1.5), 0.01)
        # nested fits with the same loss type: the inner fit cannot take the outer fit's dual
        # number for its own (overlapping fits lease different owners, part of their tags' types)
        @test_throws provenance fit(Yield.Constant(0.0), [Quote(0.94, NestedConstant())])
        # also when the quote vectors' element type hides the dual number from the loss's and the
        # quotes' types (this fitted a 0% rate)
        @test_throws provenance fit(Yield.Constant(0.0), Any[Quote(0.94, NestedAny())])
        @test_throws provenance fit(Yield.Constant(0.0), Quote[Quote(0.94, NestedQuotes())])
        # and when, in addition, the inner price depends on the outer fit only away from the inner start
        for V in (Any, Quote), optimizer in (NM(), FinanceModels.OptimizationOptimJL.LBFGS())
            @test_throws provenance fit(Yield.Constant(0.0), V[Quote(0.2, NestedConditional{V}(optimizer))])
        end
        # a dual number outside the quotes (a fixed model field, the loss function) leaves the
        # parameters primal and reaches only the loss value
        qbsm = [Quote(price, Option.EuroCall(CommonEquity(), 1.0, 1.0))]
        @test_throws provenance ForwardDiff.derivative(r -> fit(Equity.BlackScholesMerton(r, 0.02, Volatility.Constant()), qbsm).σ.σ, 0.01)
        @test_throws provenance ForwardDiff.derivative(
            a -> discount(fit(Spline.Linear(), ZCBPrice.([0.97, 0.94, 0.9], [1.0, 2.0, 3.0]), Fit.Loss(x -> a * x^2)), 1.5), 1.0
        )
        # every owner is returned, also by the fits that threw
        @test !any(FinanceModels.__FIT_OWNERS)
        # and in a primal residual whose untyped captured data holds a dual number of its own tag's type
        hidden = Ref{Any}(0.0)
        residual = FinanceModels.__PrimalResidual(x -> x - hidden[], "residual", "")
        hidden[] = ForwardDiff.Dual{ForwardDiff.Tag{typeof(residual), Float64}}(0.5, 1.0)
        @test_throws provenance residual(0.25)
        # a legitimate nesting: an inner spline fit, differentiated implicitly, inside a generic fit
        outer = fit(Yield.Constant(0.0), [Quote(0.94, NestedSpline())])
        inner(m) = fit(Spline.Linear(), [ZCBPrice(discount(m, 1.0), 1.0), ZCBPrice(discount(m, 3.0), 3.0)])
        @test discount(inner(outer), 2.0) ≈ 0.94 atol = 1.0e-8
        # the solves keep their own dual numbers: second-order optimizers and implicit roots
        @test fit(bsm, [Quote(price, Option.EuroCall(CommonEquity(), 1.05, 1.0))]; optimizer = FinanceModels.OptimizationOptimJL.IPNewton()).σ.σ ≈
            fit(bsm, [Quote(price, Option.EuroCall(CommonEquity(), 1.05, 1.0))]).σ.σ rtol = 1.0e-6
        h = 1.0e-6
        iq(r) = implied_quote(Yield.Constant(Continuous(r)), ZCBYield, 5.0)
        @test ForwardDiff.derivative(iq, 0.03) ≈ (iq(0.03 + h) - iq(0.03 - h)) / 2h rtol = 1.0e-6
    end

    @testset "fit owners: distinct while solves overlap, reused afterwards" begin
        ready, release = Channel{Int}(4), Channel{Nothing}(4)
        tasks = map(1:4) do _
            Threads.@spawn begin
                owner = FinanceModels.__lease_fit_owner()
                try
                    put!(ready, owner)
                    take!(release)
                finally
                    FinanceModels.__return_fit_owner(owner)
                end
            end
        end
        @test allunique([take!(ready) for _ in 1:4])
        foreach(_ -> put!(release, nothing), 1:4)
        foreach(wait, tasks)
        @test !any(FinanceModels.__FIT_OWNERS)
        @test FinanceModels.__lease_fit_owner() == 1
        FinanceModels.__return_fit_owner(1)
    end

    @testset "refitting a knot curve" begin
        # one knot per quote: the refit is the spline fit, derivatives included
        refit(spline) = x -> pv(fit(ZeroRateCurve(fill(0.02, length(tenors)), tenors, spline), CMTYield.(x, tenors)), cfs)
        direct(spline) = x -> pv(fit(spline, CMTYield.(x, tenors)), cfs)
        for spline in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            @test ForwardDiff.gradient(refit(spline), rates) == ForwardDiff.gradient(direct(spline), rates)
        end
        # Another number of knots leaves the repricing conditions non-square: they do not
        # determine the knot rates' derivatives, and the implicit solve cannot be formed.
        for knots in ([1.0, 3.0, 10.0], [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0])
            c0 = ZeroRateCurve(fill(0.02, length(knots)), knots, Spline.Linear())
            @test fit(c0, CMTYield.(rates, tenors)) isa Yield.Spline    # the primal fit is fine
            @test_throws DimensionMismatch ForwardDiff.gradient(x -> pv(fit(c0, CMTYield.(x, tenors)), cfs), rates)
        end
    end
end
