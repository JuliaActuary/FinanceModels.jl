# Test-only solver that deterministically exercises `fit`'s unsuccessful-retcode path.
struct FailingFitOptimizer end
function FinanceModels.Optimization.solve(
        prob::FinanceModels.Optimization.OptimizationProblem,
        ::FailingFitOptimizer
    )
    return (
        u = prob.u0,
        retcode = FinanceModels.Optimization.SciMLBase.ReturnCode.Failure,
    )
end

# Regression tests from the 2026-06 ecosystem audit
# A contract defined through the documented `asfoldable` extension point
struct AsfoldableContract <: FinanceCore.AbstractContract end
Transducers.asfoldable(p::Projection{AsfoldableContract}) = [Cashflow(1.0, 1.0), Cashflow(2.0, 2.0)]

@testset "audit regressions" begin
    @testset "fit(spline, quotes, Fit.Loss) uses the supplied loss" begin
        qs = ZCBPrice([0.9, 0.8, 0.7])
        calls = Ref(0)
        counted = Fit.Loss(x -> (calls[] += 1; x^2))
        fit(Spline.Linear(), qs, counted)
        # the user-supplied loss was previously silently replaced with the default
        @test calls[] > 0
    end

    @testset "par with a non-representable stub maturity" begin
        c = Yield.Constant(0.04)
        @test rate(par(c, 4)) ≈ 0.03960780543711406
        @test par(c, 0.2; frequency = 4).compounding == Periodic(5)
        # 1/0.3 is not an integer → informative error instead of an InexactError
        @test_throws ArgumentError par(c, 0.3)
    end

    @testset "present_value past maturity is 0, not an error" begin
        b = Bond.Fixed(0.05, Periodic(1), 3)
        m = Yield.Constant(0.03)
        @test present_value(m, b, 4.0) == 0.0
        @test present_value(m, Projection(b, m, CashflowProjection()), 5.0) == 0.0
        # (n.b. `present_value(m, ::Cashflow, t)` hits FinanceCore's Cashflow method,
        # which discounts to time 0 regardless of the third argument — so the
        # empty-fold path is exercised via a single-cashflow Projection instead)
        @test present_value(m, Projection(Cashflow(10.0, 1.0), m, CashflowProjection()), 2.0) == 0.0
    end

    @testset "Forward contract shifts cashflow times" begin
        b = Bond.Fixed(0.05, Periodic(1), 2)
        fwd = Forward(1.0, b)
        cfs = collect(Projection(fwd, NullModel(), CashflowProjection()))
        @test [cf.time for cf in cfs] ≈ [2.0, 3.0]
        @test [cf.amount for cf in cfs] ≈ [0.05, 1.05]
        # for a flat curve, shifting all cashflows by Δ scales the PV by discount(Δ)
        m = Yield.Constant(0.03)
        @test present_value(m, fwd) ≈ present_value(m, b) * discount(m, 1.0)
        # maturity is the forward time plus the instrument's, so pv_mc's default horizon covers it
        @test maturity(fwd) == 3.0
        @test maturity(Forward(0.5, fwd)) == 3.5
        v = ShortRate.Vasicek(0.1, 0.05, 0.0, Continuous(0.03))
        @test pv_mc(v, fwd; n_scenarios = 1) ≈ present_value(v, fwd) rtol = 1.0e-6
    end

    @testset "Forward contract reads its models on its own clock" begin
        curve = ZeroRateCurve([0.02, 0.025, 0.03], [1.0, 2.0, 5.0], Spline.Linear())
        # A one-year floater starting at 2 fixes on the index rate from 2 to 3; it read 0 to 1 and paid
        # 1.0202 instead of 1.0305 here
        floater = Bond.Floating(0.0, Periodic(1), 1.0, :index)
        cfs = collect(Projection(Forward(2.0, floater), Dict(:index => curve)))
        @test [cf.time for cf in cfs] == [3.0]
        @test only(cfs).amount ≈ discount(curve, 2.0) / discount(curve, 3.0) rtol = 1.0e-14
        # Projecting a future-starting contract is projecting it on the models seen from its start
        # and shifting the cashflows
        coupon_floater = Bond.Floating(0.004, Periodic(2), 3.0, :index)
        for s in (0.0, 0.5, 2.0)
            actual = collect(Projection(Forward(s, coupon_floater), Dict(:index => curve)))
            relative = collect(Projection(coupon_floater, Dict(:index => Yield.ForwardStarting(curve, s))))
            @test [cf.amount for cf in actual] ≈ [cf.amount for cf in relative] rtol = 1.0e-14
            @test [cf.time for cf in actual] == [cf.time + s for cf in relative]
        end
        # An FX model is rebased with its forward rate at the start as spot, so a converted instrument
        # converts at its shifted payment times
        pair = FX.Pair(:EUR, :USD)
        fx = FX.Forwards(pair, 1.08, curve, Yield.Constant(Continuous(0.01)))
        @test forward(FinanceModels.__rebase(fx, 2.0), 1.0) ≈ forward(fx, 3.0) rtol = 1.0e-14
        store = Dict(:index => curve, :fx => fx)
        @test only(collect(Projection(Forward(2.0, FX.Converted(Cashflow(1.0, 1.0), pair, :fx)), store))).amount ≈
            forward(fx, 3.0) rtol = 1.0e-14

        # A model is rebased only when the projection reads it: a discount rate standing in for the
        # model, a NamedTuple store and an unused entry all work as before rebasing existed
        fixed = Bond.Fixed(0.05, Periodic(1), 2.0)
        for r in (0.03, Continuous(0.03))
            @test present_value(r, Forward(2.0, fixed)) ≈ present_value(r, fixed) * discount(r, 2.0) rtol = 1.0e-14
        end
        flows(c, models) = [(cf.amount, cf.time) for cf in collect(Projection(c, models))]
        @test flows(Forward(2.0, floater), (index = curve,)) == flows(Forward(2.0, floater), Dict(:index => curve))
        @test flows(Forward(2.0, floater), Dict(:index => curve, :equity => Equity.BlackScholesMerton(0.05, 0.02, 0.2))) ==
            flows(Forward(2.0, floater), Dict(:index => curve))

        # Time translation commutes with projection, and wrappers nest in any order: a forward at 0 is
        # the contract, nested starts add, forwarding a composite forwards its parts, and conversion
        # commutes with forwarding (each of these but the first threw a MethodError)
        same(x, y) = length(x) == length(y) && all(isapprox(a[1], b[1]; rtol = 1.0e-14) && a[2] ≈ b[2] for (a, b) in zip(x, y))
        @test same(flows(Forward(0.0, coupon_floater), store), flows(coupon_floater, store))
        @test same(flows(Forward(1.0, Forward(1.5, coupon_floater)), store), flows(Forward(2.5, coupon_floater), store))
        @test same(
            flows(Forward(2.0, FinanceCore.Composite(fixed, coupon_floater)), store),
            flows(FinanceCore.Composite(Forward(2.0, fixed), Forward(2.0, coupon_floater)), store)
        )
        @test same(
            flows(FX.Converted(Forward(2.0, coupon_floater), pair, :fx), store),
            flows(Forward(2.0, FX.Converted(coupon_floater, pair, :fx)), store)
        )
        @test same(
            flows(Forward(1.0, FinanceCore.Composite(fixed, [coupon_floater, fixed])), store),
            flows(FinanceCore.Composite(Forward(1.0, fixed), FinanceCore.Composite(Forward(1.0, coupon_floater), Forward(1.0, fixed))), store)
        )
        # A contract defined by `asfoldable` folds inside any wrapper, as alone (it threw inside a composite)
        custom = [(1.0, 1.0), (2.0, 2.0)]
        @test flows(AsfoldableContract(), NullModel()) == custom
        @test flows(FinanceCore.Composite(AsfoldableContract(), Cashflow(3.0, 3.0)), NullModel()) == [custom; (3.0, 3.0)]
        @test flows(Forward(1.0, AsfoldableContract()), NullModel()) == [(1.0, 2.0), (2.0, 3.0)]
        @test flows([AsfoldableContract()], NullModel()) == custom
        # A flat rate read as the index is the same from any start (a `Rate` store threw under `Forward`)
        flat = Dict(:index => Continuous(0.03))
        @test flows(Forward(2.0, floater), flat) == [(cf.amount, cf.time + 2.0) for cf in collect(Projection(floater, flat))]
        # Wrappers fold early-terminating and stateful transducers like a materialized vector
        nested = FinanceCore.Composite(Forward(1.0, fixed), FinanceCore.Composite(fixed, Forward(2.0, fixed)))
        for c in (nested, Forward(2.0, nested), [nested, nested])
            p = Projection(c)
            all_flows = collect(p)
            @test collect(p |> Transducers.Take(1)) == all_flows[1:1]
            @test collect(p |> Transducers.Map(cf -> cf.amount) |> Transducers.Scan(+)) ==
                collect(all_flows |> Transducers.Map(cf -> cf.amount) |> Transducers.Scan(+))
        end
        # A transducer over a contract applies in the order written, alone and inside every wrapper
        # (a projection applied a chain in reverse: [2.1, 4.1] here)
        bond = Bond.Fixed(0.05, Periodic(1), 2.0)
        dbl = Transducers.Map(cf -> Cashflow(2cf.amount, cf.time))
        inc = Transducers.Map(cf -> Cashflow(cf.amount + 1, cf.time))
        amounts(c, models = NullModel()) = [cf.amount for cf in collect(Projection(c, models))]
        @test amounts(bond |> dbl |> inc) == [1.1, 3.1] == [cf.amount for cf in collect(bond |> dbl |> inc)]
        @test amounts(FinanceCore.Composite(bond |> dbl |> inc, Cashflow(1.0, 3.0))) == [1.1, 3.1, 1.0]
        @test amounts([bond |> dbl |> inc, bond]) == [1.1, 3.1, 0.05, 1.05]
        # stateful and early-terminating transducers over a contract (they threw inside a projection),
        # and an outer stop after an inner one
        @test amounts(bond |> Transducers.Take(1)) == [0.05]
        @test amounts(FinanceCore.Composite(bond |> Transducers.Take(1), bond)) == [0.05, 0.05, 1.05]
        @test collect(Projection(bond |> Transducers.Map(cf -> cf.amount) |> Transducers.Scan(+))) ≈ [0.05, 1.1]
        two = Projection(FinanceCore.Composite(bond |> Transducers.Take(1), bond)) |> Transducers.Take(2)
        @test [cf.amount for cf in collect(two)] == [0.05, 0.05]
        @test [cf.amount for cf in collect(Projection(FinanceCore.Composite(bond, bond |> Transducers.Take(1))) |> Transducers.Take(1))] == [0.05]
        # a contract defined by `asfoldable` behind a transducer
        @test amounts(AsfoldableContract() |> Transducers.Map(cf -> Cashflow(-cf.amount, cf.time))) == [-1.0, -2.0]
        @test amounts(FinanceCore.Composite(AsfoldableContract() |> Transducers.Take(1), AsfoldableContract())) == [1.0, 1.0, 2.0]
        @test collect(AsfoldableContract() |> Transducers.Map(cf -> cf.amount)) == [1.0, 2.0]
        # and its value
        @test present_value(curve, Projection(bond |> dbl |> inc, curve, CashflowProjection())) ≈
            1.1 * discount(curve, 1.0) + 3.1 * discount(curve, 2.0) rtol = 1.0e-14
        # Its value is the floater's at the start, discounted to now
        @test present_value(curve, Projection(Forward(2.0, floater), Dict(:index => curve))) ≈ discount(curve, 2.0) rtol = 1.0e-14
    end

    @testset "ParSwapYield" begin
        # Swap fixed-leg conventions differ by market, so the frequency is required.
        @test_throws UndefKeywordError ParSwapYield(0.04, 5)
        q = ParSwapYield(0.04, 5; frequency = 4)
        @test q.price ≈ 1.0
        @test q.instrument.frequency == Periodic(4)
        @test ParSwapYield(0.04, 5; frequency = Periodic(1)).instrument.frequency == Periodic(1)
        # a Periodic Rate input must agree with an explicit frequency
        q2 = ParSwapYield(Periodic(0.04, 2), 5; frequency = 2)
        @test q2.instrument.frequency == Periodic(2)
        @test_throws ArgumentError ParSwapYield(Periodic(0.04, 2), 5; frequency = 4)
        # round-trip: a curve fit to par-swap quotes reprices them to root-finder
        # precision (the fitting-time curve and the returned curve are identical
        # for local interpolants, including the pinned t=0 knot)
        swap_rates = [0.02, 0.025, 0.03]
        for frequency in (1, 4)
            qs = ParSwapYield.(swap_rates, [1, 2, 3]; frequency)
            c = fit(Spline.Linear(), qs, Fit.Bootstrap())
            for (r, t) in zip(swap_rates, [1, 2, 3])
                @test rate(par(c, t; frequency)) ≈ r atol = 1.0e-10
            end
        end
        @test_throws UndefKeywordError ParSwapYield(swap_rates)
        @test ParSwapYield(swap_rates; frequency = 1) == ParSwapYield.(swap_rates, [1.0, 2.0, 3.0]; frequency = 1)
    end

    @testset "ParYield frequency" begin
        @test ParYield(0.04, 5).instrument.frequency == Periodic(2)
        @test ParYield(0.04, 5; frequency = 1).instrument.frequency == Periodic(1)
        # A Periodic rate sets its own frequency; a conflicting frequency is an error,
        # not a silently ignored keyword.
        @test ParYield(Periodic(0.04, 4), 5).instrument.frequency == Periodic(4)
        @test ParYield(Periodic(0.04, 4), 5; frequency = 4) == ParYield(Periodic(0.04, 4), 5)
        @test_throws ArgumentError ParYield(Periodic(0.04, 4), 5; frequency = 2)
        # Other rates convert to the requested frequency.
        @test ParYield(Continuous(0.04), 5; frequency = 1).instrument.coupon_rate ≈ exp(0.04) - 1
    end

    @testset "OISYield conventions" begin
        # One year or less: a single payment. Longer: annual payments on both legs.
        @test OISYield(0.04, 0.5).instrument == Bond.Fixed(0.0, Periodic(1), 0.5)
        @test OISYield(0.04, 0.5).price ≈ 1.04^-0.5
        # annual coupons at the quoted rate (a FinanceCore older than 2.6 converts the rate to
        # itself through continuous compounding, off by an ulp)
        b = OISYield(0.04, 5).instrument
        @test b isa Bond.Fixed && b.frequency == Periodic(1) && b.maturity == 5
        @test b.coupon_rate ≈ 0.04 rtol = 1.0e-15
        @test OISYield(0.04, 5).price == 1.0
        # A one-year single payment equals a one-coupon annual par swap.
        @test pv(Yield.Constant(0.04), OISYield(0.04, 1).instrument) ≈ OISYield(0.04, 1).price
        @test pv(Yield.Constant(Periodic(0.04, 1)), OISYield(0.04, 5).instrument) ≈ 1.0 atol = 1.0e-14
    end

    @testset "PCHIP and Akima ZeroRateCurve" begin
        zrates = [0.02, 0.025, 0.03, 0.035]
        tenors = [1.0, 2.0, 5.0, 10.0]
        @testset "$spline" for spline in (Spline.PCHIP(), Spline.Akima())
            zrc = ZeroRateCurve(zrates, tenors, spline)
            for (r, t) in zip(zrates, tenors)
                @test rate(zero(zrc, t)) ≈ r atol = 1.0e-10
            end
            @test 0 < discount(zrc, 3.0) < 1
        end
    end

    @testset "PCHIP and Akima loss fits do not stall at a flat seed" begin
        tenors = [1 / 12, 2 / 12, 3 / 12, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
        rates = [
            0.0375, 0.0372, 0.0372, 0.0372, 0.0364,
            0.0368, 0.0369, 0.038, 0.04, 0.0423,
        ]

        for quote_type in (CMTYield, ParYield), spline in (Spline.PCHIP(), Spline.Akima())
            qs = quote_type.(rates, tenors)
            curve = fit(spline, qs)
            max_price_error = maximum(
                abs,
                present_value(curve, q.instrument) - q.price for q in qs
            )
            @test max_price_error < 1.0e-6
        end
    end

    @testset "PCHIP and Akima loss fits on evenly spaced maturities" begin
        # A seed linear in the knot number is collinear on an evenly spaced grid, where Akima
        # switches formula and the loss derivatives are NaN; this fit used to fail to converge.
        t = collect(1.0:8.0)
        qs = CMTYield.(0.02 .+ 0.001 .* t .^ 2, t)
        for spline in (Spline.PCHIP(), Spline.Akima())
            curve = fit(spline, qs)
            @test maximum(abs(present_value(curve, q.instrument) - q.price) for q in qs) < 1.0e-9
        end
        # The seed is away from every formula switch, and the loss has finite first and second
        # derivatives there (coupon quotes also value the curve between its knots).
        grids = (
            collect(1.0:8.0), collect(1.0:30.0), [0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0],
            [0.0, 1.0, 2.0, 5.0], [0.1, 0.35, 1.2, 4.0, 4.5, 17.0],
        )
        for spline in (Spline.PCHIP(), Spline.Akima()), tenors in grids
            seed = FinanceModels.__knot_fit_seed(spline, tenors)
            @test !Yield.__near_kink(ZeroRateCurve(seed, tenors, spline), seed)
            quotes = CMTYield.(0.03, filter(>(0), tenors))
            loss(z) = sum(q -> (present_value(ZeroRateCurve(z, tenors, spline), q.instrument) - q.price)^2, quotes)
            @test all(isfinite, ForwardDiff.gradient(loss, seed))
            @test all(isfinite, ForwardDiff.hessian(loss, seed))
        end
    end

    @testset "fit passes solve_kwargs to the optimizer" begin
        qs = CMTYield.([0.03, 0.032, 0.035, 0.037], [1.0, 2.0, 3.0, 5.0])
        @test_throws FitConvergenceError fit(Spline.MonotoneConvex(), qs; solve_kwargs = (; maxiters = 1))
        @test_throws FitConvergenceError fit(Spline.Cubic(), qs, Fit.Loss(abs2); solve_kwargs = (; maxiters = 1))
        @test_throws FitConvergenceError fit(Yield.NelsonSiegel(), qs; solve_kwargs = (; maxiters = 1))
        # a solver-specific setting reaches Optim
        tight = fit(Spline.MonotoneConvex(), qs; solve_kwargs = (; g_tol = 1.0e-14))
        @test maximum(abs(present_value(tight, q.instrument) - q.price) for q in qs) < 1.0e-10
    end

    @testset "optimizer failures are never returned as fitted models" begin
        qs = ZCBPrice.([0.97, 0.93, 0.88], [1.0, 2.0, 3.0])
        fits = (
            () -> fit(Yield.Constant(), qs; optimizer = FailingFitOptimizer()),
            () -> fit(Spline.MonotoneConvex(), qs; optimizer = FailingFitOptimizer()),
            () -> fit(Spline.Linear(), qs, Fit.Loss(abs2); optimizer = FailingFitOptimizer()),
        )

        for run_fit in fits
            err = try
                run_fit()
                nothing
            catch e
                e
            end
            @test err isa FitConvergenceError
            @test err isa Exception
            if err isa FitConvergenceError
                @test err.retcode == FinanceModels.Optimization.SciMLBase.ReturnCode.Failure
                msg = sprint(showerror, err)
                @test startswith(msg, "FitConvergenceError: ")
                @test occursin("optimizer return code Failure", msg)
            end
        end
    end

    @testset "ZeroRateCurve negative time" begin
        zrc = ZeroRateCurve([0.02, 0.03], [1.0, 2.0])
        @test discount(zrc, 0.0) == 1.0
        @test_throws DomainError discount(zrc, -1.0)
    end

    @testset "bootstrap accepts unsorted quotes" begin
        sorted = fit(Spline.Linear(), ZCBPrice.([0.9, 0.8, 0.7], [1.0, 2.0, 3.0]), Fit.Bootstrap())
        unsorted = fit(Spline.Linear(), ZCBPrice.([0.8, 0.9, 0.7], [2.0, 1.0, 3.0]), Fit.Bootstrap())
        for t in 0.5:0.5:3.0
            @test discount(sorted, t) ≈ discount(unsorted, t)
        end
    end

    @testset "bootstrap rejects duplicate maturities loudly" begin
        # previously surfaced as a cryptic root-bracketing failure
        @test_throws ArgumentError fit(Spline.Linear(), ZCBPrice.([0.9, 0.89, 0.7], [1.0, 1.0, 3.0]), Fit.Bootstrap())
    end

    @testset "MonotoneConvex fit with a single quote" begin
        # range(0.01, 0.05, length=1) used to throw before any solving happened
        c = fit(Spline.MonotoneConvex(), [ZCBPrice(0.95, 2.0)])
        @test discount(c, 2.0) ≈ 0.95 atol = 1.0e-8
    end

    @testset "second-order optimizer fit does not warn about ADtype" begin
        # The generic and MonotoneConvex loss functions declared first-order
        # `AutoForwardDiff`, so a second-order optimizer made OptimizationBase warn and
        # auto-promote to `SecondOrder` once per fit. They now declare `SecondOrder` AD.
        # The bare `@test_logs` (default `min_level = Warn`) asserts no warnings; the
        # default `LBFGS` path only requests gradients, so it is unaffected.
        qs = ZCBPrice.([0.97, 0.93, 0.88, 0.8], [1.0, 2.0, 3.0, 5.0])
        reprice(c) = maximum(abs, present_value(c, q.instrument) - q.price for q in qs)

        # generic (bounded) path: IPNewton is the bounds-compatible second-order method
        # (plain Newton is rejected by Optim's Fminbox for box-constrained problems).
        ipn = FinanceModels.OptimizationOptimJL.IPNewton()
        cp = @test_logs fit(Yield.CairnsPritchard(), qs; optimizer = ipn)
        @test reprice(cp) < 1.0e-8
        @test_logs fit(Yield.NelsonSiegel(), qs; optimizer = ipn)  # smooth 4-param fit: assert only quiet

        # MonotoneConvex (unbounded) path: Newton is the canonical second-order method.
        mc = @test_logs fit(Spline.MonotoneConvex(), qs; optimizer = FinanceModels.OptimizationOptimJL.Newton())
        @test reprice(mc) < 1.0e-10
    end
end

@testset "refitting a flat PCHIP or Akima curve" begin
    # Refitting a knot curve ran the generic optic fit, whose candidates were public
    # `reconstruct`ions: those check dual knot rates for kinks, so the optimizer's own
    # derivatives at a flat PCHIP or Akima start threw. The refit is now the spline fit on the
    # curve's knots, whose trial curves are unchecked and whose start avoids the kinks.
    tenors = [1.0, 2.0, 3.0, 5.0, 10.0]
    qs = ZCBPrice.(exp.(-[0.03, 0.032, 0.035, 0.037, 0.04] .* tenors), tenors)
    for spline in (Spline.PCHIP(), Spline.Akima())
        c = fit(ZeroRateCurve(fill(0.03, 5), tenors, spline), qs)
        @test maximum(abs(pv(c, q.instrument) - q.price) for q in qs) < 1.0e-9
    end
end
