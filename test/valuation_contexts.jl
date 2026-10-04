module ValuationContextTests

using Test, FinanceCore, FinanceModels, Transducers, ForwardDiff

# A closed form written against the context: it reads its index curve by key.
struct OnePeriodFloater <: FinanceCore.AbstractContract
    key::Symbol
end
FinanceCore.present_value(ctx, c::OnePeriodFloater) = discount(ctx, 2.0) / discount(ctx[c.key], 1.0, 2.0)
FinanceCore.maturity(::OnePeriodFloater) = 2.0

# A contract defined by its projection (`asfoldable`).
struct TwoFlows <: FinanceCore.AbstractContract end
Transducers.asfoldable(::Projection{TwoFlows}) = [Cashflow(1.0, 1.0), Cashflow(2.0, 2.0)]

@testset "valuation contexts" begin
    curve = Yield.Constant(Continuous(0.04))
    credit = Yield.Constant(Continuous(0.06))
    fixed = Bond.Fixed(0.05, Periodic(2), 5.0)
    floating = Bond.Floating(0.001, Periodic(4), 5.0, :index)
    swap = InterestRateSwap(curve, 5.0; frequency = 4, model_key = :index)
    store = Dict(:index => curve)
    ctx = Models(credit, store)

    @testset "Models" begin
        every = Models(credit; index = curve)
        @test ctx[:index] === curve
        @test every[:anything] === curve
        @test discount(ctx, 3.0) == discount(credit, 3.0)
        @test discount(ctx, 1.0, 3.0) == discount(credit, 1.0, 3.0)
        @test valuation_model(ctx) === credit
        @test valuation_model(credit) === credit
        for c in (floating, swap, floating |> Map(-) |> Map(cf -> cf * 2), Forward(1.0, floating))
            @test present_value(ctx, c) == present_value(every, c)
            @test collect(Projection(c, ctx)) == collect(Projection(c, store))
            # derivatives through the index curve
            auto(r) = present_value(Models(credit; index = Yield.Constant(Continuous(r))), c)
            explicit(r) = present_value(Models(credit, Dict(:index => Yield.Constant(Continuous(r)))), c)
            @test ForwardDiff.derivative(auto, 0.04) ≈ ForwardDiff.derivative(explicit, 0.04)
        end
        @test present_value(Models(curve; index = curve), swap) ≈ 0.0 atol = 1.0e-12
        @test present_value(Models(credit, (index = curve,)), floating) == present_value(ctx, floating)
        # a floater reads its index from the context: a plain curve holds none
        @test_throws MethodError present_value(curve, floating)
    end

    @testset "linearity: collections of self-timed contracts" begin
        a, b = Bond.Fixed(0.0, Periodic(1), 3.0), Bond.Fixed(0.0, Periodic(1), 5.0)
        r = Continuous(0.03)
        # an index was read as a valuation time: [a] was worth 0.941765 and [a, b] ≠ [b, a]
        @test present_value(r, [a]) === present_value(r, a)
        @test present_value(r, [a, b]) ≈ present_value(r, [b, a]) ≈ present_value(r, a) + present_value(r, b)
        @test present_value(r, FinanceCore.Composite(a, b)) ≈ present_value(r, a) + present_value(r, b)
        portfolio = [swap, fixed, Forward(1.0, floating)]
        @test present_value(ctx, portfolio) ≈ sum(present_value(ctx, c) for c in portfolio)
        @test present_value(ctx, reverse(portfolio)) ≈ present_value(ctx, portfolio)
        @test present_value(ctx, FinanceCore.Composite(swap, FinanceCore.Composite(fixed, floating))) ≈
            present_value(ctx, swap) + present_value(ctx, fixed) + present_value(ctx, floating)
        # a projected contract is worth its projected cashflows
        for c in (
                fixed, floating, swap, Forward(1.5, FinanceCore.Composite(fixed, floating)),
                TwoFlows(), TwoFlows() |> Map(cf -> Cashflow(-cf.amount, cf.time)),
            )
            @test present_value(ctx, c) ≈ present_value(credit, collect(Projection(c, ctx))) rtol = 1.0e-14 atol = 1.0e-14
        end
        # a cashflow before time 0 accumulates (a valuation-time filter used to drop it)
        @test present_value(curve, Forward(-2.0, Cashflow(1.0, 1.0))) ≈ 1 / discount(curve, 1.0)
        # empty: an exact zero of the valuation's type
        m32 = Yield.Constant(Continuous(0.03f0))
        @test present_value(m32, FinanceCore.AbstractContract[]) === 0.0f0
        @test present_value(m32, TwoFlows() |> Filter(cf -> false)) === 0.0f0
        @test present_value(Yield.Constant(Continuous(big"0.03")), TwoFlows() |> Filter(cf -> false)) isa BigFloat
    end

    @testset "closed forms compose" begin
        bsm = Equity.BlackScholesMerton(0.01, 0.02, 0.15)
        call = Option.EuroCall(CommonEquity(), 1.0, 1.0)
        put = Option.EuroPut(CommonEquity(), 1.0, 1.0)
        # these threw: a portfolio read an index as a valuation time, and a composite folded cashflows
        @test present_value(bsm, [call, put]) ≈ present_value(bsm, call) + present_value(bsm, put)
        @test present_value(bsm, FinanceCore.Composite(call, put)) ≈ present_value(bsm, call) + present_value(bsm, put)
        hw = ShortRate.HullWhite(0.1, 0.01, curve)
        cap = Option.Cap(0.03, 4, 3.0)
        options = (cap, Option.Floor(0.03, 4, 3.0), Option.Swaption(1.0, 2.0, 0.011, 4), Option.ZCBCall(1.0, 2.0, 0.95), Option.ZCBPut(1.0, 2.0, 0.95))
        for c in options
            # the context's model prices it, whatever else the context holds
            @test present_value(Models(hw; index = curve), c) == present_value(hw, c)
        end
        hwctx = Models(hw; index = curve)
        @test present_value(hwctx, FinanceCore.Composite(cap, floating)) ≈ present_value(hw, cap) + present_value(hwctx, floating)
        @test present_value(hwctx, [cap, floating]) ≈ present_value(hw, cap) + present_value(hwctx, floating)
        # a model that does not price the contract fails in its formula
        @test_throws MethodError present_value(curve, cap)
        @test_throws MethodError present_value(Models(curve; index = curve), call)
        # a custom closed form written against the context composes too
        floater = OnePeriodFloater(:index)
        @test present_value(ctx, floater) ≈ discount(credit, 2.0) / discount(curve, 1.0, 2.0)
        @test present_value(ctx, [floater, fixed]) ≈ present_value(ctx, floater) + present_value(ctx, fixed)
        @test present_value(ctx, FinanceCore.Composite(floater, fixed)) ≈ present_value(ctx, floater) + present_value(ctx, fixed)
        # transformers act on cashflow streams, and a closed-form-only contract has none
        @test_throws MethodError present_value(bsm, Forward(1.0, call))
    end

    @testset "value as of t: an explicit reduction" begin
        asof(ctx, c, t) = foldxl(
            +, Projection(c, ctx) |> Filter(cf -> cf.time >= t) |> Map(cf -> cf.amount * discount(ctx, t, cf.time));
            init = zero(discount(ctx, t, t))
        )
        bond = Bond.Fixed(0.05, Periodic(1), 3.0)
        @test asof(curve, bond, 1.0) ≈ 0.05 + 0.05 * exp(-0.04) + 1.05 * exp(-0.08)
        @test asof(curve, bond, 0.0) ≈ present_value(curve, bond)
        @test asof(curve, bond, 5.0) === 0.0
        # each retained cashflow takes its interval discount, which is safe where an accumulation
        # times the time-0 value overflows (Inf * 0)
        r1 = Yield.Constant(Continuous(1.0))
        far = Cashflow(1.0, 1000.0)
        @test asof(r1, far, 1000.0) == 1.0
        @test isnan(accumulation(r1, 1000.0) * present_value(r1, far))
        @test asof(Yield.Constant(Continuous(0.03f0)), bond, 5.0f0) === 0.0f0
        @test asof(Yield.Constant(Continuous(big"0.03")), bond, 5.0) isa BigFloat
    end

    @testset "FX values: one reporting currency per context" begin
        eurusd = FX.Pair(:EUR, :USD)
        usd = Yield.Constant(Continuous(0.05))
        fx_explicit = FX.Forwards(eurusd, 1.1, usd, Yield.Constant(Continuous(0.03)) + Yield.Constant(Continuous(-0.002)))
        outrights = FX.Outright.(eurusd, [1.1055, 1.1113, 1.1225, 1.1459], [0.25, 0.5, 1.0, 2.0])
        fx_absorbed = fit(FX.Forwards(eurusd, 1.1, usd, Spline.Linear()), outrights, Fit.Bootstrap())
        leg = FX.ParBasisSwap(eurusd, -0.0015, 3.0; reference = Yield.Constant(Continuous(0.03))).instrument
        fwd = FX.Forward(eurusd, 1.2, 2.0)
        conv(c) = FX.Converted(c, eurusd, :fx)
        unconverted = "is denominated in :EUR"

        # an FX model discounts in its reporting (quote) currency, and an FX.Forward is one
        # quote-currency cashflow
        @test discount(fx_explicit, 2.0) == discount(usd, 2.0)
        @test present_value(fx_explicit, fwd) ≈ (forward(fx_explicit, 2.0) - 1.2) * discount(usd, 2.0) rtol = 1.0e-14
        @test collect(Projection(fwd, fx_explicit)) == [Cashflow(forward(fx_explicit, 2.0) - 1.2, 2.0)]
        @test_throws ArgumentError present_value(inv(fx_explicit), fwd)   # a crossed pair
        @test_throws MethodError present_value(usd, fwd)                  # a curve is not an FX model

        # the acceptance matrix: a native curve values the leg, through wrappers as well
        native = present_value(fx_explicit.foreign, leg)
        @test native ≈ present_value(fx_explicit.foreign, leg.cashflows)
        @test present_value(fx_explicit.foreign, leg |> Map(identity)) ≈ native
        @test present_value(fx_explicit.foreign, Forward(0.0, leg)) ≈ native
        @test present_value(fx_explicit.foreign, Forward(0.5, leg)) ≈
            sum(cf.amount * discount(fx_explicit.foreign, cf.time + 0.5) for cf in leg.cashflows)
        # a EUR leg throws under the USD-reporting FX model, Models over it, and Models over a plain curve
        for c in (fx_explicit, Models(fx_explicit, Dict(:fx => fx_explicit)), Models(usd, Dict(:fx => fx_explicit)))
            for form in (
                    leg, FinanceCore.AbstractContract[leg], FinanceCore.Composite(leg, Forward(1.0, leg)),
                    Forward(0.5, leg), Forward(0.0, leg), leg |> Map(cf -> Cashflow(-cf.amount, cf.time)),
                )
                @test_throws unconverted present_value(c, form)
            end
        end
        for form in (FinanceCore.AbstractContract[fwd, leg], FinanceCore.Composite(fwd, leg))
            @test_throws unconverted present_value(fx_explicit, form)
        end
        # FX.Converted converts each cashflow at its forward, through wrappers inside or outside it
        for fx in (fx_explicit, fx_absorbed)
            uctx = Models(usd, Dict(:fx => fx))
            fxctx = Models(fx, Dict(:fx => fx))
            v = present_value(uctx, conv(leg))
            @test v ≈ fx.spot * present_value(fx.foreign, leg) rtol = 1.0e-12
            @test present_value(fxctx, conv(leg)) ≈ v rtol = 1.0e-12
            @test present_value(uctx, conv(leg |> Map(identity))) ≈ v rtol = 1.0e-12
            @test present_value(uctx, conv(Forward(0.5, leg))) ≈ present_value(uctx, Forward(0.5, conv(leg))) rtol = 1.0e-12
            @test present_value(fxctx, FinanceCore.Composite(fwd, conv(leg))) ≈ present_value(fx, fwd) + v rtol = 1.0e-12
            @test present_value(fxctx, FinanceCore.AbstractContract[fwd, conv(leg)]) ≈ present_value(fx, fwd) + v rtol = 1.0e-12
        end
        # an already converted leg inside FX.Converted throws: the conversion multiplies every amount by
        # its forward, so it would be converted again
        fxctx = Models(fx_explicit, Dict(:fx => fx_explicit))
        @test_throws "pays in :USD, but its context is in :EUR" present_value(fxctx, conv(conv(leg)))
        # an FX.Forward prices on the context's own FX model, which a conversion does not have (a mixed
        # composite gave 1.0382 where converting only the leg gives 1.0454)
        for form in (
                conv(fwd), conv(FinanceCore.Composite(fwd, leg)), conv(FinanceCore.AbstractContract[leg, fwd]),
                conv(Forward(0.5, fwd)), conv(fwd |> Map(identity)),
            )
            @test_throws MethodError present_value(fxctx, form)
        end
        # a conversion into JPY under a USD context throws before the forward is reached (this
        # discounted JPY on the USD curve: 11.677 where the JPY value is 12.879)
        usdjpy = FX.Pair(:USD, :JPY)
        fx_uj = FX.Forwards(usdjpy, 150.0, Yield.Constant(Continuous(0.001)), usd)
        @test_throws "pays in :JPY, but its context is in :USD" present_value(Models(fx_explicit, Dict(:uj => fx_uj)), FX.Converted(fwd, usdjpy, :uj))
        # a chain of conversions, each from its own pair's base currency, converts GBP to EUR to USD
        gbpeur = FX.Pair(:GBP, :EUR)
        fx_gbp = FX.Forwards(gbpeur, 1.15, fx_explicit.foreign, Yield.Constant(Continuous(0.04)))
        chain = FX.Converted(FX.Converted(Cashflow(100.0, 2.0), gbpeur, :gbp), eurusd, :fx)
        @test present_value(Models(usd, Dict(:fx => fx_explicit, :gbp => fx_gbp)), chain) ≈
            100.0 * forward(fx_gbp, 2.0) * forward(fx_explicit, 2.0) * discount(usd, 2.0) rtol = 1.0e-12
        # under an FX model's context a conversion must pay in its reporting currency: a GBP→EUR
        # conversion under the USD context was discounted on the USD curve (1.0643, where GBP→EUR→USD
        # is 1.2185), also beside a USD-paying forward
        store = Dict(:fx => fx_explicit, :gbp => fx_gbp)
        gbp_eur = FX.Converted(Cashflow(100.0, 2.0), gbpeur, :gbp)
        @test_throws "pays in :EUR, but its context is in :USD" present_value(Models(fx_explicit, store), gbp_eur)
        @test_throws "pays in :EUR, but its context is in :USD" present_value(Models(fx_explicit, store), FinanceCore.Composite(fwd, gbp_eur))
        @test present_value(Models(fx_explicit, store), chain) ≈
            100.0 * forward(fx_gbp, 2.0) * forward(fx_explicit, 2.0) * discount(usd, 2.0) rtol = 1.0e-12
        # a leg that pays in the context's reporting currency is valued in it
        usd_leg = FX.ParBasisSwap(usdjpy, -0.001, 3.0; reference = usd).instrument
        @test present_value(fx_explicit, usd_leg) ≈ present_value(fx_explicit.domestic, usd_leg) rtol = 1.0e-14
        @test present_value(Models(fx_explicit, store), usd_leg) ≈ present_value(fx_explicit.domestic, usd_leg) rtol = 1.0e-14
        # cashflows can still be listed without a model
        @test length(collect(Projection(leg))) == length(leg.cashflows)

        # inverting the model keeps the economic trade: a long EURUSD forward struck at K is −K units
        # of the USDEUR forward struck at 1/K (the reciprocal trade alone is a different one)
        for fx in (fx_explicit, fx_absorbed), (K, T) in ((1.2, 2.0), (0.95, 1.5))
            eur_value = present_value(fx, FX.Forward(eurusd, K, T)) / fx.spot
            @test eur_value ≈ -K * present_value(inv(fx), FX.Forward(inv(eurusd), 1 / K, T)) rtol = 1.0e-12
            @test !(eur_value ≈ present_value(inv(fx), FX.Forward(inv(eurusd), 1 / K, T)))
        end

        # par basis-swap calibration prices its native quotes on the fitted foreign curve
        estr = Yield.Constant(0.03)
        quotes = [
            FX.Outright.(eurusd, [1.1055, 1.1113], [0.25, 0.5]);
            FX.ParBasisSwap.(eurusd, [-0.0012, -0.0018], [2.0, 5.0]; reference = estr)
        ]
        m = fit(FX.Forwards(eurusd, 1.1, Yield.Constant(0.05), Spline.Linear()), quotes, Fit.Bootstrap())
        for q in quotes
            if q.instrument isa FX.BasisSwapLeg
                @test present_value(m.foreign, q.instrument) ≈ q.price atol = 1.0e-10
            else
                @test present_value(m, q.instrument) ≈ q.price atol = 1.0e-10
            end
        end
    end
end

end
