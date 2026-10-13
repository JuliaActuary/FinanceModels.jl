using ForwardDiff

@testset "quote_shift" begin
    tenors = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0]
    semiannual(y, t) = ParYield(y, t; frequency = 2)
    base = fit(Spline.MonotoneConvex(), semiannual.([0.045, 0.042, 0.04, 0.039, 0.04, 0.042, 0.046, 0.045], tenors))
    liability = Cashflow.(1.0, 0.5:0.5:25)
    bumps(c) = FinanceCore.rate.(knot_rates(c.r2))
    quotes(c, family, ts) = [implied_quote(c, family, t) for t in ts]
    function thrown(f)
        try
            f()
        catch e
            return e
        end
        return nothing
    end
    throws(f, needle) = (e = thrown(f); e isa ArgumentError && occursin(needle, sprint(showerror, e)))

    @testset "targets are met" begin
        families = (
            (y, t) -> ParYield(y, t; frequency = 1),
            semiannual,
            (y, t) -> ParSwapYield(y, t; frequency = 4),
            CMTYield,
            OISYield,
            ZCBYield,
        )
        for family in families
            q0 = quotes(base, family, tenors)
            shocked = quote_shift((q, t) -> q + 0.005, base, family, tenors)
            @test quotes(shocked, family, tenors) ≈ q0 .+ 0.005 atol = 1.0e-13
        end
        q0 = quotes(base, ZCBPrice, tenors)
        shocked = quote_shift((q, t) -> q * 0.99, base, ZCBPrice, tenors)
        @test quotes(shocked, ZCBPrice, tenors) ≈ q0 .* 0.99 rtol = 1.0e-13
        # every quoted par yield +100bp moves the zero rates by about as much
        up = quote_shift((q, t) -> q + 0.01, base, semiannual, tenors)
        @test all(b -> 0.0097 < b < 0.0102, bumps(up))
    end

    @testset "a zero shock returns the curve" begin
        for family in (semiannual, CMTYield, OISYield, ZCBPrice)
            same = quote_shift((q, t) -> q, base, family, tenors)
            @test maximum(abs(discount(same, t) - discount(base, t)) for t in 0.1:0.1:40) < 1.0e-14
        end
    end

    @testset "one tenor's bump leaves the others' quotes" begin
        q0 = quotes(base, CMTYield, tenors)
        for k in eachindex(tenors)
            shocked = quote_shift((q, t) -> t == tenors[k] ? q + 0.0001 : q, base, CMTYield, tenors)
            @test quotes(shocked, CMTYield, tenors) ≈ q0 .+ 0.0001 .* (eachindex(tenors) .== k) atol = 1.0e-13
        end
    end

    @testset "derivatives" begin
        # The quote Jacobian at a zero shock is the identity.
        J = ForwardDiff.jacobian(zeros(length(tenors))) do s
            shocked = quote_shift((q, t) -> q + s[findfirst(==(t), tenors)], base, semiannual, tenors)
            quotes(shocked, semiannual, tenors)
        end
        @test J ≈ [i == j for i in eachindex(tenors), j in eachindex(tenors)] atol = 1.0e-12

        # Against central differences, away from zero and at a shock whose value is zero.
        value(s) = pv(quote_shift((q, t) -> q + s, base, semiannual, tenors), liability)
        h = 1.0e-5
        for s in (0.01, 0.0)
            @test ForwardDiff.derivative(value, s) ≈ (value(s + h) - value(s - h)) / 2h rtol = 1.0e-7
        end
        @test ForwardDiff.derivative(value, 0.0) != 0

        # The value of a dual result is the primal result exactly.
        dual = quote_shift((q, t) -> q + ForwardDiff.Dual(0.0, 1.0), base, semiannual, tenors)
        same = quote_shift((q, t) -> q, base, semiannual, tenors)
        @test all(ForwardDiff.value(discount(dual, t)) == discount(same, t) for t in 0.25:0.25:40)

        # through the base curve's parameters
        flat(x) = Yield.Constant(Continuous(x))
        on_curve(x) = pv(quote_shift((q, t) -> q + 0.01, flat(x), semiannual, tenors), liability)
        @test ForwardDiff.derivative(on_curve, 0.03) ≈ (on_curve(0.03 + h) - on_curve(0.03 - h)) / 2h rtol = 1.0e-7

        # Partials carried only by a quote's price.
        price(s) = pv(quote_shift((q, t) -> q * (1 + s), base, ZCBPrice, tenors), liability)
        @test ForwardDiff.derivative(price, 0.0) ≈ (price(h) - price(-h)) / 2h rtol = 1.0e-7
        zcb(s) = pv(quote_shift((q, t) -> q + s, base, ZCBYield, tenors), liability)
        @test ForwardDiff.derivative(zcb, 0.0) ≈ (zcb(h) - zcb(-h)) / 2h rtol = 1.0e-7

        # Partials that a rule introduces only at the last tenor.
        last_only(s) = quote_shift((q, t) -> t == last(tenors) ? q + s : q, base, semiannual, tenors)
        db = ForwardDiff.derivative(s -> bumps(last_only(s)), 0.0)
        @test all(iszero, db[1:(end - 1)])
        @test last(db) ≈ (last(bumps(last_only(h))) - last(bumps(last_only(-h)))) / 2h rtol = 1.0e-7
        tail(s) = pv(last_only(s), liability)
        @test ForwardDiff.derivative(tail, 0.0) ≈ (tail(h) - tail(-h)) / 2h rtol = 1.0e-7
    end

    @testset "a shock and its reverse cancel" begin
        up = quote_shift((q, t) -> q + 0.01, base, semiannual, tenors)
        back = quote_shift((q, t) -> q - 0.01, up, semiannual, tenors)
        @test maximum(abs(discount(back, t) - discount(base, t)) for t in 0.1:0.1:40) < 1.0e-13
    end

    @testset "one tenor shifts in parallel" begin
        q0 = implied_quote(base, semiannual, 5.0)
        shocked = quote_shift((q, t) -> q + 0.01, base, semiannual, [5.0])
        @test implied_quote(shocked, semiannual, 5.0) ≈ q0 + 0.01 atol = 1.0e-13
        spreads = [FinanceCore.rate(zero(shocked, t)) - FinanceCore.rate(zero(base, t)) for t in 0.5:0.5:40]
        @test maximum(spreads) - minimum(spreads) < 1.0e-14
    end

    @testset "stubs and short maturities" begin
        for (family, ts) in (
                (semiannual, [0.25, 0.4, 0.75, 1.3, 2.7, 5.0]),
                (CMTYield, [0.25, 0.5, 1.0, 1.5, 2.0]),
                (OISYield, [0.5, 1.0, 2.0, 5.0]),
            )
            shocked = quote_shift((q, t) -> q + 0.002, base, family, ts)
            @test quotes(shocked, family, ts) ≈ quotes(base, family, ts) .+ 0.002 atol = 1.0e-13
        end
    end

    @testset "negative rates and coupons" begin
        negative = Yield.Constant(Continuous(-0.005))
        q0 = quotes(negative, semiannual, tenors)
        @test all(<(0), q0)
        shocked = quote_shift((q, t) -> q - 0.003, negative, semiannual, tenors)
        @test quotes(shocked, semiannual, tenors) ≈ q0 .- 0.003 atol = 1.0e-13
    end

    @testset "payments at time 0 count toward the floor" begin
        upfront(q, t) = Quote(q, [Cashflow(-0.5, 0.0), Cashflow(1.0, t)])
        flat = Yield.Constant(Continuous(0.03))
        shocked = quote_shift((q, t) -> q + 0.01, flat, upfront, [1.0, 2.0])
        @test quotes(shocked, upfront, [1.0, 2.0]) ≈ quotes(flat, upfront, [1.0, 2.0]) .+ 0.01 atol = 1.0e-13
        # The price must exceed -0.5, the value of the payment at time 0.
        @test throws(() -> quote_shift((q, t) -> -0.5, flat, upfront, [1.0, 2.0]), "does not exceed -0.5")
        shocked = quote_shift((q, t) -> -0.49, flat, upfront, [1.0, 2.0])
        @test quotes(shocked, upfront, [1.0, 2.0]) ≈ [-0.49, -0.49] atol = 1.0e-13
    end

    @testset "long tenors" begin
        # Doubling a 1000-year zero-coupon price needs a bump of -log(2)/1000, although a probe at
        # x = -1 overflows.
        flat = Yield.Constant(Continuous(0.03))
        shocked = quote_shift((q, t) -> 2q, flat, ZCBPrice, [1000.0])
        @test only(bumps(shocked)) ≈ -log(2) / 1000 rtol = 1.0e-12
    end

    @testset "an exact root at zero keeps its partials" begin
        flat = Yield.Constant(Continuous(0.03))
        ts = [1.0, 2.5, 30.0, 1000.0]
        # the residual at a zero bump is exactly zero for these quotes
        for t in ts
            q = ZCBPrice(implied_quote(flat, ZCBPrice, t), t)
            @test iszero(pv(flat + ZeroRateCurve([0.0], [t], Spline.Linear(); extrapolation = :flat_zero), q.instrument) - q.price)
        end
        # discount(t)·exp(-b·t) = discount(t)·(1 + s), so b = -log(1 + s)/t
        @test bumps(quote_shift((q, t) -> q, flat, ZCBPrice, ts)) == zeros(length(ts))
        @test ForwardDiff.derivative(s -> bumps(quote_shift((q, t) -> q * (1 + s), flat, ZCBPrice, ts)), 0.0) ≈ -1 ./ ts rtol = 1.0e-12
    end

    @testset "errors" begin
        @test throws(() -> quote_shift((q, t) -> q, base, semiannual, [0.0, 1.0]), "t = 0")
        @test throws(() -> ForwardDiff.derivative(s -> pv(quote_shift((q, t) -> q, base, semiannual, [1.0, 2.0 + s]), liability), 0.0), "tenors")
        # The grid is checked before any quote is built.
        @test throws(() -> quote_shift((q, t) -> error("rule called"), base, semiannual, [2.0, 1.0]), "quote_shift: tenors must be strictly increasing")
        @test throws(() -> quote_shift((q, t) -> error("rule called"), base, semiannual, [1.0, 1.0]), "quote_shift: tenors must be strictly increasing")

        # (S1) A quote paying after its tenor would also move with the next bump. Its target still
        # reprices at a zero shock, but its quote Jacobian is not the identity, so it is rejected,
        # also when the shock is a zero-valued dual number.
        beyond(q, t) = Quote(q, [Cashflow(1.0, t), Cashflow(1.0, t + 1)])
        @test throws(() -> quote_shift((q, t) -> q, base, beyond, [1.0, 2.0]), "after its tenor")
        @test throws(() -> ForwardDiff.derivative(s -> pv(quote_shift((q, t) -> q + s, base, beyond, [1.0, 2.0]), liability), 0.0), "after its tenor")
        # (S2) nothing to move
        early(q, t) = Quote(q, Cashflow(1.0, t / 4))
        @test throws(() -> quote_shift((q, t) -> q, base, early, [1.0, 2.0]), "pays nothing after the previous tenor")
        # (S3) amounts that change sign from positive to negative
        @test throws(() -> quote_shift((q, t) -> q, base, (q, t) -> Quote(q, [Cashflow(1.0, t - 0.5), Cashflow(-0.5, t)]), [1.0, 2.0]), "change sign at most once")
        # (S4) a +500% par yield at 30 years
        @test throws(() -> quote_shift((q, t) -> t == 30 ? q + 5.0 : q, base, semiannual, tenors), "unreachable")
        # A valid quote whose root (x ≈ -691181) lies past overflow: -1 at t - 0.001 and 1e-300 at t.
        @test throws(
            () -> quote_shift((q, t) -> 0.5, Yield.Constant(Continuous(0.03)), (q, t) -> Quote(q, [Cashflow(-1.0, t - 0.001), Cashflow(1.0e-300, t)]), [1.0]),
            "could not bracket"
        )
        @test throws(
            () -> ForwardDiff.derivative(u -> ForwardDiff.derivative(s -> pv(quote_shift((q, t) -> q + s + u, base, semiannual, tenors), liability), 0.0), 0.0),
            "first-order"
        )
    end

    @testset "the result is a CompositeYield" begin
        @test quote_shift((q, t) -> q, base, CMTYield, tenors) isa Yield.CompositeYield
    end
end
