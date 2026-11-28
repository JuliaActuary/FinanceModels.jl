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

    @testset "shift quotes, reprice them, measure bucket risk" begin
        # Stub tenors, a shock that differs by tenor in size and sign, and a liability paying
        # between the tenors and beyond the last.
        ts = [0.25, 0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0]
        owed = Cashflow.(1.0, 0.3:0.45:40)
        shock = [0.004 * sin(3k) + 0.001 for k in eachindex(ts)]
        at(t) = findfirst(==(t), ts)
        only_at(k) = eachindex(ts) .== k
        for family in (semiannual, (y, t) -> ParSwapYield(y, t; frequency = 4), CMTYield, OISYield, ZCBYield)
            # Unchanged quotes leave the liability's value.
            @test pv(quote_shift((q, t) -> q, base, family, ts), owed) ≈ pv(base, owed) rtol = 1.0e-13
            # Each shifted instrument reprices by `pv` on the result.
            shocked = quote_shift((q, t) -> q + shock[at(t)], base, family, ts)
            for t in ts
                target = family(implied_quote(base, family, t) + shock[at(t)], t)
                @test pv(shocked, target.instrument) ≈ target.price rtol = 1.0e-13
            end
            # Every bucket's sensitivity against that quote bumped up and down on its own, around
            # no shock and around the uneven one. The coupons carry the derivative on fixed dates.
            value(s) = pv(quote_shift((q, t) -> q + s[at(t)], base, family, ts), owed)
            h = 1.0e-5
            for s0 in (zero(shock), shock)
                bumped = [(value(s0 .+ h .* only_at(k)) - value(s0 .- h .* only_at(k))) / 2h for k in eachindex(ts)]
                # bucket by bucket: `≈` on the vectors would compare their norms
                @test all(isapprox.(ForwardDiff.gradient(value, s0), bumped; rtol = 1.0e-6))
            end
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

        # The bumps of a dual result have the values of the plain bumps exactly, here under a shock
        # that alternates in sign. (Discount factors between the tenors can differ in the last bit:
        # the interpolation takes another path through dual knot rates on an evenly spaced grid.)
        even = collect(3.0:3.0:36.0)
        alternating(s) = quote_shift((q, t) -> q * exp(-t * (isodd(Int(t / 3)) ? 0.005 : -0.005)) * (1 + s), base, ZCBPrice, even)
        dual, plain = alternating(ForwardDiff.Dual(0.0, 1.0)), alternating(0.0)
        @test ForwardDiff.value.(bumps(dual)) == bumps(plain)
        @test all(ForwardDiff.value(discount(dual, t)) ≈ discount(plain, t) for t in 3.0:0.01:36.0)

        # through the base curve's parameters
        flat(x) = Yield.Constant(Continuous(x))
        on_curve(x) = pv(quote_shift((q, t) -> q + 0.01, flat(x), semiannual, tenors), liability)
        @test ForwardDiff.derivative(on_curve, 0.03) ≈ (on_curve(0.03 + h) - on_curve(0.03 - h)) / 2h rtol = 1.0e-7

        # Partials carried only by a quote's price.
        price(s) = pv(quote_shift((q, t) -> q * (1 + s), base, ZCBPrice, tenors), liability)
        @test ForwardDiff.derivative(price, 0.0) ≈ (price(h) - price(-h)) / 2h rtol = 1.0e-7
        zcb(s) = pv(quote_shift((q, t) -> q + s, base, ZCBYield, tenors), liability)
        @test ForwardDiff.derivative(zcb, 0.0) ≈ (zcb(h) - zcb(-h)) / 2h rtol = 1.0e-7

        # A payment up to the tenor that is zero in value only keeps its derivative: the quote's
        # own solve carries it. The coupon at t/2 is `(q - q0)/2`, zero at a zero shock.
        q0 = [implied_quote(base, ZCBPrice, t) for t in tenors]
        early(q, t) = Quote(q, [Cashflow((q - q0[findfirst(==(t), tenors)]) / 2, t / 2), Cashflow(1.0, t)])
        coupon(s) = pv(quote_shift((q, t) -> q + s, base, early, tenors), liability)
        @test ForwardDiff.derivative(coupon, 0.0) ≈ (coupon(h) - coupon(-h)) / 2h rtol = 1.0e-7
        reprice = ForwardDiff.derivative(0.0) do s
            c = quote_shift((q, t) -> q + s, base, early, tenors)
            pv(c, early(q0[4] + s, tenors[4]).instrument)
        end
        @test reprice ≈ 1 rtol = 1.0e-12

        # Partials that a rule introduces only at the last tenor.
        last_only(s) = quote_shift((q, t) -> t == last(tenors) ? q + s : q, base, semiannual, tenors)
        db = ForwardDiff.derivative(s -> bumps(last_only(s)), 0.0)
        @test all(iszero, db[1:(end - 1)])
        @test last(db) ≈ (last(bumps(last_only(h))) - last(bumps(last_only(-h)))) / 2h rtol = 1.0e-7
        tail(s) = pv(last_only(s), liability)
        @test ForwardDiff.derivative(tail, 0.0) ≈ (tail(h) - tail(-h)) / 2h rtol = 1.0e-7
    end

    @testset "payment dates may not move under differentiation" begin
        # A unit payment at t(1 + q)/3, priced q: its date depends on the quote, and the quote on
        # the curve's rate. Without dual numbers the date is a number like any other. With them it
        # is rejected, also strictly between the tenors, away from every kink.
        dated(q, t) = Quote(q, Cashflow(1.0, t * (1 + q) / 3))
        level(r) = Yield.Constant(Continuous(r))
        same(r) = quote_shift((q, t) -> q, level(r), dated, [1.0, 2.0])
        @test bumps(same(0.1)) ≈ zeros(2) atol = 1.0e-13
        @test throws(() -> ForwardDiff.derivative(r -> bumps(same(r)), 0.1), "Payment dates must be fixed")

        # The reason: a date that moves across a tenor meets a kink in time there. Either direction
        # is a valid move with its own derivative (-0.6947 later, -0.6250 earlier, for this
        # two-year quote that costs 1.8 + s and pays at 1 + s and at 2), and ForwardDiff would take
        # the side from the seed's direction.
        crossing(q, t) = t == 1 ? ZCBPrice(q, t) : Quote(q, [Cashflow(1.0, q - 0.8), Cashflow(1.0, 2.0)])
        second(s) = bumps(quote_shift((q, t) -> t == 1 ? 1.0 : 1.8 + s, level(0.0), crossing, [1.0, 2.0]))[2]
        h = 1.0e-6
        @test (second(h) - second(0.0)) / h ≈ -0.6947 rtol = 1.0e-3
        @test (second(0.0) - second(-h)) / h ≈ -0.625 rtol = 1.0e-3
        @test throws(() -> ForwardDiff.derivative(second, 0.0), "Payment dates must be fixed")
        @test throws(() -> ForwardDiff.derivative(s -> second(-s), 0.0), "Payment dates must be fixed")
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

    @testset "portfolios" begin
        # A quote on several contracts: a vector of bonds, a tuple of cashflows, and one inside the other.
        flat = Yield.Constant(Continuous(0.03))
        ts = [1.0, 2.0, 5.0]
        families = (
            (q, t) -> Quote(q, [Bond.Fixed(0.03, Periodic(2), t)]),
            (q, t) -> Quote(q, [Bond.Fixed(0.03, Periodic(2), t), Bond.Fixed(0.05, Periodic(1), t)]),
            (q, t) -> Quote(q, (Cashflow(0.02, t / 2), Cashflow(1.0, t))),
            (q, t) -> Quote(q, ([Bond.Fixed(0.03, Periodic(2), t)], Cashflow(0.1, t))),
        )
        for family in families
            shocked = quote_shift((q, t) -> 0.98q, flat, family, ts)
            @test quotes(shocked, family, ts) ≈ 0.98 .* quotes(flat, family, ts) rtol = 1.0e-13
        end
        # An empty member adds no payments.
        hollow(q, t) = Quote(q, (Cashflow{Float64, Float64}[], Cashflow(1.0, t)))
        @test quotes(quote_shift((q, t) -> 0.98q, flat, hollow, ts), hollow, ts) ≈ 0.98 .* quotes(flat, hollow, ts) rtol = 1.0e-13
        # Integer amounts on one date are summed as floating-point numbers: two of 1.5e9 overflow an Int32.
        whole(q, t) = Quote(q, [Cashflow(Int32(1_500_000_000), t), Cashflow(Int32(1_500_000_000), t)])
        @test quotes(quote_shift((q, t) -> 0.98q, flat, whole, ts), whole, ts) ≈ 0.98 .* quotes(flat, whole, ts) rtol = 1.0e-13
        # The checks see every member's payments: one bond of the two matures after the tenor.
        longer(q, t) = Quote(q, [Bond.Fixed(0.03, Periodic(2), t), Bond.Fixed(0.03, Periodic(2), t + 1)])
        @test throws(() -> quote_shift((q, t) -> q, flat, longer, ts), "after its tenor")
    end

    @testset "payments that offset beyond the precision of the valuation" begin
        # +1e17 today, a nested leg of -1e17 and -0.5 today, and +1 in a year. The portfolio's
        # value rounds the 0.5 away, so its implied price is the discount factor, 1, where the
        # payments sum to 0.5. No bump is right for both, and the quote is rejected (S5), also
        # under a rule that returns the quote unchanged.
        par = Yield.Constant(Continuous(0.0))
        legs(q, t) = Quote(q, (Cashflow(1.0e17, 0.0), [Cashflow(-1.0e17, 0.0), Cashflow(-0.5, 0.0)], Cashflow(1.0, t)))
        @test implied_quote(par, legs, 1.0) == 1
        @test throws(() -> quote_shift((q, t) -> q, par, legs, [1.0]), "does not determine its bump")
        # A large payment the bump does not move, and a price it dominates: 2^53 today and 1 in a
        # year, at a price of 2^53 + 2. The year's discount factor is 2, but the value meets the
        # price exactly anywhere from 1 to 3 (the first probe, e, used to be returned).
        heavy(s, t) = Quote(2.0^53 + 2s, [Cashflow(2.0^53, 0.0), Cashflow(1.0, t)])
        @test pv(par + ZeroRateCurve([-1.0], [1.0], Spline.Linear(); extrapolation = :flat_zero), heavy(1.0, 1.0).instrument) == heavy(1.0, 1.0).price
        @test throws(() -> quote_shift((q, t) -> 1.0, par, heavy, [1.0]), "does not determine its bump")
        # Two payments on one date that offset each other hide a third from the sum by date,
        # though not from the value: 3.5 + 2^60 - 2^60 at one year, summed in that order, is 0.
        # With it the payments change sign twice, and a price of 1 has three roots (discount
        # factors of 1/2, 1 and 2 at one year); without it they pass (S3).
        hidden(a, t) = (Cashflow(a, t), (Cashflow(2.0^60, t), Cashflow(-2.0^60, t)))
        thrice(q, t) = Quote(q, (hidden(3.5, 1.0), Cashflow(-3.5, 2.0), Cashflow(1.0, t)))
        @test pv(par, thrice(1.0, 3.0).instrument) == 1
        for price in (1 - 1.0e-6, 1.0, 1 + 1.0e-6)
            @test throws(() -> quote_shift((q, t) -> price, par, thrice, [3.0]), "does not determine its bump")
        end
        # The same pair hides a payment of 1 today from the floor (S4): a price of 0.25 is below
        # it, and has two roots.
        sunk(q, t) = Quote(q, (hidden(1.0, 0.0), Cashflow(-2.0, 1.0), Cashflow(1.0, t)))
        @test throws(() -> quote_shift((q, t) -> 0.25, par, sunk, [2.0]), "does not determine its bump")
        # A payment just after the previous tenor hardly moves with the bump, whatever it is
        # worth: 1 at 8 + 2^-21 years responds to the 24-year bump 2.4e-7 times as much as to its
        # own rate. Beside 2^24 today, the value met its price at the first probe, 13% away from
        # the 24-year discount factor the price implies.
        near(q, t) = t == 8 ? ZCBPrice(q, t) : Quote(2.0^24 + q, [Cashflow(2.0^24, 0.0), Cashflow(1.0, 8 + 2.0^-21)])
        @test throws(() -> quote_shift((q, t) -> t == 8 ? 1.0 : 1 - 3 * eps(2.0^24), par, near, [8.0, 24.0]), "does not determine its bump")
        # Sizes at the ends of the floating-point range. Two payments of 2^1023 on one date sum
        # to `Inf`, against which any value passed for the price (twice it, here).
        vast(q, t) = Quote((1 + q) * 2.0^1023, [Cashflow(2.0^1023, t), Cashflow(2.0^1023, t)])
        @test throws(() -> quote_shift((q, t) -> -0.5, Yield.Constant(Continuous(log(2.0) / 10)), vast, [10.0]), "does not determine its bump")
        # Five of the smallest number, discounted by e, round to two of it: the price of two is
        # met at the first probe, where it implies a discount factor of 2/5.
        speck(q, t) = Quote((5 + q) * nextfloat(0.0), Cashflow(5 * nextfloat(0.0), t))
        @test throws(() -> quote_shift((q, t) -> -3.0, par, speck, [1.0]), "does not determine its bump")
        # Offsetting payments of ordinary size are resolved: a million each way around a unit.
        modest(q, t) = Quote(q, (Cashflow(1.0e6, 0.0), [Cashflow(-1.0e6, 0.0), Cashflow(-0.5, 0.0)], Cashflow(1.0, t)))
        @test quotes(quote_shift((q, t) -> q + 0.01, par, modest, [1.0]), modest, [1.0]) ≈ [0.51] rtol = 1.0e-8
    end

    @testset "a target price of zero" begin
        # Payments that offset one another: the quote's size is their gross value, 1.53 here, not
        # its price, so a target of zero is reached and checked like any other.
        level = Yield.Constant(Continuous(0.0))
        netted(q, t) = Quote(q, [Cashflow(-0.33, 0.0), Cashflow(0.2, t / 2), Cashflow(1.0, t)])
        @test implied_quote(level, netted, 1.0) ≈ 0.87 atol = 1.0e-14
        shocked = quote_shift((q, t) -> 0.0, level, netted, [1.0])
        @test pv(shocked, netted(0.0, 1.0).instrument) ≈ 0 atol = 1.0e-15
        # 0.2·exp(-x/2) + exp(-x) = 0.33, a quadratic in exp(-x/2)
        @test only(bumps(shocked)) ≈ -2log((-0.2 + sqrt(0.2^2 + 4 * 0.33)) / 2) rtol = 1.0e-13
        slope = ForwardDiff.derivative(s -> only(bumps(quote_shift((q, t) -> s, level, netted, [1.0]))), 0.0)
        @test slope ≈ (only(bumps(quote_shift((q, t) -> 1.0e-6, level, netted, [1.0]))) - only(bumps(quote_shift((q, t) -> -1.0e-6, level, netted, [1.0])))) / 2.0e-6 rtol = 1.0e-7
    end

    @testset "offsetting payments the bump does not move" begin
        # +1e6 today and -1e6 at one year cancel, and the two-year bump moves neither. The quote's
        # gross size is 2e6, but the bump's slope is that of the unit payment at two years: with
        # the one-year quote held, d(bump)/d(price) = -1/2.
        level = Yield.Constant(Continuous(0.0))
        offset(q, t) = t == 1 ? ZCBPrice(q, t) : Quote(q, [Cashflow(1.0e6, 0.0), Cashflow(-1.0e6, 1.0), Cashflow(1.0, 2.0)])
        second(s) = bumps(quote_shift((q, t) -> t == 2 ? q + s : q, level, offset, [1.0, 2.0]))[2]
        @test ForwardDiff.derivative(second, 0.0) ≈ -0.5 rtol = 1.0e-7
        @test second(1.0e-3) ≈ -log(1.001) / 2 rtol = 1.0e-6
    end

    @testset "long tenors" begin
        # Doubling a 1000-year zero-coupon price needs a bump of -log(2)/1000, although a probe at
        # x = -1 overflows.
        flat = Yield.Constant(Continuous(0.03))
        shocked = quote_shift((q, t) -> 2q, flat, ZCBPrice, [1000.0])
        @test only(bumps(shocked)) ≈ -log(2) / 1000 rtol = 1.0e-12
        # Large shocks at long tenors: the residual, relative to the quote's size before the shock,
        # is far below any absolute tolerance long before the root (the 1000-year target price is
        # about 4e-31), so the root is solved in its bracket and checked against the price itself.
        for (t, b) in ((1000.0, 0.04), (100.0, 0.4), (50.0, 0.5), (30.0, 0.5))
            shocked = quote_shift((q, t) -> q * exp(-b * t), flat, ZCBPrice, [t])
            @test only(bumps(shocked)) ≈ b rtol = 1.0e-13
            @test discount(shocked, t) ≈ exp(-(0.03 + b) * t) rtol = 1.0e-10
        end
        # A target price hundreds of orders above today's value: exp(-690) to exp(334). Both are
        # finite numbers, and so is the bump, but their ratio is not.
        steep = Yield.Constant(Continuous(0.69))
        raised = quote_shift((q, t) -> exp(334.0), steep, ZCBPrice, [1000.0])
        @test only(bumps(raised)) ≈ -1.024 rtol = 1.0e-13
        @test discount(raised, 1000.0) ≈ exp(334.0) rtol = 1.0e-10
        # And one hundreds of orders below it: exp(50) to exp(-700), a yield of -5% raised by 75%.
        # Measured against today's value, that target is nothing: a bump of 1.024, which discounts
        # the price to zero, used to pass for the root.
        negative = Yield.Constant(Continuous(-0.05))
        spread(s, t) = ZCBYield(Continuous(s - 0.05), t)
        lowered = quote_shift((q, t) -> 0.75, negative, spread, [1000.0])
        @test only(bumps(lowered)) ≈ 0.75 rtol = 1.0e-13
        @test discount(lowered, 1000.0) ≈ exp(-700.0) rtol = 1.0e-10
        # Its bump moves one for one with the quoted spread. (The slope of a residual divided by a
        # unit that carried the bump's derivative was not a number here: the unit's square is zero.)
        @test ForwardDiff.derivative(s -> only(bumps(quote_shift((q, t) -> s, negative, spread, [1000.0]))), 0.75) ≈ 1 rtol = 1.0e-12
        # A tenor so short that 1/t overflows its type: the first probe is the largest finite number.
        tiny = quote_shift((q, t) -> 0.999, Yield.Constant(Continuous(0.0)), ZCBPrice, Float16[1.0e-5])
        @test discount(tiny, Float16(1.0e-5)) ≈ 0.999 rtol = 1.0e-3
        # A probe that lands on the root exactly is a bracket end.
        @test only(bumps(quote_shift((q, t) -> q * exp(-1.0), flat, ZCBPrice, [1.0]))) ≈ 1.0 rtol = 1.0e-14
    end

    @testset "partials far from the quote's size" begin
        # The derivative step divides the price's partials by the value's slope, not by the
        # quote's unit first. A partial of 1e308 on a price of 0.5 is finite (it overflowed to
        # `Inf` over the unit), and one of 5e-324 on a tenor of 5e-324 years survives (it
        # rounded to zero).
        half = Yield.Constant(Continuous(log(2) / 8))
        @test ForwardDiff.derivative(s -> only(bumps(quote_shift((q, t) -> 0.5 + 1.0e308 * s, half, ZCBPrice, [8.0]))), 0.0) ≈ -2.5e307 rtol = 1.0e-12
        tiny = nextfloat(0.0)
        four(q, t) = Quote(q, Cashflow(4.0, t))
        @test ForwardDiff.derivative(s -> discount(quote_shift((q, t) -> 4.0 + tiny * s, Yield.Constant(Continuous(0.0)), four, [tiny]), 1.0), 0.0) ≈ 0.25 rtol = 1.0e-12
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
        # (S1) also for a payment of nothing after the tenor. Under differentiation its derivative
        # with respect to the shock is not zero, and the two-year bump would move it after the
        # one-year solve: the one-year quote's own sensitivity would come out as 0.9, not 1.
        late(q, t) = t == 1 ? Quote(q, [Cashflow(1.0, 1.0), Cashflow((q - 1) / 2, 2.0)]) : ZCBPrice(q, t)
        late_rule(s) = (q, t) -> t == 1 ? 1 + s : 0.8
        level = Yield.Constant(Continuous(0.0))
        @test throws(() -> ForwardDiff.derivative(s -> discount(quote_shift(late_rule(s), level, late, [1.0, 2.0]), 1.5), 0.0), "after its tenor")
        @test throws(() -> quote_shift(late_rule(0.0), level, late, [1.0, 2.0]), "after its tenor")
        @test throws(() -> quote_shift(late_rule(0.01), level, late, [1.0, 2.0]), "after its tenor")
        # Nor is a plain zero harmless there: a later bump can raise a discount factor far beyond
        # the tenor to `Inf` (here exp(800)), and the quote's value would be 0·Inf.
        padded(q, t) = t == 1 ? Quote(q, [Cashflow(1.0, t), Cashflow(0.0, 40_000.0)]) : ZCBPrice(q, t)
        @test throws(() -> quote_shift((q, t) -> t == 1 ? 1.0 : exp(0.04), level, padded, [1.0, 2.0]), "after its tenor")
        # An error reports the target's value, not its dual number.
        @test !occursin("Dual", sprint(showerror, thrown(() -> ForwardDiff.derivative(s -> discount(quote_shift(late_rule(s), level, late, [1.0, 2.0]), 1.5), 0.0))))
        # A payment date that moves with the shock is rejected. Without dual numbers the date is
        # a number, here the tenor itself.
        moving(q, t) = t == 1 ? Quote(q, Cashflow(1.0, (1 + q) / 2)) : ZCBPrice(q, t)
        moving_rule(s) = (q, t) -> t == 1 ? 1 + s : 0.8
        @test throws(() -> ForwardDiff.derivative(s -> discount(quote_shift(moving_rule(s), level, moving, [1.0, 2.0]), 1.5), 0.0), "Payment dates must be fixed")
        @test bumps(quote_shift(moving_rule(0.0), level, moving, [1.0, 2.0])) ≈ [0.0, -log(0.8) / 2] atol = 1.0e-14
        # A payment of nothing may sit on a moving date, inside the tenor on both sides.
        idle(q, t) = t == 1 ? Quote(q, [Cashflow(1.0, 0.5), Cashflow(0.0, (1 + q) / 3)]) : ZCBPrice(q, t)
        @test ForwardDiff.derivative(s -> bumps(quote_shift(moving_rule(s), level, idle, [1.0, 2.0]))[1], 0.0) ≈ -2 rtol = 1.0e-12
        # (S1) is checked one payment at a time. Two payments after the tenor that offset each
        # other are still paid after it, where the next bump moves both.
        paired(q, t) = t == 1 ? Quote(q, [Cashflow(1.0, 1.0), Cashflow(1.0, 2.0), Cashflow(-1.0, 2.0)]) : ZCBPrice(q, t)
        @test throws(() -> quote_shift(moving_rule(0.0), level, paired, [1.0, 2.0]), "after its tenor")
        # Summed first they would pass for nothing, and the quote's value would not survive the
        # next solve: +2^53 at two years, and a nested leg of -2^53 at two years and 1 at one. At
        # par the leg is -2^53 + 1, exactly, and the quote is worth 1. Once the two-year discount
        # factor is 33/32 the leg's 1 is rounded away, and the quote is worth 0.
        big = 2.0^53
        nested(q, t) = t == 1 ? Quote(q, (Cashflow(big, 2.0), [Cashflow(-big, 2.0), Cashflow(1.0, 1.0)])) : ZCBPrice(q, t)
        @test pv(level, nested(1.0, 1.0).instrument) == 1
        @test throws(() -> quote_shift((q, t) -> t == 1 ? q : 33 / 32, level, nested, [1.0, 2.0]), "after its tenor")
        # A target the valuation cannot represent: -2^54 today, and a nested leg of 2^54 at one
        # year and 1 at two. With the one-year discount factor held at 1, the leg's value is 2^54
        # or 2^54 + 4, and the quote's 0 or 4, never the target of 1. (Measured against the gross
        # size of its payments, a two-year discount factor of 2 used to pass for the root.)
        coarse(q, t) = t == 1 ? ZCBPrice(q, t) : Quote(q, (Cashflow(-2.0^54, 0.0), [Cashflow(2.0^54, 1.0), Cashflow(1.0, 2.0)]))
        @test throws(() -> quote_shift((q, t) -> 1.0, level, coarse, [1.0, 2.0]), "does not reprice")
        # (S5) A target of 4 it can represent, but for any two-year discount factor from 2 to 6.
        @test throws(() -> quote_shift((q, t) -> t == 1 ? 1.0 : 4.0, level, coarse, [1.0, 2.0]), "does not determine its bump")
        # (S2) nothing to move
        early(q, t) = Quote(q, Cashflow(1.0, t / 4))
        @test throws(() -> quote_shift((q, t) -> q, base, early, [1.0, 2.0]), "pays nothing after the previous tenor")
        # (S3) amounts that change sign from positive to negative
        @test throws(() -> quote_shift((q, t) -> q, base, (q, t) -> Quote(q, [Cashflow(1.0, t - 0.5), Cashflow(-0.5, t)]), [1.0, 2.0]), "change sign at most once")
        # (S4) At 25% the 30-year discount factor is 0.0006, so the 50-year par yield can sit at
        # most 2bp above the 30-year one. A shock asking for 80bp has no curve.
        steep = Yield.Constant(Continuous(0.25))
        @test throws(() -> quote_shift((q, t) -> t == 50 ? q + 0.005 : q - 0.003, steep, semiannual, [30.0, 50.0]), "unreachable")
        # (S4) a +500% par yield at 30 years
        @test throws(() -> quote_shift((q, t) -> t == 30 ? q + 5.0 : q, base, semiannual, tenors), "unreachable")
        # A valid quote whose root (x ≈ -691181) lies past overflow: -1 at t - 0.001 and 1e-300 at t.
        @test throws(
            () -> quote_shift((q, t) -> 0.5, Yield.Constant(Continuous(0.03)), (q, t) -> Quote(q, [Cashflow(-1.0, t - 0.001), Cashflow(1.0e-300, t)]), [1.0]),
            "could not bracket"
        )
        # A residual that keeps its sign and stays finite at every probe: the probes stop at the
        # largest finite one. (They used to double to `Inf` and repeat it without end.)
        probes = Ref(0)
        @test throws(() -> Yield.__quote_bracket(x -> (probes[] += 1; probes[] > 5000 && error("still probing at $x"); 5.0e-18), 1.0), "could not bracket")
        @test probes[] < 1100
        # also from the negative side, and where the last probes are too large to average by their sum
        probes[] = 0
        @test throws(() -> Yield.__quote_bracket(x -> (probes[] += 1; probes[] > 5000 && error("still probing at $x"); -5.0e-18), 1.0), "could not bracket")
        @test Yield.__quote_bracket(x -> x < 1.0e308 ? 1.0 : x < 1.5e308 ? -1.0 : Inf, 1.0) == (2.0^1023, 2.0^1022 + floatmax(Float64) / 2)
        @test throws(
            () -> ForwardDiff.derivative(u -> ForwardDiff.derivative(s -> pv(quote_shift((q, t) -> q + s + u, base, semiannual, tenors), liability), 0.0), 0.0),
            "first-order"
        )
    end

    @testset "the result is a CompositeYield" begin
        @test quote_shift((q, t) -> q, base, CMTYield, tenors) isa Yield.CompositeYield
    end
end
