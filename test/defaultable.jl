using ForwardDiff

@testset "defaultable contracts" begin
    r, λ, R = 0.04, 0.02, 0.4
    curve = Yield.Constant(Continuous(r))
    credit(h, rec) = Credit.Assumptions(h; recovery = rec)
    ctx(h, rec; c = curve) = Models(c, Dict("ACME" => credit(h, rec)))
    flat = ctx(Credit.ConstantHazard(λ), R)
    risky(a) = Credit.Defaultable(a, "ACME")

    # At constant r and λ: Σ cᵢ e^{-(r+λ)tᵢ} + R·λ·Σ Bᵢ₋₁ ∫ e^{-(r+λ)u} du over each payment interval.
    function analytic(asset, r, λ, R)
        rows = Bond.__schedule(asset)
        k = r + λ
        total, t_prev = 0.0, 0.0
        for row in rows
            total += (row.interest + row.principal) * exp(-k * row.time)
            total += R * λ * row.opening * (exp(-k * t_prev) - exp(-k * row.time)) / k
            t_prev = row.time
        end
        return total
    end

    @testset "closed forms" begin
        for asset in (
                Bond.Fixed(0.05, Periodic(2), 10), Bond.Fixed(0.0, Periodic(1), 7), Bond.Fixed(0.05, Periodic(2), 2.3),
                Bond.Amortizing(0.05, [1.0, 2.0, 3.0, 4.0], [25.0, 25.0, 25.0, 25.0]), Bond.LevelPayment(0.06, 12, 15),
            )
            @test present_value(flat, risky(asset)) ≈ analytic(asset, r, λ, R) rtol = 1.0e-14
        end
        @test maturity(risky(Bond.Fixed(0.05, Periodic(2), 10))) == 10
    end

    @testset "limits" begin
        b = Bond.Fixed(0.05, Periodic(2), 10)
        # no default risk is the bond itself
        @test present_value(ctx(Credit.ConstantHazard(0.0), R), risky(b)) ≈ present_value(curve, b)
        # no recovery is the survival-weighted payments
        @test present_value(ctx(Credit.ConstantHazard(λ), 0.0), risky(b)) ≈ sum(cf.amount * exp(-(r + λ) * cf.time) for cf in collect(b))
        # full recovery of principal is not default-free: it loses coupons above the curve...
        full = ctx(Credit.ConstantHazard(0.05), 1.0)
        @test present_value(full, risky(b)) < present_value(curve, b)
        # ...and gains on a bond below it, which it repays early at par
        low = Bond.Fixed(0.01, Periodic(2), 10)
        @test present_value(full, risky(low)) > present_value(curve, low)
        # the value falls as default becomes more likely, for this bond and recovery
        @test present_value(ctx(Credit.ConstantHazard(0.05), R), risky(b)) < present_value(flat, risky(b))
    end

    @testset "quadrature against a fine reference" begin
        # a piecewise hazard and a curve whose forward has kinks: compare with Simpson's rule on a fine grid
        h = Credit.PiecewiseHazard([0.01, 0.03, 0.02], [1.5, 4.25, 8.0])
        zc = ZeroRateCurve([0.02, 0.035, 0.03, 0.045], [1.0, 3.0, 7.0, 12.0])
        c = ctx(h, R; c = zc)
        b = Bond.Fixed(0.05, Periodic(4), 10)
        payments = sum(cf.amount * Credit.survival(h, cf.time) * discount(zc, cf.time) for cf in collect(b))
        # Simpson's rule loses its order across the hazard's jumps, so the reference splits there
        reference = 0.0
        cuts = [0.0, 1.5, 4.25, 8.0, 10.0]
        for (a, z) in zip(cuts[1:(end - 1)], cuts[2:end])
            m = 100_000
            xs = range(a, z; length = m + 1)
            g(u) = R * Credit.hazard(h, (a + z) / 2) * Credit.survival(h, u) * discount(zc, u)
            reference += (z - a) / (3m) * sum((i == 0 || i == m ? 1 : (isodd(i) ? 4 : 2)) * g(xs[i + 1]) for i in 0:m)
        end
        @test present_value(c, risky(b)) ≈ payments + reference rtol = 1.0e-11
    end

    @testset "projection" begin
        b = Bond.Fixed(0.05, Periodic(2), 2)
        cfs = collect(Projection(risky(b), flat))
        @test issorted([cf.time for cf in cfs])
        # four payments and, before each, 6 monthly parts of 4 nodes of expected recovery
        @test length(cfs) == 4 + 4 * 6 * 4
        @test sum(cf.amount * discount(curve, cf.time) for cf in cfs) ≈ present_value(flat, risky(b))
        # the expected payments are the payments times survival
        pays = [cf for cf in cfs if cf.time in (0.5, 1.0, 1.5, 2.0)]
        @test [cf.amount for cf in pays] ≈ [cf.amount * exp(-λ * cf.time) for cf in collect(b)]
        # a hazard tenor inside a payment interval splits it
        split = collect(Projection(risky(b), ctx(Credit.PiecewiseHazard([0.01, 0.02], [0.75, 2.0]), R)))
        @test length(split) == 4 + (6 + 3 + 3 + 6 + 6) * 4
        # portfolios and composites add
        a1, a2 = risky(b), risky(Bond.Fixed(0.03, Periodic(1), 5))
        @test present_value(flat, [a1, a2]) ≈ present_value(flat, a1) + present_value(flat, a2)
        @test present_value(flat, FinanceCore.Composite(a1, a2)) ≈ present_value(flat, a1) + present_value(flat, a2)
        # a missing key is a KeyError
        @test_throws KeyError present_value(Models(curve, Dict("OTHER" => credit(Credit.ConstantHazard(λ), R))), a1)
    end

    @testset "pools that default" begin
        loan = Bond.LevelPayment(0.06, 12, 10)
        speed = Prepayment.CPR(0.08)
        c(h) = Models(curve, Dict("ACME" => credit(h, R), "speed" => speed))
        pool = Bond.Pool(loan, "speed")
        @test present_value(c(Credit.ConstantHazard(0.0)), risky(pool)) ≈ present_value(c(Credit.ConstantHazard(0.0)), pool)
        # by hand: survival and the continuing share weight each payment; recovery uses each unit's balance
        rows = FinanceModels.__asset_rows(pool, c(Credit.ConstantHazard(λ)))
        k, alive, t_prev, total = r + λ, 1.0, 0.0, 0.0
        for row in rows
            total += alive * R * λ * row.opening * (exp(-k * t_prev) - exp(-k * row.time)) / k * exp(λ * t_prev)
            total += alive * exp(-λ * (row.time - t_prev)) * row.payment * exp(-r * row.time)
            alive *= exp(-λ * (row.time - t_prev)) * row.continuing
            t_prev = row.time
        end
        @test present_value(c(Credit.ConstantHazard(λ)), risky(pool)) ≈ total rtol = 1.0e-13
    end

    @testset "derivatives" begin
        b = Bond.Fixed(0.05, Periodic(2), 10)
        v(θ) = present_value(Models(Yield.Constant(Continuous(θ[1])), Dict("ACME" => credit(Credit.ConstantHazard(θ[2]), θ[3]))), risky(b))
        θ = [r, λ, R]
        exact(θ) = analytic(b, θ...)
        @test ForwardDiff.gradient(v, θ) ≈ ForwardDiff.gradient(exact, θ) rtol = 1.0e-12
        @test ForwardDiff.hessian(v, θ) ≈ ForwardDiff.hessian(exact, θ) rtol = 1.0e-11
        # a value with dual numbers keeps the primal value exactly
        d = ForwardDiff.Dual{:t}(λ, 1.0)
        @test ForwardDiff.value(v([r, d, R])) === v(θ)
    end

    @testset "hazard calibration" begin
        truth = Credit.PiecewiseHazard([0.01, 0.015, 0.025], [2.0, 5.0, 10.0])
        bonds = [risky(Bond.Fixed(0.05, Periodic(2), t)) for t in (5.0, 2.0, 10.0)]   # unsorted on purpose
        prices = [present_value(ctx(truth, R), x) for x in bonds]
        template = Credit.Assumptions(Credit.PiecewiseHazard(); recovery = R)
        fitted = fit(template, Quote.(prices, bonds), Fit.Bootstrap(); discount = curve)
        @test fitted.hazard isa Credit.PiecewiseHazard
        @test collect(fitted.hazard.tenors) == [2.0, 5.0, 10.0]
        @test collect(fitted.hazard.rates) ≈ [0.01, 0.015, 0.025] rtol = 1.0e-10
        @test fitted.recovery == R
        for (p, x) in zip(prices, bonds)
            @test present_value(Models(curve, Dict("ACME" => fitted)), x) ≈ p rtol = 1.0e-12
        end
        # amortizing bonds calibrate too
        amort = [risky(Bond.LevelPayment(0.05, 4, t)) for t in (2, 5, 10)]
        aprices = [present_value(ctx(truth, R), x) for x in amort]
        @test collect(fit(template, Quote.(aprices, amort), Fit.Bootstrap(); discount = curve).hazard.rates) ≈ [0.01, 0.015, 0.025] rtol = 1.0e-9

        # first and second derivatives through the calibration
        target = risky(Bond.Fixed(0.06, Periodic(1), 7.0))
        value(p) = present_value(Models(curve, Dict("ACME" => fit(template, Quote.(p, bonds), Fit.Bootstrap(); discount = curve))), target)
        g = ForwardDiff.gradient(value, prices)
        e(i) = 1:3 .== i
        @test g ≈ [(value(prices .+ 1.0e-6e(i)) - value(prices .- 1.0e-6e(i))) / 2.0e-6 for i in 1:3] rtol = 1.0e-7
        H = ForwardDiff.hessian(value, prices)
        Hfd = reduce(hcat, [(ForwardDiff.gradient(value, prices .+ 1.0e-5e(i)) - ForwardDiff.gradient(value, prices .- 1.0e-5e(i))) / 2.0e-5 for i in 1:3])
        @test H ≈ Hfd atol = 1.0e-8
        @test H ≈ H'
        # the same derivatives whatever the chunk size, and from nested calls with their own tags
        chunked(n) = ForwardDiff.gradient(value, prices, ForwardDiff.GradientConfig(value, prices, ForwardDiff.Chunk{n}()))
        @test chunked(1) ≈ chunked(3) rtol = 1.0e-14
        mixed = ForwardDiff.derivative(x -> ForwardDiff.derivative(y -> value([prices[1], x, y]), prices[3]), prices[2])
        @test mixed ≈ H[2, 3] rtol = 1.0e-12
        # derivatives with respect to the recovery assumption and the discount curve
        by_recovery(x) = present_value(
            Models(curve, Dict("ACME" => fit(Credit.Assumptions(Credit.PiecewiseHazard(); recovery = x), Quote.(prices, bonds), Fit.Bootstrap(); discount = curve))), target
        )
        @test ForwardDiff.derivative(by_recovery, R) ≈ (by_recovery(R + 1.0e-6) - by_recovery(R - 1.0e-6)) / 2.0e-6 rtol = 1.0e-6
        by_rate(x) = (c = Yield.Constant(Continuous(x)); present_value(Models(c, Dict("ACME" => fit(template, Quote.(prices, bonds), Fit.Bootstrap(); discount = c))), target))
        @test ForwardDiff.derivative(by_rate, r) ≈ (by_rate(r + 1.0e-6) - by_rate(r - 1.0e-6)) / 2.0e-6 rtol = 1.0e-6
        @test ForwardDiff.derivative(x -> ForwardDiff.derivative(by_rate, x), r) ≈
            (ForwardDiff.derivative(by_rate, r + 1.0e-5) - ForwardDiff.derivative(by_rate, r - 1.0e-5)) / 2.0e-5 rtol = 1.0e-6

        # errors
        dup = [risky(Bond.Fixed(0.05, Periodic(2), 5.0)), risky(Bond.Fixed(0.04, Periodic(2), 5.0))]
        @test_throws "two quotes mature at 5.0" fit(template, Quote.([0.95, 0.9], dup), Fit.Bootstrap(); discount = curve)
        mixed = [risky(Bond.Fixed(0.05, Periodic(2), 2.0)), Credit.Defaultable(Bond.Fixed(0.05, Periodic(2), 5.0), "OTHER")]
        @test_throws "one issuer's bonds" fit(template, Quote.([0.99, 0.97], mixed), Fit.Bootstrap(); discount = curve)
        # a price above the default-free value needs a negative intensity
        rich = [risky(Bond.Fixed(0.05, Periodic(2), 2.0))]
        @test_throws "negative default intensity" fit(template, Quote.([present_value(curve, Bond.Fixed(0.05, Periodic(2), 2.0)) + 0.01], rich), Fit.Bootstrap(); discount = curve)
        # the template has no tenors to evaluate
        @test_throws BoundsError Credit.survival(Credit.PiecewiseHazard(), 1.0)
    end
end
