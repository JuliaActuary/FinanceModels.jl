using Random
using DataInterpolations

@testset "asset-liability projections" begin
    curve = Yield.Constant(Continuous(0.04))
    acme(λ; R = 0.5) = Credit.Assumptions(Credit.ConstantHazard(λ); recovery = R)
    store = Dict("speed" => Prepayment.CPR(0.08), "ACME" => acme(0.01))
    ctx = Models(curve, store)
    reconciles(rows) = all(r -> isapprox(r.opening - r.principal - r.prepayment - r.defaulted - r.redeemed, r.closing; atol = 1.0e-14), rows) &&
        all(i -> rows[i].closing == rows[i + 1].opening, 1:(length(rows) - 1))

    @testset "expected runoff" begin
        loan = Bond.LevelPayment(0.06, 12, 10)
        pool = Credit.Defaultable(Bond.Pool(loan, "speed"), "ACME")
        rows = collect(Projection(pool, ctx, AssetProjection()))
        @test length(rows) == 120
        @test reconciles(rows)
        @test first(rows).opening == 1.0
        @test last(rows).closing == 0.0
        @test all(r -> r.recovery ≈ 0.5 * r.defaulted, rows)
        # the payments are the expected cashflows' payments...
        cfs = collect(Projection(pool, ctx))
        pays = [cf.amount for cf in cfs if any(r -> r.time == cf.time, rows)]
        @test [r.interest + r.principal + r.prepayment for r in rows] ≈ pays rtol = 1.0e-13
        # ...and each interval's recovery is the sum of its quadrature recoveries, undiscounted
        times = [0.0; [r.time for r in rows]]
        rec = [sum(cf.amount for cf in cfs if times[i] < cf.time < times[i + 1]) for i in eachindex(rows)]
        @test [r.recovery for r in rows] ≈ rec rtol = 1.0e-12
        # without default risk or prepayment, the rows are the schedule
        plain = collect(Projection(loan, ctx, AssetProjection()))
        sched = Bond.__schedule(loan)
        @test [r.interest for r in plain] == [s.interest for s in sched]
        @test [r.principal for r in plain] == [s.principal for s in sched]
        @test all(r -> iszero(r.defaulted) && iszero(r.prepayment), plain)
        @test reconciles(plain)
    end

    @testset "sampled defaults" begin
        b1 = Credit.Defaultable(Bond.Fixed(0.05, Periodic(2), 10), "ACME")
        b2 = Credit.Defaultable(Bond.Fixed(0.03, Periodic(1), 15), "ACME")
        # an explicit default time: payments before it, then the recovery of the outstanding principal
        rows = collect(Projection(b1, ctx, AssetProjection(; default_times = Dict("ACME" => 4.2))))
        @test [r.time for r in rows] == [0.5:0.5:4.0; 4.2]
        @test all(r -> iszero(r.defaulted), rows[1:(end - 1)])
        @test last(rows).defaulted == 1.0
        @test last(rows).recovery == 0.5
        @test last(rows).closing == 0.0
        @test reconciles(rows)
        # no default
        @test length(collect(Projection(b1, ctx, AssetProjection(; default_times = Dict("ACME" => Inf))))) == 20
        # a default exactly at a payment date comes before that payment
        at = collect(Projection(b1, ctx, AssetProjection(; default_times = Dict("ACME" => 4.0))))
        @test last(at).time == 4.0 && last(at).defaulted == 1.0 && iszero(last(at).interest)

        # an issuer's bonds default together, at a time drawn once
        risky = Models(curve, Dict("ACME" => acme(0.3)))
        kind = AssetProjection(; rng = Xoshiro(3))
        both = collect(Projection([b1, b2], risky, kind))
        τ = kind.default_times["ACME"]
        @test τ < 10
        @test [r.time for r in both if r.defaulted > 0] == [τ, τ]
        # reproducible from the seed, and reused by later projections with the same kind
        @test Credit.__default_time(Credit.ConstantHazard(0.3), randexp(Xoshiro(3))) == τ
        @test last(collect(Projection(b2, risky, kind))).time == τ
        @test_throws KeyError collect(Projection(b1, risky, AssetProjection(; default_times = Dict("OTHER" => 1.0))))

        # sampled defaults average to the expected runoff
        expected = collect(Projection(b1, risky, AssetProjection()))
        n = 20_000
        rng = Xoshiro(11)
        totals = map(1:n) do _
            rs = collect(Projection(b1, risky, AssetProjection(; rng)))
            (sum(r.interest for r in rs), sum(r.recovery for r in rs))
        end
        for (k, field) in ((1, :interest), (2, :recovery))
            xs = getindex.(totals, k)
            m, se = sum(xs) / n, sqrt(sum(abs2, xs .- sum(xs) / n) / (n - 1) / n)
            @test abs(m - sum(getfield(r, field) for r in expected)) < 4se
        end

        # pools use expected defaults in a sampled projection
        pool = Credit.Defaultable(Bond.Pool(Bond.LevelPayment(0.06, 12, 10), "speed"), "ACME")
        @test collect(Projection(pool, ctx, AssetProjection(; rng = Xoshiro(1)))) == collect(Projection(pool, ctx, AssetProjection()))

        # the inverse cumulative hazard
        h = Credit.PiecewiseHazard([0.1, 0.0, 0.2], [1.0, 2.0, 3.0])
        for E in (0.05, 0.1, 0.15, 0.3, 1.0)
            @test Credit.cumhazard(h, Credit.__default_time(h, E)) ≈ E
        end
        @test Credit.__default_time(Credit.PiecewiseHazard([0.1, 0.0], [1.0, 2.0]), 0.5) == Inf
        @test Credit.__default_time(Credit.ConstantHazard(0.0), 0.5) == Inf
    end

    @testset "exercise along rate paths" begin
        hw = ShortRate.HullWhite(0.1, 0.01, ZeroRateCurve([0.03, 0.035, 0.04, 0.042], [1.0, 3.0, 7.0, 12.0]))
        call = Bond.Callable(Bond.Fixed(0.05, Periodic(2), 10), 3.0:0.5:9.5, 1.0)
        pricing = Models(ShortRate.HullWhiteTree(hw))
        price = present_value(pricing, call)
        paths = simulate(hw; n_scenarios = 4000, timestep = 1 / 12, horizon = 10.0, rng = Xoshiro(7))
        kind = AssetProjection(; pricing)
        projected = [collect(Projection(call, Models(path), kind)) for path in paths]
        @test all(reconciles, projected)
        # the policy was computed once, for this asset
        @test length(kind.policies) == 1
        # paths of the pricing model value the bond at its lattice price
        vals = [sum((r.interest + r.principal + r.recovery + r.redemption) * discount(path, r.time) for r in rows) for (rows, path) in zip(projected, paths)]
        m = sum(vals) / length(vals)
        se = sqrt(sum(abs2, vals .- m) / (length(vals) - 1) / length(vals))
        @test abs(m - price) < 4se
        # a called bond redeems at par after its coupon and pays nothing more
        called = first(filter(rows -> last(rows).redeemed > 0, projected))
        @test last(called).redemption == last(called).redeemed
        @test last(called).closing == 0 && last(called).time in 3.0:0.5:9.5

        # a decision reads only the rate observed then: changing the path after a date changes no earlier row
        path = paths[1]
        ts = FinanceModels.simulation_times(path)
        rates = [rate(short_rate(path, t)) for t in ts]
        k = findfirst(>=(6.0), ts)
        altered = copy(rates)
        altered[(k + 1):end] .= -0.05
        integral = [0.0; cumsum([(altered[j] + altered[j + 1]) / 2 * (ts[j + 1] - ts[j]) for j in 1:(length(ts) - 1)])]
        other = RatePath(LinearInterpolation(integral, ts; extrapolation = ExtrapolationType.None), altered)
        before(rows) = [r for r in rows if r.time <= 6.0]
        @test before(collect(Projection(call, Models(other), kind))) == before(projected[1])

        # an exercisable asset needs the pricing context
        @test_throws "pass `pricing`" collect(Projection(call, Models(path), AssetProjection()))
        # a defaultable callable bond along a path, with expected defaults
        rstore = Dict("ACME" => acme(0.02; R = 0.4))
        rkind = AssetProjection(; pricing = Models(ShortRate.HullWhiteTree(hw), rstore))
        rrows = collect(Projection(Credit.Defaultable(call, "ACME"), Models(path, rstore), rkind))
        @test reconciles(rrows)
        @test all(r -> r.recovery ≈ 0.4 * r.defaulted, rrows)
    end
end
