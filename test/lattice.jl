using ForwardDiff
using Random

@testset "callable and putable assets on a Hull–White lattice" begin
    curve = ZeroRateCurve([0.03, 0.035, 0.04, 0.042], [1.0, 3.0, 7.0, 12.0])
    hw = ShortRate.HullWhite(0.1, 0.01, curve)
    tree = ShortRate.HullWhiteTree(hw)
    bond = Bond.Fixed(0.05, Periodic(2), 10)
    straight = present_value(curve, bond)
    call = Bond.Callable(bond, 3.0:0.5:9.5, 1.0)

    @testset "HullWhiteTree forwards to its model" begin
        @test tree.max_step == 1 / 12
        @test_throws ArgumentError ShortRate.HullWhiteTree(hw; max_step = 0.0)
        for (f, n) in FinanceModels.Yield.__FORWARDED_CAPABILITIES
            args = n == 0 ? () : n == 1 ? (3.5,) : (1.0, 3.5)
            @test f(tree, args...) == f(hw, args...)
        end
        @test discount(tree, 2.0, 5.0, 0.03) == discount(hw, 2.0, 5.0, 0.03)
        swaption = Option.Swaption(2.0, 7.0, 0.04, 1)
        @test present_value(tree, swaption) == present_value(hw, swaption)
        @test present_value(Models(tree), Option.Cap(0.04, 4, 5.0)) == present_value(hw, Option.Cap(0.04, 4, 5.0))
        paths(m) = simulate(m; n_scenarios = 3, horizon = 2.0, rng = Xoshiro(1))
        @test all(discount(p, 2.0) == discount(q, 2.0) for (p, q) in zip(paths(tree), paths(hw)))
        # straight assets value on its curve without a lattice
        @test present_value(tree, bond) == straight
    end

    @testset "exercise that never pays reproduces the asset" begin
        # the lattice reprices the curve at every level, so a straight asset is exact
        @test present_value(tree, Bond.Callable(bond, 3.0:0.5:9.5, 100.0)) ≈ straight rtol = 1.0e-14
        @test present_value(tree, Bond.Putable(bond, 3.0:0.5:9.5, 0.0)) ≈ straight rtol = 1.0e-14
        loan = Bond.LevelPayment(0.05, 12, 10)
        @test present_value(tree, Bond.Callable(loan, 1.0:1.0:9.0, 100.0)) ≈ present_value(curve, loan) rtol = 1.0e-14
    end

    @testset "calls lower and puts raise the value" begin
        v = present_value(tree, call)
        @test v < straight
        @test present_value(tree, Bond.Putable(bond, 3.0:0.5:9.5, 1.0)) > straight
        # a lower call price, or more call dates, can only lower it
        @test present_value(tree, Bond.Callable(bond, 3.0:0.5:9.5, 0.99)) < v
        @test present_value(tree, Bond.Callable(bond, 1.0:0.5:9.5, 1.0)) < v
        # negative rates
        neg = ShortRate.HullWhiteTree(ShortRate.HullWhite(0.05, 0.005, Yield.Constant(Continuous(-0.01))))
        vn = present_value(neg, call)
        @test isfinite(vn) && vn < present_value(Yield.Constant(Continuous(-0.01)), bond)
    end

    @testset "one exercise date against Jamshidian's formula" begin
        # callable at par once, after the coupon at 3: the bond less a call on its remainder struck at 1,
        # which is a receiver swaption; a putable bond is the bond plus a payer swaption
        receiver = present_value(hw, Option.Swaption(3.0, 10.0, 0.05, 2; payer = false))
        payer = present_value(hw, Option.Swaption(3.0, 10.0, 0.05, 2; payer = true))
        fine = ShortRate.HullWhiteTree(hw; max_step = 1 / 96)
        @test present_value(fine, Bond.Callable(bond, [3.0], 1.0)) ≈ straight - receiver atol = 5.0e-5
        @test present_value(fine, Bond.Putable(bond, [3.0], 1.0)) ≈ straight + payer atol = 5.0e-5
        # the values converge as the steps shrink
        err(ms) = abs(present_value(ShortRate.HullWhiteTree(hw; max_step = ms), Bond.Callable(bond, [3.0], 1.0)) - (straight - receiver))
        @test err(1 / 96) < err(1 / 12) / 4
    end

    @testset "QuantLib's tree" begin
        # test/reference/quantlib_callable.py: QuantLib 1.44's TreeCallableFixedRateBondEngine with 1000
        # steps over 10 years, on the same conventions. Equal step counts give the same lattice (the
        # script also shows 4000 steps agreeing with max_step = 1/400 to 1e-13).
        cases = [
            (0.04, 0.1, 0.01, 0.05, 3.0, Bond.Callable, 1.01432023176),
            (0.04, 0.1, 0.01, 0.05, 3.0, Bond.Putable, 1.092517089406),
            (0.02, 0.03, 0.015, 0.03, 2.0, Bond.Callable, 0.976477712818),
        ]
        for (r, a, σ, c, first, kind, ql1000) in cases
            asset = kind(Bond.Fixed(c, Periodic(2), 10), first:0.5:9.5, 1.0)
            m = ShortRate.HullWhite(a, σ, Yield.Constant(Continuous(r)))
            @test present_value(ShortRate.HullWhiteTree(m; max_step = 1 / 100), asset) ≈ ql1000 atol = 1.0e-11
        end
    end

    @testset "zero volatility and no mean reversion" begin
        flat_hw = ShortRate.HullWhite(0.1, 0.0, curve)
        v0 = present_value(ShortRate.HullWhiteTree(flat_hw), call)
        # by hand: roll back on the curve's forward discount factors, redeeming where cheaper
        rows = Bond.__schedule(bond)
        V = 0.0
        for i in length(rows):-1:1
            r = rows[i]
            cont = V
            if r.time in 3.0:0.5:9.5
                cont = min(cont, 1.0)
            end
            V = discount(curve, i == 1 ? 0.0 : rows[i - 1].time, r.time) * (r.interest + r.principal + cont)
        end
        @test v0 ≈ V rtol = 1.0e-13
        # continuous as σ → 0
        @test present_value(ShortRate.HullWhiteTree(ShortRate.HullWhite(0.1, 1.0e-9, curve)), call) ≈ v0 atol = 1.0e-7
        # without mean reversion (Ho–Lee)
        ho_lee = present_value(ShortRate.HullWhiteTree(ShortRate.HullWhite(0.0, 0.01, curve)), call)
        @test ho_lee ≈ present_value(ShortRate.HullWhiteTree(ShortRate.HullWhite(1.0e-8, 0.01, curve)), call) rtol = 1.0e-7
    end

    @testset "defaultable and pooled assets" begin
        acme = Credit.Assumptions(Credit.PiecewiseHazard([0.01, 0.02], [4.25, 8.0]); recovery = 0.4)
        risky_straight = present_value(Models(curve, Dict("ACME" => acme)), Credit.Defaultable(bond, "ACME"))
        ctx(ms) = Models(ShortRate.HullWhiteTree(hw; max_step = ms), Dict("ACME" => acme))
        never = Bond.Callable(bond, 3.0:0.5:9.5, 100.0)
        # within a step the lattice's recovery is exact at the node's rate: second order in the step
        @test present_value(ctx(1 / 12), Credit.Defaultable(never, "ACME")) ≈ risky_straight atol = 1.0e-7
        @test present_value(ctx(1 / 48), Credit.Defaultable(never, "ACME")) ≈ risky_straight atol = 1.0e-8
        # the wrappers compose in either order
        @test present_value(ctx(1 / 12), Credit.Defaultable(call, "ACME")) ==
            present_value(ctx(1 / 12), Bond.Callable(Credit.Defaultable(bond, "ACME"), 3.0:0.5:9.5, 1.0))
        @test present_value(ctx(1 / 12), Credit.Defaultable(call, "ACME")) < present_value(tree, call)
        # a pool redeemed after its prepayment
        loan = Bond.LevelPayment(0.06, 12, 10)
        store = Dict("speed" => Prepayment.CPR(0.1))
        pool = Bond.Pool(loan, "speed")
        @test present_value(Models(tree, store), Bond.Callable(pool, 2.0:1.0:9.0, 100.0)) ≈ present_value(Models(curve, store), pool) rtol = 1.0e-13
        @test present_value(Models(tree, store), Bond.Callable(pool, 2.0:1.0:9.0, 1.0)) < present_value(Models(curve, store), pool)
    end

    @testset "derivatives" begin
        fine(m) = ShortRate.HullWhiteTree(m; max_step = 1 / 48)
        f(s) = present_value(fine(ShortRate.HullWhite(0.1, 0.01, curve + Continuous(s))), call)
        @test ForwardDiff.derivative(f, 0.0) ≈ (f(1.0e-6) - f(-1.0e-6)) / 2.0e-6 rtol = 1.0e-5
        g(σ) = present_value(fine(ShortRate.HullWhite(0.1, σ, curve)), call)
        @test ForwardDiff.derivative(g, 0.01) ≈ (g(0.01 + 1.0e-7) - g(0.01 - 1.0e-7)) / 2.0e-7 rtol = 1.0e-4
        acme(λ) = Models(fine(hw), Dict("ACME" => Credit.Assumptions(Credit.ConstantHazard(λ); recovery = 0.4)))
        h(λ) = present_value(acme(λ), Credit.Defaultable(call, "ACME"))
        @test ForwardDiff.derivative(h, 0.02) ≈ (h(0.02 + 1.0e-7) - h(0.02 - 1.0e-7)) / 2.0e-7 rtol = 1.0e-5
        # derivatives with respect to the redemption price and to the schedule's principal
        k(price) = present_value(fine(hw), Bond.Callable(bond, 3.0:0.5:9.5, price))
        @test ForwardDiff.derivative(k, 1.0) ≈ (k(1.0 + 1.0e-7) - k(1.0 - 1.0e-7)) / 2.0e-7 rtol = 1.0e-5
        amort(x) = present_value(fine(hw), Bond.Callable(Bond.Amortizing(0.05, [1.0, 2.0, 3.0], [0.0, x, 1.0 - x]), [1.0], 1.0))
        @test ForwardDiff.derivative(amort, 0.3) ≈ (amort(0.3 + 1.0e-6) - amort(0.3 - 1.0e-6)) / 2.0e-6 rtol = 1.0e-6
        # a dual value keeps the primal value
        @test ForwardDiff.value(ForwardDiff.Dual{:t}(1.0, 0.0) * f(ForwardDiff.Dual{:t}(0.0, 1.0))) === f(0.0)
        # second derivatives through exercise decisions leave out the moving boundary: they throw
        @test_throws "second derivatives through an exercise decision" ForwardDiff.derivative(s -> ForwardDiff.derivative(f, s), 0.0)
        # the lattice has no derivative in σ at σ = 0
        @test_throws "not differentiable in σ at σ = 0" ForwardDiff.derivative(g, 0.0)
    end

    @testset "errors" begin
        @test_throws "not a payment date" present_value(tree, Bond.Callable(bond, [3.25], 1.0))
        @test_throws "exercise at maturity" present_value(tree, Bond.Callable(bond, [10.0], 1.0))
        @test_throws "HullWhiteTree" present_value(curve, call)
        @test_throws "HullWhiteTree" present_value(hw, call)
        @test_throws "two exercise schedules" present_value(tree, Bond.Putable(call, [4.0], 1.0))
        @test_throws "Credit.Defaultable twice" present_value(Models(tree), Bond.Callable(Credit.Defaultable(Credit.Defaultable(bond, "A"), "B"), [3.0], 1.0))
        @test_throws ArgumentError Bond.Callable(bond, [4.0, 3.0], 1.0)
        @test_throws DimensionMismatch Bond.Callable(bond, [3.0, 4.0], [1.0])
        @test maturity(call) == 10
    end

    @testset "exercise policy" begin
        ctx = Models(tree)
        policy = exercise_policy(ctx, call)
        @test policy.times == collect(3.0:0.5:9.5)
        @test policy.call
        # the issuer calls when rates are low
        @test policy(3.0, 0.0)
        @test !policy(3.0, 0.12)
        @test policy(5.0, Continuous(0.0))
        @test_throws "not an exercise time" policy(3.25, 0.03)
        @test_throws "beyond the exercise policy's lattice" policy(3.0, 5.0)
        # the policy's wider lattice values the asset the same: added nodes carry no weight
        inputs = FinanceModels.__exercise_inputs(ctx, call, tree.max_step)
        narrow = FinanceModels.__build_lattice(hw, inputs.times, 0.0)
        wide = FinanceModels.__build_lattice(hw, inputs.times, 0.3)
        v(lat) = first(FinanceModels.__rollback(lat, inputs.rows, inputs.events, inputs.exercise, inputs.opening, inputs.credit))
        @test v(wide) ≈ v(narrow) rtol = 1.0e-14
        @test v(narrow) == present_value(ctx, call)
        # its nodes span ±8 standard deviations of the state at maturity from the first exercise time
        sd = 0.01 * sqrt((1 - exp(-2 * 0.1 * 10)) / (2 * 0.1))
        @test policy.lo[1] * policy.dx[1] <= -8sd
        # a put is exercised when rates are high
        put_policy = exercise_policy(ctx, Bond.Putable(bond, 3.0:1.0:9.0, 1.0))
        @test !put_policy.call
        @test put_policy(3.0, 0.12)
        @test !put_policy(3.0, 0.0)
        @test_throws "HullWhiteTree" exercise_policy(Models(curve), call)
    end
end
