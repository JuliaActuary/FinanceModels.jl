using ForwardDiff

@testset "amortizing loans" begin
    @testset "Amortizing" begin
        loan = Bond.Amortizing(0.05, [1.0, 2.0, 3.0, 4.0], [25.0, 25.0, 25.0, 25.0])
        @test maturity(loan) == 4.0
        @test [cf.amount for cf in collect(loan)] ≈ [30.0, 28.75, 27.5, 26.25]
        rows = Bond.__schedule(loan)
        @test [r.opening for r in rows] == [100.0, 75.0, 50.0, 25.0]
        @test last(rows).closing === 0.0
        @test all(r -> r.opening - r.principal ≈ r.closing, rows)
        @test sum(r.principal for r in rows) == first(rows).opening

        # a balloon, and explicit accruals (a day count)
        balloon = Bond.Amortizing(0.06, [0.5, 1.0, 1.5], [0.0, 0.0, 1.0]; accruals = [0.5, 0.49, 0.51])
        @test [cf.amount for cf in collect(balloon)] ≈ [0.03, 0.0294, 1.0306]
        # negative amortization adds to the balance
        negam = Bond.Amortizing(0.05, [1.0, 2.0], [-0.1, 1.1])
        @test [r.opening for r in Bond.__schedule(negam)] ≈ [1.0, 1.1]
        @test [cf.amount for cf in collect(negam)] ≈ [0.05 - 0.1, 1.1 * 0.05 + 1.1]

        # valued like its cashflows
        curve = Yield.Constant(Continuous(0.04))
        @test present_value(curve, loan) ≈ sum(cf.amount * exp(-0.04 * cf.time) for cf in collect(loan))
        # a fixed bond's schedule matches its cashflows, including a short first stub
        for b in (Bond.Fixed(0.05, Periodic(2), 3), Bond.Fixed(0.05, Periodic(2), 2.3))
            @test [r.interest + r.principal for r in Bond.__schedule(b)] ≈ [cf.amount for cf in collect(b)]
            @test [r.time for r in Bond.__schedule(b)] == [cf.time for cf in collect(b)]
        end

        @test_throws ArgumentError Bond.Amortizing(0.05, [2.0, 1.0], [0.5, 0.5])
        @test_throws ArgumentError Bond.Amortizing(0.05, [0.0, 1.0], [0.5, 0.5])
        @test_throws ArgumentError Bond.Amortizing(0.05, Float64[], Float64[])
        @test_throws DimensionMismatch Bond.Amortizing(0.05, [1.0, 2.0], [1.0])
    end

    @testset "LevelPayment" begin
        i, n = 0.06 / 12, 360
        m = Bond.LevelPayment(0.06, 12, 30)
        pmts = [cf.amount for cf in collect(m)]
        @test length(pmts) == n
        @test all(p -> p ≈ i / (1 - (1 + i)^-n), pmts)
        @test maximum(pmts) - minimum(pmts) <= 1.0e-14
        @test sum(m.principal) ≈ 1.0
        @test last(Bond.__schedule(m)).closing == 0
        # level payments at the loan's own rate are worth par
        @test present_value(Yield.Constant(Periodic(0.06, 12)), m) ≈ 1.0
        # a Periodic rate of the loan's frequency is its nominal rate
        @test Bond.LevelPayment(Periodic(0.06, 12), 12, 30).principal == m.principal
        @test_throws ArgumentError Bond.LevelPayment(Periodic(0.06, 1), 12, 30)
        @test_throws ArgumentError Bond.LevelPayment(0.06, 12, 30.01)

        # a zero rate repays 1/n each period, and rates near zero keep their precision
        zero_rate = Bond.LevelPayment(0.0, 4, 5)
        @test all(cf -> cf.amount ≈ 1 / 20, collect(zero_rate))
        for r in (1.0e-12, 1.0e-6, 0.0099, 0.0101, 0.05, -0.002)
            @test Bond.__level_payment(r, n) ≈ r / -expm1(-n * log1p(r)) rtol = 1.0e-12
        end
        # derivatives at a zero rate are the limits: (n + 1) / (2n), then (n² - 1) / (6n)
        @test ForwardDiff.derivative(r -> Bond.__level_payment(r, n), 0.0) ≈ (n + 1) / (2n)
        @test ForwardDiff.derivative(r -> ForwardDiff.derivative(x -> Bond.__level_payment(x, n), r), 0.0) ≈ (n^2 - 1) / (6n)
        # and they are continuous across the switch to the closed form
        d(r) = ForwardDiff.derivative(x -> Bond.__level_payment(x, 10), r)
        @test d(prevfloat(0.01)) ≈ d(0.01) rtol = 1.0e-12
        # the schedule differentiates in its rate
        pv_of(r) = present_value(Yield.Constant(0.05), Bond.LevelPayment(r, 12, 10))
        @test ForwardDiff.derivative(pv_of, 0.04) ≈ (pv_of(0.04 + 1.0e-6) - pv_of(0.04 - 1.0e-6)) / 2.0e-6 rtol = 1.0e-7
    end

    @testset "Pool and prepayment" begin
        loan = Bond.LevelPayment(0.06, 12, 10)
        curve = Yield.Constant(0.05)
        pool = Bond.Pool(loan, "prepay")
        ctx(p) = Models(curve, Dict("prepay" => p))
        @test maturity(pool) == maturity(loan)

        # no prepayment is the loan itself
        @test [cf.amount for cf in collect(Projection(pool, ctx(Prepayment.CPR(0.0))))] ≈ [cf.amount for cf in collect(loan)]
        @test present_value(ctx(Prepayment.CPR(0.0)), pool) ≈ present_value(curve, loan)

        # full prepayment repays everything at the first payment
        full = collect(Projection(pool, ctx(Prepayment.CPR(1.0))))
        @test first(full).amount ≈ 1.0 + 0.06 / 12
        @test all(cf -> iszero(cf.amount), full[2:end])

        # the pool's principal, scheduled and prepaid, reconciles to its initial balance
        rows = Bond.__schedule(loan)
        cpr = Prepayment.CPR(0.1)
        smms = FinanceModels.__prepayment_fractions(cpr, rows)
        F = 1.0
        principal_paid = 0.0
        for (r, s) in zip(rows, smms)
            principal_paid += F * (r.principal + s * r.closing)
            F *= 1 - s
        end
        @test principal_paid ≈ 1.0
        # each interval's fraction is 1 - (1 - CPR)^Δt
        @test first(smms) ≈ 1 - 0.9^(1 / 12)

        # an SMM of m and a CPR of 1 - (1 - m)^12 agree on a monthly schedule
        m = 0.007
        @test present_value(ctx(Prepayment.SMM(m)), pool) ≈ present_value(ctx(Prepayment.CPR(1 - (1 - m)^12)), pool) rtol = 1.0e-12
        # SMM applies to monthly payments only
        quarterly = Bond.Pool(Bond.LevelPayment(0.06, 4, 10), "prepay")
        @test_throws ArgumentError present_value(ctx(Prepayment.SMM(m)), quarterly)
        # one rate per payment
        ramp = Prepayment.CPR([min(0.06, 0.002 * k) for k in 1:120])
        # a 6% loan on a 5% curve is a premium loan, which loses value when it prepays at par
        @test present_value(ctx(ramp), pool) < present_value(ctx(Prepayment.CPR(0.0)), pool)
        @test_throws DimensionMismatch present_value(ctx(Prepayment.CPR(fill(0.05, 119))), pool)
        @test_throws DomainError Prepayment.CPR(1.5)
        @test_throws DomainError Prepayment.SMM([0.01, -0.01])

        # a pool scales with its loan, and a fixed bond pools too
        big = Bond.Pool(Bond.Amortizing(0.06, loan.times, 100 .* loan.principal; accruals = loan.accruals), "prepay")
        @test present_value(ctx(cpr), big) ≈ 100 * present_value(ctx(cpr), pool)
        @test present_value(ctx(Prepayment.CPR(0.0)), Bond.Pool(Bond.Fixed(0.05, Periodic(2), 5), "prepay")) ≈
            present_value(curve, Bond.Fixed(0.05, Periodic(2), 5))

        # forward-starting: the assumption counts payments, the same from any start
        @test present_value(ctx(cpr), Forward(2.0, pool)) ≈ present_value(ctx(cpr), pool) * discount(curve, 2.0)

        # derivatives with respect to the prepayment rate
        pv_cpr(c) = present_value(ctx(Prepayment.CPR(c)), pool)
        @test ForwardDiff.derivative(pv_cpr, 0.1) ≈ (pv_cpr(0.1 + 1.0e-6) - pv_cpr(0.1 - 1.0e-6)) / 2.0e-6 rtol = 1.0e-6
    end
end
