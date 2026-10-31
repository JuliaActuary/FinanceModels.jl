using TOML

@testset "Official EIOPA May 2026 workbook" begin
    fixture = TOML.parsefile(joinpath(@__DIR__, "data", "eiopa_2026.toml"))
    p, v, expected = fixture["inputs"], fixture["va"], fixture["expected"]
    i = EIOPA.RFRInputs(
        p["rates"], p["tenors"]; fsp = p["fsp"], llfr_weights = p["llfr_weights"],
        ufr = Periodic(p["ufr"], 1), α = p["alpha"], cra_bp = p["cra_bp"], frequency = p["frequency"]
    )
    va = EIOPA.VAInputs(;
        rcs_currency = v["rcs_currency"], usa = v["usa_selected"] ? v["usa"] : nothing,
        rcs_country = v["rcs_country"], debt_to_assets = v["debt_to_assets"], use_macro_va = v["macro_selected"],
        valuations = (;
            assets = v["assets"], assets_notional = v["assets_notional"],
            liabilities = v["liabilities"], liabilities_notional = v["liabilities_notional"],
        )
    )
    shocks = EIOPA.ShockParameters()
    result = EIOPA.publication(i; va_inputs = va, shock_parameters = shocks)
    precise = EIOPA.calculate(i; va_inputs = va, shock_parameters = shocks)
    @test rate.(zero.(Ref(precise.bootstrapped), 1:50)) ≈ expected["bootstrap"] atol = 1.0e-13
    for name in (:llfr, :up_llfr, :down_llfr, :cssr_llfr, :va_llfr, :up_va_llfr, :down_va_llfr, :notional_va, :currency_va, :country_weight)
        @test getproperty(result.diagnostics, name) ≈ expected[string(name)] atol = 1.0e-14 rtol = 0
    end
    for name in (:assets_pvbp, :liabilities_pvbp)
        @test getproperty(result.diagnostics, name) ≈ expected[string(name)] rtol = 2.0e-14
    end
    @test result.diagnostics.cssr == expected["cssr_value"]
    @test result.diagnostics.total_va ≈ expected["total_va"] atol = 1.0e-15
    for name in keys(result.rates)
        @test getproperty(result.rates, name) ≈ expected[string(name)] atol = 1.0e-14 rtol = 0
    end
    for name in (:up_relative, :up_absolute, :down_relative, :down_absolute, :floor)
        @test getproperty(shocks, name) ≈ fixture["shocks"][string(name)] atol = 1.0e-15
    end
    @test result.rates.basic == EIOPA._round_rate.(EIOPA._annual.(Ref(precise.curves.basic), 1:150))
    @test rate(zero(precise.curves.basic, Inf)) ≈ log1p(p["ufr"])
    @test EIOPA.publication(i).rates.va === nothing
    @test EIOPA.calculate(i; va_inputs = EIOPA.VAInputs(rcs_currency = 0.01)).curves.cssr !== nothing
    @test EIOPA.calculate(i; va_inputs = EIOPA.VAInputs(rcs_currency = 0.01)).curves.va === nothing

    # Article 51a shocks assets in parallel even beyond FSP. Using the fading
    # liability CSSR curve for both sides would incorrectly return CSSR=0.5.
    n = EIOPA.calculate(i; va_inputs = EIOPA.VAInputs(rcs_currency = 0.01))
    asset_curve = EIOPA.asset_spread_shock(n.curves.basic, Continuous(log1p(n.diagnostics.notional_va)))
    for t in (5, 20, 60, 150)
        @test discount(asset_curve, t) ≈ discount(n.curves.basic, t) * (1 + n.diagnostics.notional_va)^(-t)
    end
    assets = 100discount(n.curves.basic, 60)
    liabilities = 200discount(n.curves.basic, 60)
    assets_notional = 100discount(asset_curve, 60)
    liabilities_notional = 200discount(n.curves.cssr, 60)
    correct = EIOPA.VAInputs(rcs_currency = 0.01; valuations = (; assets, assets_notional, liabilities, liabilities_notional))
    @test EIOPA.calculate(i; va_inputs = correct).diagnostics.cssr == 1
    @test EIOPA.cssr(assets - 100discount(n.curves.cssr, 60), liabilities - liabilities_notional) ≈ 0.5
end

@testset "Asset spread shock conventions" begin
    base = Yield.Constant(Periodic(0.04, 1))
    va = 0.0085
    annual = @inferred EIOPA.asset_spread_shock(base, Periodic(va, 1))
    continuous = @inferred EIOPA.asset_spread_shock(base, Continuous(log1p(va)))
    semiannual = @inferred EIOPA.asset_spread_shock(base, Periodic(va, 2))
    @test annual isa Yield.TenorShift
    # Price a 10-year 4% coupon bond with the bumped yields by hand. A nominal
    # semiannual bump uses the supplied number, rather than an annual conversion.
    annual_yield = 0.04 + va
    multiplicative_yield = 1.04 * 1.0085 - 1
    semiannual_yield = 2 * (sqrt(1.04) - 1) + va
    hand_annual(y) = sum(4 / (1 + y)^t for t in 1:10) + 100 / (1 + y)^10
    hand_semiannual(y) = sum(4 / (1 + y / 2)^(2t) for t in 1:10) + 100 / (1 + y / 2)^20
    price(c) = sum(4discount(c, t) for t in 1:10) + 100discount(c, 10)
    @test price(annual) ≈ hand_annual(annual_yield) rtol = 2.0e-14
    @test price(continuous) ≈ hand_annual(multiplicative_yield) rtol = 2.0e-14
    @test price(semiannual) ≈ hand_semiannual(semiannual_yield) rtol = 2.0e-14
    @test price(annual) > price(continuous)
    old = base + Yield.Constant(Periodic(va, 1))
    for t in (0.0, 1.0, 10.0, 60.0, Inf)
        @test discount(continuous, t) == discount(old, t)
    end
    @test discount(annual, Inf) == 0.0
    @test discount(semiannual, Inf) == 0.0
    @test discount(EIOPA.asset_spread_shock(base, Periodic(-0.05, 1)), Inf) == Inf
    # A zero asymptotic yield can conceal a nonzero log-discount intercept.
    tail = Yield.ExponentialForward(base; start = 20, initial_forward = Continuous(0.03), ufr = Continuous(0.0), α = 0.11)
    @test isfinite(discount(tail, Inf))
    @test_throws ArgumentError discount(EIOPA.asset_spread_shock(tail, Periodic(0.0, 1)), Inf)
    @test discount(EIOPA.asset_spread_shock(tail, Continuous(0.0)), Inf) == discount(tail, Inf)
    @test_throws ArgumentError discount(EIOPA.asset_spread_shock(Yield.Constant(Periodic(-0.04, 1)), Periodic(-0.99, 1)), 10)
    for shift in (Continuous(Inf), Continuous(NaN), Periodic(-1.0, 1))
        @test_throws ArgumentError EIOPA.asset_spread_shock(base, shift)
    end
    value(v) = discount(EIOPA.asset_spread_shock(base, Periodic(v, 2)), 10)
    @test ForwardDiff.derivative(value, va) ≈ (value(va + 1.0e-6) - value(va - 1.0e-6)) / 2.0e-6 rtol = 1.0e-7
end

@testset "Publication decimal rounding" begin
    for bp in 1:500
        expected = Float64(round(BigInt, BigInt(85) * bp // 100, RoundNearestTiesAway)) / 10000
        @test EIOPA._round_va(0.85 * (bp / 10000)) == expected
        @test EIOPA._round_va(-0.85 * (bp / 10000)) == -expected
    end
    # Native Excel worksheet ROUND results, also checked in a scratch workbook.
    for (x, digits, expected) in (
            (1.005, 2, 1.01), (2.675, 2, 2.68), (-1.475, 2, -1.48),
            (2.15, 1, 2.2), (0.00005, 4, 0.0001), (-0.00005, 4, -0.0001),
        )
        @test EIOPA._excel_round(x, digits) == expected
    end
    @test EIOPA._excel_round(1.0e16, -2) == 1.0e16
    @test iszero(EIOPA._round_va(nextfloat(0.0)))
    @test signbit(EIOPA._round_va(-nextfloat(0.0)))
    i = EIOPA.RFRInputs([0.02], [20]; fsp = 20, llfr_weights = [1.0], ufr = Periodic(0.033, 1), α = 0.11)
    r = EIOPA.publication(i; va_inputs = EIOPA.VAInputs(rcs_currency = 0.007, cssr = 1.0))
    @test r.diagnostics.notional_va == 0.006
    @test r.diagnostics.total_va == 0.006
    @test r.rates.va[1] == EIOPA._round_rate((1 + r.rates.basic[1]) * 1.006 - 1)
end

@testset "Publication rounding fast path and fallback" begin
    for digits in (0, 2, 4, 5, 22, 23), n in (-100, -3, -1, 0, 1, 3, 100)
        tie = Float64((big(n) + 1 // 2) / big(10)^digits)
        for x in (prevfloat(tie), tie, nextfloat(tie))
            @test isequal(EIOPA._excel_round(x, digits), EIOPA._excel_round_exact(x, digits))
        end
    end
    # Either side of the scaled-magnitude guard, both signs and adjacent floats.
    for point in (999.99999, 1000.00001), sign in (-1, 1)
        for x in (prevfloat(sign * point), sign * point, nextfloat(sign * point))
            @test isequal(EIOPA._excel_round(x, 5), EIOPA._excel_round_exact(x, 5))
        end
    end
    # Either side of the 1e-6 half-integer margin, including negative values.
    for distance in (0.9e-6, 1.1e-6), side in (-1, 1), sign in (-1, 1)
        x = sign * (100.5 + side * distance) / 1.0e5
        @test isequal(EIOPA._excel_round(x, 5), EIOPA._excel_round_exact(x, 5))
    end
    for x in (0.0, -0.0, Inf, -Inf, NaN, nextfloat(0.0), -nextfloat(0.0), 1.0e15), digits in (0, 4, 5)
        @test isequal(EIOPA._excel_round(x, digits), EIOPA._excel_round_exact(x, digits))
    end
    @test EIOPA._excel_round(1.0e-23, 23) == 1.0e-23
    @test EIOPA._excel_round(1.0e16, -2) == EIOPA._excel_round_exact(1.0e16, -2)
    for x in (Float32(0.00595), big"0.00595")
        @test isequal(EIOPA._excel_round(x, 4), EIOPA._excel_round_exact(x, 4))
    end
end

@testset "Downward VA ordering and exact floor tails" begin
    make(r, ufr = 0.033; α = 0.11) = EIOPA.RFRInputs(
        [r, r], [1, 20]; instrument = :government,
        fsp = 20, llfr_weights = [0.0, 1.0], ufr = Periodic(ufr, 1), α
    )
    p = EIOPA.ShockParameters()
    v = EIOPA.VAInputs(rcs_currency = 0.01, cssr = 1.0)
    for r in (-0.005, -0.04)
        i = make(r)
        raw_down = r * (1 - p.down_relative[1]) - p.down_absolute[1]
        raw_with_va = (1 + raw_down) * 1.0085 - 1
        published = EIOPA.publication(i; va_inputs = v, shock_parameters = p)
        @test published.rates.down[1] == p.floor[1]
        @test published.rates.down_va[1] == max(EIOPA._round_rate(raw_with_va), p.floor[1])
        @test published.rates.down_va[1] != EIOPA._round_rate((1 + p.floor[1]) * 1.0085 - 1)
        precise = EIOPA.calculate(i; va_inputs = v, shock_parameters = p)
        @test EIOPA._annual(precise.curves.down_va, 1) ≈ max(raw_with_va, p.floor[1]) atol = 1.0e-15
    end
    i = make(0.03, 0.0015)
    raw = EIOPA._calculate(i, v, p, identity, identity)
    r = EIOPA.calculate(i; va_inputs = v, shock_parameters = p)
    for name in (:down, :down_va)
        c, unfloored = getproperty(r.curves, name), getproperty(raw.curves, name)
        @test discount(c, Inf) ≈ discount(unfloored, Inf)
        @test discount(c, 0, Inf) ≈ discount(c, Inf)
        @test discount(c, Inf, Inf) == 1
        @test discount(c, 1.0e6) ≈ discount(c, Inf) rtol = 1.0e-5
    end
    negative = EIOPA.calculate(make(0.03, -0.01); shock_parameters = p).curves.down
    @test rate(zero(negative, Inf)) ≈ log1p(last(p.floor))
    @test discount(negative + Yield.Constant(Continuous(-log1p(last(p.floor)))), Inf) ≈ 1
    @test discount(EIOPA.calculate(make(0.03); shock_parameters = p).curves.down, Inf) == 0
    at_infinity(a) = discount(EIOPA.calculate(make(0.03, 0.0015; α = a); shock_parameters = p).curves.down, Inf)
    @test ForwardDiff.derivative(at_infinity, 0.11) ≈ (at_infinity(0.110001) - at_infinity(0.109999)) / 0.000002 rtol = 1.0e-7
    @test_throws ArgumentError ForwardDiff.derivative(u -> discount(EIOPA.calculate(make(0.03, u); shock_parameters = p).curves.down, Inf), 0.0015)
end

@testset "EIOPA conventions and edge cases" begin
    input(r = [0.02, 0.03]; kwargs...) = EIOPA.RFRInputs(
        r, [2, 5]; fsp = 5,
        llfr_weights = [0.0, 1.0], ufr = Periodic(0.033, 1), α = 0.11, kwargs...
    )
    i = input(; instrument = :government, cra_bp = 10)
    r = EIOPA.calculate(i)
    @test rate(zero(r.bootstrapped, 0)) ≈ log1p(0.019)
    @test discount(r.bootstrapped, 3) ≈ exp(-(2log1p(0.019) + (5log1p(0.029) - 2log1p(0.019)) / 3))
    @test r.diagnostics.llfr ≈ (5log1p(0.029) - 2log1p(0.019)) / 3
    @test EIOPA.cssr(1.0, 0.0) == 1
    @test EIOPA.cssr(-2.0, -1.0) == 1
    @test EIOPA.cssr(-1.0, 2.0) == 0
    @test EIOPA.cssr(3.0, 2.0) == 1
    for sign in (-1, 0, 1)
        v = EIOPA.VAInputs(rcs_currency = sign * 0.01, cssr = 0.5)
        precise = EIOPA.calculate(i; va_inputs = v)
        @test precise.diagnostics.total_va ≈ sign * 0.00425
        @test rate(Yield.instantaneous_forward(precise.curves.va, 5)) ≈ r.diagnostics.llfr + log1p(sign * 0.00425)
    end
    for (debt, weight) in ((0.2, 0.0), (0.35, 1 / 3), (0.5, 1.0))
        v = EIOPA.VAInputs(rcs_currency = 0.01, rcs_country = 0.02, debt_to_assets = debt, use_macro_va = true, cssr = 0.5)
        @test EIOPA.calculate(i; va_inputs = v).diagnostics.macro_va ≈ 0.85 * 0.5 * 0.007 * weight atol = 1.0e-15
    end
    @test_throws ArgumentError EIOPA.VAInputs(rcs_currency = 0.01, use_macro_va = true, usa = 1.0)
    @test_throws ArgumentError EIOPA.VAInputs(rcs_currency = 0.01, usa = 1.051)
    @test EIOPA.calculate(i; va_inputs = EIOPA.VAInputs(rcs_currency = 0.01, usa = 1.05, cssr = 1.0)).diagnostics.total_va ≈ 0.008925
    @test_throws ArgumentError EIOPA.VAInputs(rcs_currency = 0.01, rcs_country = Inf)
    @test_throws ArgumentError EIOPA.VAInputs(rcs_currency = 0.01, debt_to_assets = 1.01)
    @test_throws MethodError EIOPA.asset_spread_shock(r.curves.basic, 0.01)
    # The DLT assessment can locate FSP before the first positive LLFR weight.
    zero_weight = EIOPA.RFRInputs(
        [0.02, 0.03], [5, 10]; fsp = 5,
        llfr_weights = [0.0, 1.0], ufr = Periodic(0.033, 1), α = 0.11
    )
    without_va = EIOPA.calculate(zero_weight)
    with_va = EIOPA.calculate(zero_weight; va_inputs = EIOPA.VAInputs(rcs_currency = 0.01, cssr = 1.0))
    @test with_va.diagnostics.va_llfr == without_va.diagnostics.llfr
    @test discount(with_va.curves.va, 5) ≈ discount(without_va.curves.basic, 5) * 1.0085^(-5)
    @test rate(Yield.instantaneous_forward(with_va.curves.va, 60)) == rate(Yield.instantaneous_forward(without_va.curves.basic, 60))
    @test_throws ArgumentError EIOPA.calculate(i; va_inputs = EIOPA.VAInputs(rcs_currency = 0.0, valuations = (assets = 1.0, assets_notional = 1.0, liabilities = 2.0, liabilities_notional = 2.0)))
    @test EIOPA._round_va(0.00005) == 0.0001
    @test EIOPA._round_va(-0.00005) == -0.0001
    @test EIOPA._round_rate(0.000005) == 0.00001
    @test EIOPA._round_rate(-0.000005) == -0.00001
    floor = EIOPA.ShockParameters(
        [1, 150]; up_relative = [0.0, 0.0], up_absolute = [0.0, 0.0],
        down_relative = [0.0, 0.0], down_absolute = [0.5, 0.5], floor = [-0.0123456, -0.0123456]
    )
    pub = EIOPA.publication(i; shock_parameters = floor)
    @test all(==(-0.0123456), pub.rates.down)
    @test EIOPA._annual(EIOPA.calculate(i; shock_parameters = floor).curves.down, 1) ≈ -0.0123456
    @test_throws ArgumentError EIOPA.RFRInputs([0.02, 0.03], [2, 2]; fsp = 2, llfr_weights = [0.0, 1.0], ufr = Periodic(0.03, 1), α = 0.11)
    @test_throws ArgumentError input([NaN, 0.03])
    @test_throws ArgumentError EIOPA.publication(i; maturities = [2, 1])
    val(x) = discount(EIOPA.calculate(input(x)).curves.basic, 60)
    x = [0.02, 0.03]
    fd = [(val(x + 1.0e-6 * (1:2 .== k)) - val(x - 1.0e-6 * (1:2 .== k))) / 2.0e-6 for k in 1:2]
    @test ForwardDiff.gradient(val, x) ≈ fd rtol = 1.0e-7
    @test_throws ArgumentError ForwardDiff.derivative(x -> EIOPA.publication(input([x, 0.03])).rates.basic[60], 0.02)
end
