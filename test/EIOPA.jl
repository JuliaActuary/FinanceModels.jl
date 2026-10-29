using TOML

@testset "Official EIOPA May 2026 workbook" begin
    fixture = TOML.parsefile(joinpath(@__DIR__, "data", "eiopa_2026.toml"))
    p, v, expected = fixture["inputs"], fixture["va"], fixture["expected"]
    i = EIOPA.RFRInputs(
        p["rates"], p["tenors"]; fsp = p["fsp"], llfr_weights = p["llfr_weights"],
        ufr = Periodic(p["ufr"], 1), α = p["alpha"], cra_bp = p["cra_bp"], frequency = p["frequency"]
    )
    va = EIOPA.VAInputs(;
        rcs_currency = v["rcs_currency"], usa = v["usa"],
        valuations = (;
            assets = v["assets"], assets_notional = v["assets_notional"],
            liabilities = v["liabilities"], liabilities_notional = v["liabilities_notional"],
        )
    )
    shocks = EIOPA.ShockParameters()
    result = EIOPA.publication(i; va_inputs = va, shock_parameters = shocks)
    precise = EIOPA.calculate(i; va_inputs = va, shock_parameters = shocks)
    @test rate.(zero.(Ref(precise.bootstrapped), 1:50)) ≈ expected["bootstrap"] atol = 1.0e-13
    for name in (:llfr, :up_llfr, :down_llfr)
        @test getproperty(result.diagnostics, name) ≈ expected[string(name)] atol = 1.0e-13
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
