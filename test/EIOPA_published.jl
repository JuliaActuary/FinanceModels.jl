using TOML

@testset "Published EIOPA basic inputs" begin
    f = TOML.parsefile(joinpath(@__DIR__, "data", "eiopa_2026.toml"))
    p, v, e = f["inputs"], f["va"], f["expected"]
    wf = p["llfr_weights"][findfirst(==(p["fsp"]), p["tenors"])]
    make(r = e["basic"][1:p["fsp"]]; kwargs...) = EIOPA.PublishedRFR(
        r;
        fsp = p["fsp"], llfr = Continuous(e["llfr"]), w_fsp = wf,
        ufr = Periodic(p["ufr"], 1), α = p["alpha"], kwargs...
    )
    i = make()
    va = EIOPA.VAInputs(
        rcs_currency = v["rcs_currency"], usa = v["usa"],
        valuations = (
            assets = v["assets"], assets_notional = v["assets_notional"],
            liabilities = v["liabilities"], liabilities_notional = v["liabilities_notional"],
        )
    )
    r = @inferred EIOPA.calculate(i; va_inputs = va)
    pub = @inferred EIOPA.publication(i; va_inputs = va)
    @test r.bootstrapped === nothing
    @test pub.rates.basic[1:i.fsp] == collect(i.rates)
    @test pub.rates.cssr ≈ e["cssr"] atol = 1.0e-14 rtol = 0
    @test pub.rates.va ≈ e["va"] atol = 1.0e-14 rtol = 0
    @test pub.diagnostics.llfr == e["llfr"]
    @test pub.diagnostics.cssr == e["cssr_value"]
    @test pub.diagnostics.total_va == e["total_va"]
    @test rate(Yield.instantaneous_forward(r.curves.basic, i.fsp)) == e["llfr"]
    u, h = log1p(p["ufr"]), 40 - i.fsp
    expected_log = i.fsp * log1p(last(i.rates)) + u * h + (i.llfr - u) * (-expm1(-i.α * h)) / i.α
    @test discount(r.curves.basic, 40) ≈ exp(-expected_log)
    @test r.curves.up === nothing

    rates = collect(i.rates)
    copied = make(rates)
    rates[1] += 0.01
    @test first(copied.rates) == first(i.rates)
    @test_throws Base.CanonicalIndexError copied.rates[1] = 0.5
    for bad in (fill(-1.0, i.fsp), fill(NaN, i.fsp), Float64[], zeros(i.fsp - 1), zeros(i.fsp + 1))
        @test_throws ArgumentError make(bad)
    end
    controls = (; fsp = 2, llfr = Continuous(0.03), w_fsp = 0.4, ufr = Periodic(0.033, 1), α = 0.11)
    for patch in (
            (; fsp = 0), (; fsp = 1.5), (; w_fsp = 1.01), (; w_fsp = -0.1),
            (; llfr = Continuous(Inf)), (; α = 0.0), (; α = NaN),
        )
        @test_throws ArgumentError EIOPA.PublishedRFR([0.02, 0.03]; merge(controls, patch)...)
    end
    @test_throws TypeError EIOPA.PublishedRFR([0.02, 0.03]; merge(controls, (; llfr = 0.03))...)
    @test_throws ArgumentError EIOPA.calculate(i; shock_parameters = EIOPA.ShockParameters())

    # No sensitivity is lost to a hidden bootstrap, CRA deduction or reconstructed LLFR.
    function value(x)
        rates = [k == 1 ? x[1] : i.rates[k] for k in eachindex(i.rates)]
        input = EIOPA.PublishedRFR(
            rates; fsp = i.fsp, llfr = Continuous(x[2]),
            w_fsp = x[3], ufr = Periodic(x[4], 1), α = x[5]
        )
        result = EIOPA.calculate(input; va_inputs = EIOPA.VAInputs(rcs_currency = x[6], cssr = 0.5))
        return discount(result.curves.va, 60) + discount(result.curves.basic, 1)
    end
    x = [first(i.rates), i.llfr, i.w_fsp, p["ufr"], i.α, 0.01]
    finite = [(value(x + 1.0e-6 * (1:6 .== k)) - value(x - 1.0e-6 * (1:6 .== k))) / 2.0e-6 for k in 1:6]
    @test ForwardDiff.gradient(value, x) ≈ finite rtol = 1.0e-7 atol = 1.0e-9
    @test_throws ArgumentError ForwardDiff.derivative(
        q -> EIOPA.publication(
            EIOPA.PublishedRFR(
                collect(i.rates); fsp = i.fsp, llfr = Continuous(q), w_fsp = i.w_fsp, ufr = i.ufr, α = i.α
            )
        ).rates.basic[60], i.llfr
    )
end

@testset "Published and market basic/CSSR/VA equivalence" begin
    v = EIOPA.VAInputs(rcs_currency = 0.01, cssr = 1.0)
    for level in (0.03, -0.005, -0.04)
        market = EIOPA.RFRInputs(
            [level, level], [1, 20]; instrument = :government,
            fsp = 20, llfr_weights = [0.0, 1.0], ufr = Periodic(0.033, 1), α = 0.11
        )
        raw = EIOPA.calculate(market)
        basic = EIOPA.PublishedRFR(
            EIOPA._annual.(Ref(raw.curves.basic), 1:20);
            fsp = 20, llfr = Continuous(raw.diagnostics.llfr), w_fsp = 1.0,
            ufr = Periodic(0.033, 1), α = 0.11
        )
        actual = @inferred EIOPA.calculate(basic; va_inputs = v)
        expected = EIOPA.calculate(market; va_inputs = v)
        for name in (:basic, :cssr, :va), t in (0.0, 1.0, 10.5, 20.0, 60.0, 150.0, Inf)
            @test discount(getproperty(actual.curves, name), t) ≈ discount(getproperty(expected.curves, name), t) rtol = 2.0e-14
        end
        pub = @inferred EIOPA.publication(basic; va_inputs = v)
        reference = EIOPA.publication(market; va_inputs = v)
        for name in (:basic, :cssr, :va)
            @test getproperty(pub.rates, name) ≈ getproperty(reference.rates, name) atol = 1.0e-14 rtol = 0
        end
        @test actual.bootstrapped === nothing
        for calculate in (EIOPA.calculate, EIOPA.publication)
            @test_throws ArgumentError calculate(basic; shock_parameters = EIOPA.ShockParameters())
        end
    end
end
