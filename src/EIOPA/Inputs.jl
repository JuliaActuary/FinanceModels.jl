_require(ok, message) = ok || throw(ArgumentError("EIOPA: " * message))
_finite(xs) = all(isfinite, xs)

"""
    EIOPA.RFRInputs(rates, tenors; fsp, llfr_weights, ufr::Rate, α,
        instrument=:swaps, frequency=1, cra_bp=0)

Complete, selected DLT market observations. `rates` are annual nominal swap coupons
for `:swaps`, or annual-effective government zero rates for `:government`.
`frequency` is the swap coupon frequency; CRA is in basis points. Supply positive,
increasing integer maturities, an observed FSP, and aligned nonnegative LLFR weights
summing to one (zero before the FSP). All vectors are copied and read-only.
Supply the assessed FSP explicitly; it is not inferred from positive weights.
An assessed FSP with zero weight is supported and has no VA contribution to LLFR.
UFR compounding is explicit; α is supplied, not calibrated. Currency adjustments,
if needed, must already be reflected in the observations supplied to this API.
"""
struct RFRInputs{R, T, W, U, A, C}
    rates::R
    tenors::T
    llfr_weights::W
    fsp::Int
    ufr::U
    α::A
    instrument::Symbol
    frequency::Int
    cra_bp::C
    function RFRInputs(rates, tenors; fsp, llfr_weights, ufr::Rate, α, instrument = :swaps, frequency = 1, cra_bp = 0)
        r, t, w = collect(rates), collect(tenors), collect(llfr_weights)
        _require(!isempty(t) && length(r) == length(t) == length(w), "rates, tenors and weights must have equal nonzero lengths")
        _require(_finite(t) && all(x -> x > 0 && isinteger(x), t) && all(>(0), diff(t)), "DLT maturities must be positive, strictly increasing integers")
        _require(fsp in t && isinteger(fsp), "FSP must be an observed DLT maturity")
        _require(_finite(w) && all(>=(0), w) && isapprox(sum(w), 1; atol = 1.0e-12, rtol = 0), "LLFR weights must be nonnegative and sum to one")
        _require(all(i -> t[i] >= fsp || iszero(w[i]), eachindex(t)), "LLFR weights must be zero before the FSP")
        _require(instrument in (:swaps, :government), "instrument must be :swaps or :government")
        _require(frequency isa Integer && frequency > 0, "coupon frequency must be a positive integer")
        _require(_finite(r) && isfinite(cra_bp) && isfinite(α) && α > 0 && isfinite(rate(Continuous(ufr))), "rates and parameters must be finite, with α positive")
        corrected = r .- cra_bp / 10000
        _require(all(>(-1), instrument === :swaps ? corrected ./ frequency : corrected), "adjusted rates must have positive compounding factors")
        rr, tt, ww = ReadOnlyVector(float.(r)), ReadOnlyVector(Int.(t)), ReadOnlyVector(float.(w))
        return new{typeof(rr), typeof(tt), typeof(ww), typeof(ufr), typeof(α), typeof(cra_bp)}(rr, tt, ww, Int(fsp), ufr, α, instrument, frequency, cra_bp)
    end
end

"""
    EIOPA.PublishedRFR(rates; fsp, llfr::Rate, w_fsp, ufr::Rate, α)

Import a 2027-methodology basic curve: `rates[t]` is the published annual-effective
zero rate at year `t`, for exactly `1:fsp`. Supply the official LLFR, FSP weight,
UFR and convergence parameter. LLFR and UFR have explicit compounding; the LLFR
published in continuous terms should be wrapped in `Continuous`.

No bootstrapping, CRA deduction or LLFR estimation is performed. The basic tail
is reconstructed from the supplied FSP rate and controls; rounding that FSP rate
can affect the last published digit beyond FSP. Inputs are copied and read-only.
Use the matching currency, date and 2027 methodology, rather than legacy
Smith–Wilson publications. Optional stressed bases are supplied separately to
`calculate`/`publication` through `shocked_inputs=(up=..., down=...)`.
"""
struct PublishedRFR{R, L, W, U, A}
    rates::R
    fsp::Int
    llfr::L
    w_fsp::W
    ufr::U
    α::A
    function PublishedRFR(rates; fsp, llfr::Rate, w_fsp, ufr::Rate, α)
        r = collect(rates)
        _require(
            isfinite(fsp) && isinteger(fsp) && fsp > 0 && length(r) == fsp,
            "published rates must cover exactly the integer years 1:FSP"
        )
        _require(_finite(r) && all(>(-1), r), "published annual zero rates must be finite and greater than -1")
        l, u = rate(Continuous(llfr)), rate(Continuous(ufr))
        _require(
            _finite((l, u, α, w_fsp)) && α > 0 && 0 <= w_fsp <= 1,
            "published controls must be finite, with α positive and FSP weight in [0,1]"
        )
        rr, w = ReadOnlyVector(float.(r)), float(w_fsp)
        return new{typeof(rr), typeof(l), typeof(w), typeof(ufr), typeof(α)}(rr, Int(fsp), l, w, ufr, α)
    end
end

"""
    EIOPA.VAInputs(; rcs_currency, rcs_country=nothing, debt_to_assets=nothing,
        use_macro_va=false, usa=nothing, valuations=nothing, sensitivities=nothing, cssr=nothing)

Spreads are decimal annual rates. `usa` is the undertaking-specific multiplier
(at most 1.05); it is mutually exclusive with macro VA. For CSSR supply one of:
`cssr`, `(assets=..., liabilities=...)` PVBP `sensitivities`, or four `valuations`
named `assets`, `assets_notional`, `liabilities`, `liabilities_notional`.
With spreads alone the calculation returns the notional CSSR curve for a second
valuation pass, leaving VA-dependent outputs unavailable (`nothing`).
`assets_notional` must use a parallel asset spread shock at all maturities (see
`asset_spread_shock`), while `liabilities_notional` uses the fading CSSR curve.
Euro-area undertakings without USA must select macro VA, except the peer-country
cases prescribed to have zero macro VA. Eligibility is supplied by the caller.
The notional VA uses unadjusted currency RCS even when USA is selected, matching
the demonstration workbook; USA affects the final currency VA.
"""
struct VAInputs{R, C, D, U, V, S, Q}
    rcs_currency::R
    rcs_country::C
    debt_to_assets::D
    use_macro_va::Bool
    usa::U
    valuations::V
    sensitivities::S
    cssr::Q
    function VAInputs(; rcs_currency, rcs_country = nothing, debt_to_assets = nothing, use_macro_va = false, usa = nothing, valuations = nothing, sensitivities = nothing, cssr = nothing)
        _require(isfinite(rcs_currency), "currency RCS must be finite")
        _require(usa === nothing || (isfinite(usa) && usa <= 1.05), "USA must be finite and at most 1.05")
        _require(!(use_macro_va && usa !== nothing), "macro VA and USA are mutually exclusive")
        _require(rcs_country === nothing || isfinite(rcs_country), "country RCS must be finite")
        _require(debt_to_assets === nothing || (isfinite(debt_to_assets) && 0 <= debt_to_assets <= 1), "debt-to-assets ratio must be in [0,1]")
        _require(!use_macro_va || (rcs_country !== nothing && debt_to_assets !== nothing && isfinite(rcs_country) && isfinite(debt_to_assets) && 0 <= debt_to_assets <= 1), "macro VA requires country RCS and a debt-to-assets ratio in [0,1]")
        _require(count(!isnothing, (valuations, sensitivities, cssr)) <= 1, "supply only one CSSR source")
        if valuations !== nothing
            _require(valuations isa NamedTuple && Set(keys(valuations)) == Set((:assets, :assets_notional, :liabilities, :liabilities_notional)) && _finite(values(valuations)), "valuations require four finite named amounts")
        end
        if sensitivities !== nothing
            _require(sensitivities isa NamedTuple && Set(keys(sensitivities)) == Set((:assets, :liabilities)) && _finite(values(sensitivities)), "sensitivities require finite assets and liabilities PVBP")
        end
        _require(cssr === nothing || (isfinite(cssr) && 0 <= cssr <= 1), "CSSR must be in [0,1]")
        return new{typeof(rcs_currency), typeof(rcs_country), typeof(debt_to_assets), typeof(usa), typeof(valuations), typeof(sensitivities), typeof(cssr)}(rcs_currency, rcs_country, debt_to_assets, use_macro_va, usa, valuations, sensitivities, cssr)
    end
end

"""
    EIOPA.Result

`bootstrapped` retains all DLT observations for `RFRInputs`, and is `nothing` for
`PublishedRFR`. `curves` has `basic`, `cssr`, `va`,
`up`, `down`, `up_va`, `down_va`; optional unavailable curves are `nothing`.
`diagnostics` exposes basic/adjusted LLFRs, CSSR, asset/liability PVBP (in the
workbook's per-unit-rate units for valuation inputs), country weight, and VA amounts.
"""
struct Result{B, C, D}
    bootstrapped::B
    curves::C
    diagnostics::D
end

"""
    EIOPA.Publication

Annual `maturities`, a `rates` named tuple with the same seven names as
`Result.curves`, and calculation `diagnostics`. Rates are decimal annual-effective
numbers. Downward floors are applied after rounding and can retain more digits.
"""
struct Publication{T, R, D}
    maturities::T
    rates::R
    diagnostics::D
end
