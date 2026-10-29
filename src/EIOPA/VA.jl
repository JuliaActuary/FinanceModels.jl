"""
    EIOPA.cssr(assets_pvbp, liabilities_pvbp)

Statutory credit spread sensitivity ratio: clamp the PVBP ratio to [0,1], or
return one when liability PVBP is nonpositive (Article 51a). Both inputs must be
finite. PVBP inputs use the same units and sign convention (value decreases for
a positive spread increase). This exception is omitted from RFR TD §9.3.2.
"""
function cssr(assets_pvbp, liabilities_pvbp)
    _require(_finite((assets_pvbp, liabilities_pvbp)), "PVBP values must be finite")
    liabilities_pvbp <= 0 && return one(assets_pvbp + liabilities_pvbp)
    return clamp(assets_pvbp / liabilities_pvbp, 0, 1)
end

function _cssr(v::VAInputs, notional)
    v.cssr !== nothing && return v.cssr
    v.sensitivities !== nothing && return cssr(v.sensitivities.assets, v.sensitivities.liabilities)
    v.valuations === nothing && return nothing
    _require(!iszero(notional), "zero notional VA requires supplied sensitivities or CSSR")
    a = v.valuations
    return cssr((a.assets - a.assets_notional) / notional, (a.liabilities - a.liabilities_notional) / notional)
end

function _va_amount(v::VAInputs, sensitivity, round_va)
    currency = 0.85 * sensitivity * v.rcs_currency * (v.usa === nothing ? 1 : v.usa)
    macro_va = if v.use_macro_va
        ω = clamp((v.rcs_country * v.debt_to_assets - 0.006) / 0.003, 0, 1)
        0.85 * sensitivity * max(v.rcs_country - 1.3 * v.rcs_currency, 0) * ω
    else
        zero(currency)
    end
    return currency, macro_va, round_va(currency + macro_va)
end

function _va_curve(prefix, llfr, adjustment, i; ufr = i.ufr)
    _require(isfinite(adjustment) && adjustment > -1, "VA must be finite and greater than -1")
    shift = log1p(adjustment)
    # A parallel continuous-zero shift through the FSP, but only the first
    # weighted term changes in the LLFR (RFR TD §10.2.1).
    base = ZeroRateCurve(log1p.(prefix) .+ shift, 1:i.fsp, Spline.ConstantForward())
    wf = i.llfr_weights[findfirst(==(i.fsp), i.tenors)]
    return _extrapolate(base, i, llfr + wf * shift; ufr)
end
