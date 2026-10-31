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
    v.cssr !== nothing && return (v.cssr, nothing, nothing)
    if v.sensitivities !== nothing
        assets, liabilities = v.sensitivities.assets, v.sensitivities.liabilities
    elseif v.valuations !== nothing
        _require(!iszero(notional), "zero notional VA requires supplied sensitivities or CSSR")
        a = v.valuations
        assets = (a.assets - a.assets_notional) / notional
        liabilities = (a.liabilities - a.liabilities_notional) / notional
    else
        return nothing, nothing, nothing
    end
    return cssr(assets, liabilities), assets, liabilities
end

function _va_amount(v::VAInputs, sensitivity, round_va)
    # Match the worksheet order, including the position of USA before RCS.
    currency = 0.85 * sensitivity * (v.usa === nothing ? 1 : v.usa) * v.rcs_currency
    ω = v.rcs_country === nothing || v.debt_to_assets === nothing ? nothing :
        clamp((v.rcs_country * v.debt_to_assets - 0.006) / 0.003, 0, 1)
    macro_va = if v.use_macro_va
        0.85 * sensitivity * max(v.rcs_country - 1.3 * v.rcs_currency, 0) * ω
    else
        zero(currency)
    end
    return currency, macro_va, round_va(currency + macro_va), ω
end

"""
    EIOPA.asset_spread_shock(asset_curve, notional_va)

Apply the Article 51a(2)(c) parallel spread increase to an asset's valuation curve
at **all** maturities, including beyond the FSP. `notional_va` is an annual-effective
decimal spread; its continuous shift is `log1p(notional_va)`. Use each asset's own
valuation curve, including its existing spread. Liabilities instead use the
notional CSSR curve returned by `calculate` or `publication`: its VA impact fades
after FSP. For instruments with spread-dependent cashflows, apply the parallel
spread shock within the asset valuation model as well.
"""
function asset_spread_shock(asset_curve::Yield.AbstractYieldModel, notional_va::Real)
    _require(isfinite(notional_va) && notional_va > -1, "notional asset spread must be finite and greater than -1")
    return asset_curve + Yield.Constant(Periodic(notional_va, 1))
end

function _va_curve(prefix, llfr, adjustment, i; ufr = i.ufr)
    _require(isfinite(adjustment) && adjustment > -1, "VA must be finite and greater than -1")
    shift = log1p(adjustment)
    # A parallel continuous-zero shift through the FSP, but only the first
    # weighted term changes in the LLFR (RFR TD §10.2.1).
    base = ZeroRateCurve(log1p.(prefix) .+ shift, 1:i.fsp, Spline.ConstantForward())
    return _extrapolate(base, i, llfr + _fsp_weight(i) * shift; ufr)
end
