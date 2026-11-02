_annual(c, t) = expm1(rate(zero(c, t)))

function _bootstrap(i::RFRInputs)
    r = i.rates .- i.cra_bp / 10000
    if i.instrument === :government
        return ZeroRateCurve(Periodic.(r, 1), i.tenors, Spline.ConstantForward())
    end
    quotes = [ParSwapYield(r[k], i.tenors[k]; frequency = i.frequency) for k in eachindex(r)]
    return fit(Spline.ConstantForward(), quotes, Fit.Bootstrap())
end

function _llfr(c, i::RFRInputs)
    return sum(eachindex(i.tenors)) do k
        t, w = i.tenors[k], i.llfr_weights[k]
        t < i.fsp && return zero(w)
        # Article 46(1c)(a) and the workbook's GetLLFR use the preceding DLT observation.
        # On a constant-forward bootstrap this equals the one-year forward into
        # FSP. After maturity-dependent shocks they differ; retain the workbook
        # convention rather than silently switching its LLFR calculation.
        a = t == i.fsp ? (k == 1 ? 0 : i.tenors[k - 1]) : i.fsp
        return w * rate(forward(c, a, t))
    end
end

_extrapolate(base, i, llfr; ufr = i.ufr) =
    Yield.ExponentialForward(base; start = i.fsp, initial_forward = Continuous(llfr), ufr, α = i.α)

# Annual nodes are also the documented interpolation convention for shocked curves.
_annual_curve(rates) = ZeroRateCurve(Periodic.(rates, 1), 1:length(rates), Spline.ConstantForward())

function _basic(i::RFRInputs)
    bootstrap = _bootstrap(i)
    llfr = _llfr(bootstrap, i)
    return bootstrap, llfr, _extrapolate(bootstrap, i, llfr)
end
_basic(i::PublishedRFR) = (nothing, i.llfr, _extrapolate(_annual_curve(i.rates), i, i.llfr))

_prefix(i::RFRInputs, basic, round_rate) = [round_rate(_annual(basic, t)) for t in 1:i.fsp]
_prefix(i::PublishedRFR, basic, round_rate) = round_rate.(i.rates)
_fsp_weight(i::RFRInputs) = i.llfr_weights[findfirst(==(i.fsp), i.tenors)]
_fsp_weight(i::PublishedRFR) = i.w_fsp
