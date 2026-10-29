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
