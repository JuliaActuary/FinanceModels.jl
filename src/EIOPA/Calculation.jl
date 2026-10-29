# Excel ROUND: decimal rounding with ties away from zero, rather than Julia's
# default ties-to-even. Publication inputs must not silently lose AD partials.
_round_rate(x) = round(x, RoundNearestTiesAway; digits = 5)
_round_va(x) = round(x, RoundNearestTiesAway; digits = 4)

function _calculate(i::RFRInputs, v, p, round_rate, round_va)
    bootstrap = _bootstrap(i)
    llfr = _llfr(bootstrap, i)
    basic = _extrapolate(bootstrap, i, llfr)
    cssr_curve = va_curve = nothing
    notional = sensitivity = currency = macro_va = total_va = nothing
    if v !== nothing
        prefix = [round_rate(_annual(basic, t)) for t in 1:i.fsp]
        notional = round_va(0.85 * v.rcs_currency)
        cssr_curve = _va_curve(prefix, llfr, notional, i)
        sensitivity = _cssr(v, notional)
        if sensitivity !== nothing
            currency, macro_va, total_va = _va_amount(v, sensitivity, round_va)
            va_curve = _va_curve(prefix, llfr, total_va, i)
        end
    end
    up = down = up_va = down_va = up_llfr = down_llfr = nothing
    if p !== nothing
        up, up_va, up_llfr = _shocked(bootstrap, i, p, 1, total_va)
        down, down_va, down_llfr = _shocked(bootstrap, i, p, -1, total_va)
    end
    curves = (; basic, cssr = cssr_curve, va = va_curve, up, down, up_va, down_va)
    diagnostics = (; llfr, cssr = sensitivity, notional_va = notional, currency_va = currency, macro_va, total_va, up_llfr, down_llfr)
    return Result(bootstrap, curves, diagnostics)
end

"""
    EIOPA.calculate(inputs::RFRInputs; va_inputs=nothing, shock_parameters=nothing)

Build full-precision basic, CSSR, VA and shocked curves. Optional stages run only
when their inputs are supplied. `VAInputs` containing only spreads returns the
CSSR curve for a second valuation pass. Downward floors apply to the final curves.
No publication rounding is applied; use `publication` for workbook reproduction.
First-order quote and parameter sensitivities propagate through the precise path,
away from clipping/floor boundaries. Maturities and FSP are fixed inputs.
"""
function calculate(i::RFRInputs; va_inputs::Union{Nothing, VAInputs} = nothing, shock_parameters::Union{Nothing, ShockParameters} = nothing)
    r = _calculate(i, va_inputs, shock_parameters, identity, identity)
    shock_parameters === nothing && return r
    curves = merge(r.curves, (; down = _floored(r.curves.down, shock_parameters), down_va = _floored(r.curves.down_va, shock_parameters)))
    return Result(r.bootstrapped, curves, r.diagnostics)
end

_has_dual(x) = false
_has_dual(x::ForwardDiff.Dual) = true
_has_dual(x::Rate) = _has_dual(rate(x))
_has_dual(x::Union{AbstractArray, Tuple, NamedTuple}) = any(_has_dual, x)
_has_dual(x::Union{RFRInputs, VAInputs, ShockParameters}) = any(k -> _has_dual(getfield(x, k)), fieldnames(typeof(x)))

"""
    EIOPA.publication(inputs::RFRInputs; va_inputs=nothing, shock_parameters=nothing,
        maturities=1:150)

Reproduce the May 2026 demonstration workbook: basic rates rounded to five decimal
places feed the CSSR/VA stage; notional and total VA are rounded to integer basis
points. Shocks use unrounded bootstrapped annual rates and their own LLFRs. Output
rates are rounded to five decimals, then downward floors are applied last.
Returns a `Publication`, with unavailable optional series set to `nothing`.
`maturities` must be increasing positive integers. AD inputs are rejected because
the prescribed rounding is nondifferentiable.
"""
function publication(i::RFRInputs; va_inputs::Union{Nothing, VAInputs} = nothing, shock_parameters::Union{Nothing, ShockParameters} = nothing, maturities = 1:150)
    _require(!_has_dual((i, va_inputs, shock_parameters)), "publication rounding does not support automatic differentiation; use calculate")
    t = collect(maturities)
    _require(!isempty(t) && _finite(t) && all(x -> x > 0 && isinteger(x), t) && all(>(0), diff(t)), "publication maturities must be increasing positive integers")
    r = _calculate(i, va_inputs, shock_parameters, _round_rate, _round_va)
    output(c) = c === nothing ? nothing : [_round_rate(_annual(c, m)) for m in t]
    series = map(output, r.curves)
    if shock_parameters !== nothing
        floors = [_floor(shock_parameters, m) for m in t]
        series = merge(series, (; down = max.(series.down, floors), down_va = series.down_va === nothing ? nothing : max.(series.down_va, floors)))
    end
    return Publication(ReadOnlyVector(Int.(t)), map(x -> x === nothing ? nothing : ReadOnlyVector(x), series), r.diagnostics)
end
