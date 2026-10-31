"""
    EIOPA.ShockParameters()
    EIOPA.ShockParameters(tenors; up_relative, up_absolute, down_relative,
        down_absolute, floor, ufr_shift_bp=15)

Versioned 2027 standard-formula shocks (Articles 165–166), or caller-supplied
piecewise-linear parameters. Relative shocks are fractions; absolute shocks and
downward floors are annual-effective decimals. Beyond the parameter grid its
endpoint values are held constant. All vectors are copied and read-only.
"""
struct ShockParameters{T, U, A, D, B, F, S}
    tenors::T
    up_relative::U
    up_absolute::A
    down_relative::D
    down_absolute::B
    floor::F
    ufr_shift_bp::S
    function ShockParameters(tenors; up_relative, up_absolute, down_relative, down_absolute, floor, ufr_shift_bp = 15)
        t = collect(tenors)
        vs = map(collect, (up_relative, up_absolute, down_relative, down_absolute, floor))
        _require(!isempty(t) && _finite(t) && all(>(0), t) && all(>(0), diff(t)), "shock tenors must be finite, positive and strictly increasing")
        _require(all(v -> length(v) == length(t) && _finite(v), vs), "shock vectors must be finite and aligned to tenors")
        _require(all(v -> all(>=(0), v), vs[1:4]) && all(>(-1), vs[5]), "shocks must be nonnegative and floors greater than -1")
        _require(isfinite(ufr_shift_bp) && ufr_shift_bp >= 0, "UFR shock must be finite and nonnegative")
        tt = ReadOnlyVector(float.(t))
        u, a, d, b, f = map(v -> ReadOnlyVector(float.(v)), vs)
        return new{typeof(tt), typeof(u), typeof(a), typeof(d), typeof(b), typeof(f), typeof(ufr_shift_bp)}(tt, u, a, d, b, f, ufr_shift_bp)
    end
end

function ShockParameters()
    t = 1:150
    relative(v, m) = m <= 50 ? v[m] : m >= 90 ? 0.2 : v[50] + (0.2 - v[50]) * ((m - 50) / 40)
    absolute(v, m) = m <= 50 ? v[m] : m >= 60 ? 0.0 : v[50] * ((60 - m) / 10)
    floor(m) = m <= 7 ? -0.0125 : m >= 20 ? -0.00893 : -0.0125 + (0.0125 - 0.00893) * ((m - 7) / 13)
    return ShockParameters(
        t;
        up_relative = [relative(_UP_RELATIVE, m) for m in t], up_absolute = [absolute(_UP_ABSOLUTE, m) for m in t],
        down_relative = [relative(_DOWN_RELATIVE, m) for m in t], down_absolute = [absolute(_DOWN_ABSOLUTE, m) for m in t],
        floor = floor.(t)
    )
end

function _parameter(p::ShockParameters, v, t)
    Yield.__primal(t) < first(p.tenors) && return first(v) + zero(t)
    Yield.__primal(t) >= last(p.tenors) && return last(v) + zero(t)
    k = searchsortedlast(p.tenors, Yield.__primal(t))
    return v[k] + (v[k + 1] - v[k]) * ((t - p.tenors[k]) / (p.tenors[k + 1] - p.tenors[k]))
end
_floor(p, t) = _parameter(p, p.floor, t)

_validate_shocked_inputs(i::RFRInputs, p, s) =
    _require(s === nothing, "shocked_inputs are only supported with PublishedRFR")

function _validate_shocked_inputs(i::PublishedRFR, p, s)
    if s === nothing
        _require(p === nothing, "PublishedRFR shocks require shocked_inputs=(up=PublishedRFR(...), down=PublishedRFR(...)); stressed bases cannot be recovered from the basic prefix")
        return nothing
    end
    _require(p !== nothing, "shocked_inputs require shock_parameters for the final downward floors")
    _require(
        s isa NamedTuple && Set(keys(s)) == Set((:up, :down)) && all(x -> x isa PublishedRFR, values(s)),
        "shocked_inputs must contain up and down PublishedRFR bases with unfloored annual zero rates"
    )
    for (direction, j) in ((1, s.up), (-1, s.down))
        _require(
            j.fsp == i.fsp && j.w_fsp == i.w_fsp && j.α == i.α,
            "published stressed bases must share FSP, FSP weight and α with the basic inputs"
        )
        expected = rate(Periodic(1)(i.ufr)) + direction * p.ufr_shift_bp / 10000
        _require(
            isapprox(rate(Periodic(1)(j.ufr)), expected; atol = 1.0e-12, rtol = 0),
            "published stressed UFRs must already include the specified annual-effective UFR shift"
        )
    end
    return nothing
end

function _shocked(bootstrap, i::RFRInputs, p, direction, total_va, shocked_inputs = nothing)
    sr = direction == 1 ? p.up_relative : p.down_relative
    ab = direction == 1 ? p.up_absolute : p.down_absolute
    rates = [
        _annual(bootstrap, t) * (1 + direction * _parameter(p, sr, t)) + direction * _parameter(p, ab, t)
            for t in 1:last(i.tenors)
    ]
    _require(_finite(rates) && all(>(-1), rates), "shocked annual rates must be finite and greater than -1")
    base = _annual_curve(rates)
    llfr = _llfr(base, i)
    u = rate(Periodic(1)(i.ufr)) + direction * p.ufr_shift_bp / 10000
    _require(u > -1, "shocked UFR must be greater than -1")
    ufr = Periodic(u, 1)
    curve = _extrapolate(base, i, llfr; ufr)
    with_va = total_va === nothing ? nothing : _va_curve(rates[1:i.fsp], llfr, total_va, i; ufr)
    return curve, with_va, llfr
end

function _shocked(::Nothing, i::PublishedRFR, p, direction, total_va, shocked_inputs)
    j = direction == 1 ? shocked_inputs.up : shocked_inputs.down
    _, llfr, curve = _basic(j)
    with_va = total_va === nothing ? nothing : _va_curve(j.rates, llfr, total_va, j)
    return curve, with_va, llfr
end
