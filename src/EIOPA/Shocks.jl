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
    t <= first(p.tenors) && return first(v) + zero(t)
    t >= last(p.tenors) && return last(v) + zero(t)
    k = searchsortedlast(p.tenors, Yield.__primal(t))
    return v[k] + (v[k + 1] - v[k]) * ((t - p.tenors[k]) / (p.tenors[k + 1] - p.tenors[k]))
end
_floor(p, t) = _parameter(p, p.floor, t)
_floored(c, p) = c === nothing ? nothing : Yield.TenorShift(c, (z, t) -> Periodic(max(expm1(rate(z)), _floor(p, t)), 1))

function _shocked(bootstrap, i, p, direction, total_va)
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
