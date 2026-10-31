"""
    Yield.ConstantForward(rates, tenors; extrapolation=:flat_forward)

A knot curve with linear cumulative log-discounts and constant forwards between knots.
Rates follow `ZeroRateCurve` conventions: numbers are continuous, and `Rate`s are
converted. Owns copies of its inputs and supports `knot_rates`, `knot_tenors`,
`reconstruct`, and first-order differentiation through calibration. Knot tenors must
be positive. The short-end zero rate is flat; time derivatives at knots are from
the right. Supports the standard tail policies except `:extension`.
"""
struct ConstantForward{R, T, P, E} <: AbstractInterpolatedZeroCurve
    spline::Sp.ConstantForward
    _rates::ReadOnlyVector{R, R}
    tenors::ReadOnlyVector{T, T}
    extrapolation::P
    _tail::E
    function ConstantForward(g::KnotGrid; extrapolation = :flat_forward)
        __primal(first(g.tenors)) > 0 || throw(ArgumentError("ConstantForward requires positive knot tenors"))
        extrapolation === :extension && throw(ArgumentError("ConstantForward does not support :extension"))
        t, z = last(g.tenors), last(g.rates)
        f = __last_discrete_forward(g.rates, g.tenors)
        tail = __build_tail(extrapolation, t, z, () -> f, () -> (f - z) / t)
        return new{eltype(g.rates), eltype(g.tenors), typeof(extrapolation), typeof(tail)}(
            Sp.ConstantForward(), ReadOnlyVector(g.rates), ReadOnlyVector(g.tenors), extrapolation, tail
        )
    end
end

ConstantForward(rates, tenors; extrapolation = :flat_forward) =
    __build_public(Sp.ConstantForward(), KnotGrid(rates, tenors, Sp.ConstantForward()); extrapolation)
__build(::Sp.ConstantForward, interpolation, g::KnotGrid; extrapolation = :flat_forward) = ConstantForward(g; extrapolation)
__tail(c::ConstantForward) = c._tail
__extends(::ConstantForward) = false
__in_tail(c::ConstantForward, t) = __primal(t) >= c._tail.last_tenor
function force_of_interest(c::ConstantForward, t)
    __check_time(t, "instantaneous_forward")
    return __in_tail(c, t) ? __tail_forward(c._tail, t) : __interior_forward(c, t)
end

function __constant_forward_segment(c::ConstantForward, t)
    i = searchsortedlast(c.tenors, __primal(t))
    i == 0 && return (zero(first(c.tenors)), zero(first(c._rates)), first(c._rates))
    i = min(i, length(c.tenors) - 1)
    i == 0 && return (zero(first(c.tenors)), zero(first(c._rates)), first(c._rates))
    a, b = c.tenors[i], c.tenors[i + 1]
    La, Lb = a * c._rates[i], b * c._rates[i + 1]
    return a, La, (Lb - La) / (b - a)
end
function __interior_log_discount(c::ConstantForward, t)
    a, La, f = __constant_forward_segment(c, t)
    return La + (t - a) * f
end
__interior_zero(c::ConstantForward, t) =
    __primal(t) < first(c.tenors) ? first(c._rates) + zero(t) : __interior_log_discount(c, t) / t
__interior_forward(c::ConstantForward, t) = last(__constant_forward_segment(c, t)) + zero(t)
