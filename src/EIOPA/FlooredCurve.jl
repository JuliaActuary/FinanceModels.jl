# A final annual-zero floor. Unlike an opaque TenorShift, this view retains the
# exponential curve's exact asymptotic intercept, including finite limits at UFR=0.
struct _FlooredCurve{C <: Yield.ExponentialForward, P} <: Yield.AbstractYieldShift
    base::C
    parameters::P
end
_floored(c, p) = c === nothing ? nothing : _FlooredCurve(c, p)

function Base.zero(c::_FlooredCurve, t)
    Yield.__check_time(t, "EIOPA downward floor")
    return Continuous(max(rate(zero(c.base, t)), log1p(_floor(c.parameters, t))))
end

function Yield.__log_tail(c::_FlooredCurve)
    tail = Yield.__log_tail(c.base)
    floor = log1p(last(c.parameters.floor))
    u, f = Yield.__primal(tail.a1), Yield.__primal(floor)
    u > f && return tail
    u < f && return Yield.__LogTail(zero(floor), floor, zero(floor))
    # At a zero common slope, perturbing either slope can change a finite limit
    # into zero/Inf. Do not drop a derivative by selecting the constant floor.
    # ForwardDiff 1.x `!=` compares partials as well as primal values.
    if iszero(u) && tail.a1 != floor
        throw(ArgumentError("discount(floored curve, Inf) has no derivative where a zero tail slope can move; differentiate at a finite time"))
    end
    return Yield.__LogTail(zero(floor), tail.a1, max(tail.a0, zero(tail.a0)))
end
