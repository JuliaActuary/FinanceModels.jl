"""
    Yield.ExponentialForward(base; start, initial_forward::Rate, ufr::Rate, α)

Preserve `base` through `start` and extrapolate with instantaneous forward
`ufr + (initial_forward - ufr) * exp(-α * (t-start))` afterwards, in continuous
compounding. The discount factor is continuous; the forward at `start` is the
right-hand value `initial_forward`. `α` and `start` must be finite and positive.
The base must have normalized, positive discount factors on its retained domain.
Supports negative rates, first-order sensitivities, and exact long-horizon limits.
"""
struct ExponentialForward{C, T, F, U, A} <: AbstractYieldModel
    base::C
    start::T
    initial_forward::F
    ufr::U
    α::A
    function ExponentialForward(base::AbstractYieldModel; start, initial_forward::FinanceCore.Rate, ufr::FinanceCore.Rate, α)
        all(isfinite, (start, α)) && __primal(start) > 0 && __primal(α) > 0 ||
            throw(ArgumentError("ExponentialForward requires finite positive start and α"))
        l, u = __continuous(initial_forward), __continuous(ufr)
        all(isfinite, (l, u)) || throw(ArgumentError("ExponentialForward requires finite forward rates"))
        __log_discount(base, zero(start)) == 0 && isfinite(__log_discount(base, start)) ||
            throw(ArgumentError("ExponentialForward requires normalized positive base discount factors"))
        return new{typeof(base), typeof(start), typeof(l), typeof(u), typeof(α)}(base, start, l, u, α)
    end
end

# Stable even when α*dt underflows: dividing expm1(-α*dt) by α alone
# would then return zero instead of dt, and lose the initial-forward limit.
function __decay_integral(α, dt)
    x = α * dt
    threshold = sqrt(eps(one(float(__primal(x)))))
    abs(__primal(x)) < threshold && return dt * (1 - x / 2 + x^2 / 6 - x^3 / 24)
    return -expm1(-x) / α
end

# Integrate over a tail interval directly, avoiding cancellation of large cumulative logs.
function __exponential_interval(c::ExponentialForward, a, b)
    h, dt = a - c.start, b - a
    return c.ufr * dt + (c.initial_forward - c.ufr) * exp(-c.α * h) * __decay_integral(c.α, dt)
end
function __log_discount(c::ExponentialForward, t)
    __check_time(t, "ExponentialForward")
    __at_infinity(t) && return __log_discount_at_infinity(c)
    t <= c.start && return __log_discount(c.base, t)
    return __log_discount(c.base, c.start) + __exponential_interval(c, c.start, t)
end
function __log_interval(c::ExponentialForward, a, b)
    __check_time(a, "ExponentialForward")
    __check_time(b, "ExponentialForward")
    a == b && return zero(c.ufr)
    b < a && return -__log_interval(c, b, a)
    __at_infinity(b) && return __log_discount(c, b) - __log_discount(c, a)
    b <= c.start && return __log_interval(c.base, a, b)
    a >= c.start && return __exponential_interval(c, a, b)
    return __log_interval(c.base, a, c.start) + __exponential_interval(c, c.start, b)
end
FinanceCore.discount(c::ExponentialForward, t) = exp(-__log_discount(c, t))
FinanceCore.discount(c::ExponentialForward, a, b) = exp(-__log_interval(c, a, b))
function force_of_interest(c::ExponentialForward, t)
    __check_time(t, "ExponentialForward")
    t < c.start && return force_of_interest(c.base, t)
    return c.ufr + (c.initial_forward - c.ufr) * exp(-c.α * (t - c.start))
end
Base.zero(c::ExponentialForward, t) = __at_infinity(t) ? Continuous(c.ufr) : __zero_from_log_discount(c, t)
function __log_tail(c::ExponentialForward)
    a0 = __log_discount(c.base, c.start) - c.ufr * c.start + (c.initial_forward - c.ufr) / c.α
    return __LogTail(zero(a0), c.ufr, a0)
end
