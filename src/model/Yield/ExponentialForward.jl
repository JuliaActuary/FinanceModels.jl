"""
    Yield.ExponentialForward(base; start, initial_forward::Rate, ufr::Rate, α)

Preserve `base` through `start` and extrapolate with instantaneous forward
`ufr + (initial_forward - ufr) * exp(-α * (t-start))` afterwards, in continuous
compounding. The discount factor is continuous; the forward at `start` is the
right-hand value `initial_forward`. `α` and `start` must be finite and positive.
The base must have normalized, positive discount factors on its retained domain.
Supports negative rates, first-order sensitivities, and exact long-horizon limits.
"""
struct ExponentialForward{C, T, F, U, A, L} <: AbstractYieldModel
    base::C
    start::T
    initial_forward::F
    ufr::U
    α::A
    _log_start::L
    function ExponentialForward(base::AbstractYieldModel; start, initial_forward::FinanceCore.Rate, ufr::FinanceCore.Rate, α)
        all(isfinite, (start, α)) && __primal(start) > 0 && __primal(α) > 0 ||
            throw(ArgumentError("ExponentialForward requires finite positive start and α"))
        l, u = __continuous(initial_forward), __continuous(ufr)
        all(isfinite, (l, u)) || throw(ArgumentError("ExponentialForward requires finite forward rates"))
        log_start = __log_discount(base, start)
        __log_discount(base, zero(start)) == 0 && isfinite(log_start) ||
            throw(ArgumentError("ExponentialForward requires normalized positive base discount factors"))
        return new{typeof(base), typeof(start), typeof(l), typeof(u), typeof(α), typeof(log_start)}(base, start, l, u, α, log_start)
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
    __primal(t) < __primal(c.start) && return __log_discount(c.base, t)
    return c._log_start + __exponential_interval(c, c.start, t)
end
function __log_interval(c::ExponentialForward, a, b)
    __check_time(a, "ExponentialForward")
    __check_time(b, "ExponentialForward")
    a == b && return zero(c.ufr)
    __primal(b) < __primal(a) && return -__log_interval(c, b, a)
    __at_infinity(b) && return __log_discount(c, b) - __log_discount(c, a)
    __primal(b) < __primal(c.start) && return __log_interval(c.base, a, b)
    __primal(a) >= __primal(c.start) && return __exponential_interval(c, a, b)
    return __log_interval(c.base, a, c.start) + __exponential_interval(c, c.start, b)
end
FinanceCore.discount(c::ExponentialForward, t) = exp(-__log_discount(c, t))
FinanceCore.discount(c::ExponentialForward, a, b) = exp(-__log_interval(c, a, b))
function force_of_interest(c::ExponentialForward, t)
    __check_time(t, "ExponentialForward")
    __primal(t) < __primal(c.start) && return force_of_interest(c.base, t)
    return c.ufr + (c.initial_forward - c.ufr) * exp(-c.α * (t - c.start))
end
function Base.zero(c::ExponentialForward, t)
    __at_infinity(t) && return Continuous(c.ufr)
    __primal(t) < __primal(c.start) && return zero(c.base, t)
    return __zero_from_log_discount(c, t)
end
function __log_tail(c::ExponentialForward)
    a0 = c._log_start - c.ufr * c.start + (c.initial_forward - c.ufr) / c.α
    return __LogTail(zero(a0), c.ufr, a0)
end

# Fields store continuous numbers, while the validating constructor takes Rates.
# Reconstruction re-enters that constructor and refreshes the cached anchor.
Base.propertynames(c::ExponentialForward, private::Bool = false) =
    private ? fieldnames(typeof(c)) : (:base, :start, :initial_forward, :ufr, :α)
Accessors.ConstructionBase.getproperties(c::ExponentialForward) =
    (; c.base, c.start, c.initial_forward, c.ufr, c.α)
function __reconstruct_exponential(; base, start, initial_forward, ufr, α)
    return ExponentialForward(base; start, initial_forward = Continuous(initial_forward), ufr = Continuous(ufr), α)
end
Accessors.ConstructionBase.setproperties(c::ExponentialForward, patch::NamedTuple) =
    __reconstruct_exponential(; merge(Accessors.ConstructionBase.getproperties(c), patch)...)
Accessors.ConstructionBase.constructorof(::Type{<:ExponentialForward}) =
    (base, start, initial_forward, ufr, α, _...) -> __reconstruct_exponential(; base, start, initial_forward, ufr, α)
