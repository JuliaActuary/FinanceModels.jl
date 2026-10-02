# ── Internal curve capabilities ────────────────────────────────────────────────────────
#
# A curve needs only `discount(c, t)`. Each method below has a default that is safe for any such
# curve; a type adds its own where it knows more.
#
# | method                     | meaning                                        | default                        |
# |----------------------------|------------------------------------------------|--------------------------------|
# | `__log_discount(c, t)`     | L(t) = −log D(t)                               | `-log(discount(c, t))`         |
# | `__log_native(c)`          | L is direct, L(0) = 0, and D = exp(−L)         | `false`: the ratio D(b)/D(a)   |
# | `__log_interval(c, a, b)`  | L(b) − L(a), real where the interval factor is positive | from `__log_native`   |
# | `__log_tail(c)`            | the coefficients of L as t → ∞ (`__LogTail`)   | from `zero(c, Inf)`, below     |
#
# `__log_native` says that the log-difference route is valid for intervals: the curve's L is
# computed directly, with L(0) = 0, and its discount factor is exp(-L), so positive. It doesn't just
# mean that a type defines `__log_discount`: Smith–Wilson, `ForwardStarting`, Hull–White and the
# composite and scaled wrappers define L but take their intervals from `__log_interval`.
#
# - `__log_discount` + `__log_native`: `Constant`, `Spline`, `MonotoneConvex`, the zero-native curves
#   (`NelsonSiegel(Svensson)`, `CairnsPritchard(Extended)`, the yield shifts), `Vasicek`,
#   `CoxIngersollRoss`, `RatePath`.
# - `__log_discount` + `__log_interval`: `SmithWilson` (signed factors: from its interval ratio),
#   `ForwardStarting` (from the base curve's interval), `HullWhite` (its initial curve's),
#   `CompositeYield` and `ScaledYield` (their components' intervals, combined).
# - `__log_tail`, exact: `Constant`, `Spline`, `MonotoneConvex`, and the wrappers `CompositeYield`,
#   `ScaledYield`, `ForwardStarting`, `HullWhite`, which combine their components' tails. Every other
#   curve uses the fallback, which knows only the growth rate its zero rate at infinity implies.
# - A curve that forwards to another (`HullWhite`, `__PrimalCurve`) forwards every capability in
#   `__FORWARDED_CAPABILITIES` (at the end of this file), so the wrapped curve's own rules (its
#   intervals, its closed-form forward) apply. `__log_native` is for leaf curves only.

# Cumulative log-discount L(t) = −log D(t): the continuously compounded force accumulated from
# valuation time 0 to `t`. `forward` and the generic `zero` are built from it. The fallback needs
# only `discount(c, t)`; built-in curves compute L directly.
__log_discount(c, t) = -log(FinanceCore.discount(c, t))

# Whether a curve's discount factor is exp(-L) for a directly computed L, and hence positive. Such
# curves take interval factors as exp(L(from) − L(to)): a difference of log-discounts stays finite
# where both discount factors underflow, and an interval from 0 reproduces `discount(c, t)` bit for
# bit. A curve that defines only `discount` may have D(0) ≠ 1 or a nonpositive discount factor, so
# its interval factors stay the ratio D(to)/D(from). The answer depends only on the curve's type.
__log_native(c) = false

# L(to) − L(from), the log-discount accumulated over an interval; `ForwardStarting` rebases a curve
# with it. A curve that isn't log-native takes the log of its own interval factor, which is real
# wherever that factor is positive (for a signed Smith-Wilson fit, also where both discount factors
# are negative).
__log_interval(c, from, to) = __log_native(c) ? __log_native_interval(c, from, to) :
    -log(FinanceCore.discount(c, from, to))

# L(to) − L(from) of a log-native curve. L(0) = 0, so an interval from 0 (every present value) skips
# it and is L(to) exactly. Under ForwardDiff 1.x `iszero` also requires zero partials, so a time
# derivative at 0 still evaluates L(from).
function __log_native_interval(c, from, to)
    L_to = __log_discount(c, to)
    return L_to - (iszero(from) ? zero(L_to) : __log_discount(c, from))
end

# The tail of L (see `__LogTail`) for a curve without a closed form: the leading slope is its zero
# rate at infinity and the intercept is not known. When those rates are finite, a composite of such
# curves has the same limit as the sum of their zero rates, and an exact cancellation of the slopes
# is NaN (0·Inf). An infinite zero rate means L grows faster than t, at an order the curve doesn't
# say; the forward slot then holds ±Inf, which the limit reads as that unknown growth.
function __log_tail(c)
    z = FinanceCore.rate(Base.zero(c, Inf))
    return __LogTail(zero(z), z, oftype(z, NaN))
end

# Knot curves and wrappers evaluate L at t = +Inf from their tail (a wrapper's combines its
# components'); a negative time still reaches the components, which reject it. The limits are out
# of line so the finite-time path stays small.
__at_infinity(t) = isinf(t) && t > 0
@noinline __discount_at_infinity(c::AbstractYieldModel) = __discount_at_infinity(__log_tail(c))
@noinline __log_discount_at_infinity(c::AbstractYieldModel) = __log_discount_at_infinity(__log_tail(c))

# A time whose value is 0 but that carries partials: a time derivative taken at t = 0, where a zero
# rate L(t)/t is 0/0. It is `false` for a plain number, so the check compiles away (and keeps the
# zero rates small enough to inline) outside differentiation.
__dual_at_origin(t) = false
__dual_at_origin(t::ForwardDiff.Dual) = iszero(__primal(t)) && !iszero(t)

# L of a zero-native curve: z(t)·t with z the continuous zero rate, which is exactly the exponent
# `discount(zero(c, t), t)` uses. At an exact zero time L is 0 whatever the zero rate there: a
# wrapper over a discount-native curve has a removable 0/0 zero rate at t = 0. Under ForwardDiff
# 1.x `iszero` also requires zero partials, so a time derivative at 0 flows through z(t)·t, which
# needs the zero rate's limit at 0 and its derivatives: the built-in curves' own `zero` is exact
# there, and the generic `zero` of a curve that defines only its discount factor or L throws.
function __zero_log_discount(c, t)
    L = FinanceCore.rate(Base.zero(c, t)) * t
    return iszero(t) ? zero(L) : L
end

# Discount factor of the curves defined by their zero rate (NelsonSiegel(Svensson),
# CairnsPritchard, the yield shifts), shared through the one-line `discount`
# stubs at each curve definition. `exp(-L)` gives DF(0) = 1 exactly and keeps the promoted
# curve/time numeric type.
_discount_from_zero(c, t) = exp(-__zero_log_discount(c, t))

# Generic callable fallback: `curve(t) ≡ discount(curve, t)`. Covers every
# AbstractYieldModel subtype (Constant, Spline, CompositeYield, ScaledYield,
# TenorShift, ProjectedShift, NelsonSiegel, MonotoneConvex, …);
# each one routes through its own `discount`, so no per-type callable is needed.
(yc::AbstractYieldModel)(t) = FinanceCore.discount(yc, t)

# ── Public generic methods ─────────────────────────────────────────────────────────────

"""
    discount(yc, to)
    discount(yc, from,to)

The discount factor for the yield curve `yc` for times `from` through `to`.
"""
function FinanceCore.discount(yc::T, from, to) where {T <: AbstractYieldModel}
    # A log-native interval from 0 is `discount(yc, to)` exactly (see `__log_native_interval`).
    d = __log_native(yc) ? exp(-__log_native_interval(yc, from, to)) :
        FinanceCore.discount(yc, to) / FinanceCore.discount(yc, from)
    # The empty interval is the identity, also where L is infinite at both ends (from = to = Inf).
    # Under ForwardDiff 1.x `==` also compares partials, so this fires only where the
    # derivative is zero anyway.
    return from == to ? one(d) : d
end

"""
    forward(yc, from, to)

The forward `Rate` implied by the yield curve `yc` between times `from` and `to`.
"""
function FinanceCore.forward(yc::T, from, to = from + 1) where {T <: AbstractYieldModel}
    # forward = log(DF(from)/DF(to)) / (to-from): the interval's log-discount L(to) − L(from), by
    # the curve's own interval rule, per unit time.
    return Continuous(__log_interval(yc, from, to) / (to - from))
end

"""
    zero(curve,time)

Return the zero rate for the curve at the given time. At `time = 0` it is the limit of the zero
rate, the instantaneous forward rate at 0 (the short rate).
"""
function Base.zero(c::YC, time) where {YC <: AbstractYieldModel}
    # L/t; for a curve that defines only `discount`, L is -log(discount(c, time)). At t = 0 that is
    # 0/0, whose limit is L′(0), the instantaneous forward there. A time derivative at 0 would need
    # L″(0) as well, so it throws instead of returning a NaN.
    __dual_at_origin(time) && __throw_zero_rate_derivative_at_origin(time)
    iszero(time) && return Continuous(instantaneous_forward(c, time))
    return Continuous(__log_discount(c, time) / time)
end
@noinline __throw_zero_rate_derivative_at_origin(t) = throw(
    DomainError(
        t, "the zero rate of a curve without its own `zero` (L(t)/t, 0/0 at t = 0) has no derivative " *
            "in time at t = 0, and neither does a zero-rate transformation of it. Differentiate at a positive time."
    )
)

"""
    instantaneous_forward(curve, t)

The instantaneous (continuously compounded) forward rate of `curve` at time `t`,
``f(t) = -\\frac{d}{dt} \\log D(t)``. At `t = 0` it is the curve's short rate.

The built-in curves compute it in closed form. Any other curve differentiates its cumulative
log-discount with ForwardDiff, which stays finite where its discount factors underflow.

Note this is distinct from `forward(curve, from, to)`, which is the *discrete* forward `Rate`
between two times.
"""
instantaneous_forward(c::AbstractYieldModel, t) = __log_discount_derivative(c, t)
__log_discount_derivative(c, t) = ForwardDiff.derivative(s -> __log_discount(c, s), t)

"""
    accumulation(yc, from, to)

The accumulation factor for the yield curve `yc` for times `from` through `to`.
"""
function FinanceCore.accumulation(yc::AbstractYieldModel, time)
    return 1 ./ discount(yc, time)
end

# the reversed interval, so each curve's own interval method applies
FinanceCore.accumulation(yc::AbstractYieldModel, from, to) = FinanceCore.discount(yc, to, from)

# The capabilities a curve that forwards to another (`HullWhite`, `__PrimalCurve`) forwards, with
# their numbers of time arguments; the test suite checks every forwarding type against this list.
const __FORWARDED_CAPABILITIES = (
    (FinanceCore.discount, 1), (FinanceCore.discount, 2), (__log_discount, 1), (__log_interval, 2),
    (__log_tail, 0), (Base.zero, 1), (instantaneous_forward, 1),
)
