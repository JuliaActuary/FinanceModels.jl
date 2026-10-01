## Generic and Fallbacks
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
    par(curve,time;frequency=2)

Calculate the par yield for maturity `time` for the given `curve` and `frequency`. Returns a `Rate` object with periodicity corresponding to the `frequency`.

If `time` is shorter than one regular coupon period (e.g. `time=0.5` with `frequency=1`), the single stub payment implies a compounding frequency of `1/time`: the result is quoted as `Periodic(1/time)` when `1/time` is a (near-)integer, and otherwise an `ArgumentError` is thrown because the implied frequency cannot be represented as a `Periodic` rate.

If `time` is longer than one coupon period but not a whole number of periods, the schedule has a short first stub which accrues its actual length (see `Bond.coupon_times`); the result is the internal rate of return of the par-priced true-accrual schedule, quoted as `Periodic(frequency)` — so `par` of a flat curve recovers the curve's rate at any maturity. On such stub schedules this yield quote differs from the annualized par *coupon* `c` solving `c·Σᵢ δᵢ·DF(tᵢ) + DF(T) = 1`, which is what `InterestRateSwap` uses for its fixed leg.

# Examples

```julia-repl
julia> c = Yield.Constant(0.04);

julia> par(c,4)
Periodic(0.03960780543711406, 2)

julia> par(c,4;frequency=1)
Periodic(0.040000000000000036, 1)

julia> par(c,0.6;frequency=4)
Periodic(0.039413626195875295, 4)

julia> par(c,0.2;frequency=4)
Periodic(0.039374942589460726, 5)

julia> par(c,2.5)
Periodic(0.03960780543711406, 2)
```
"""
function par(curve, time; frequency = 2)
    coup_times = coupon_times(time, frequency)
    mat_disc = discount(curve, time)
    coupon_pv = sum(discount(curve, t) for t in coup_times)
    Δt = step(coup_times)
    r = (1 - mat_disc) / coupon_pv

    # A maturity that is not a whole number of periods leaves a short first stub
    # which accrues its actual length: such schedules pay `c·δᵢ` with the
    # annualized coupon `c` solving the true-accrual par condition, so the IRR
    # below is that of the correctly-accrued par-priced bond. Whole-period
    # schedules pay `r` every period (the historical construction, unchanged).
    stub = length(coup_times) > 1 && !__regular_schedule(time, frequency)
    c = stub ? __par_coupon(curve, time, frequency) : r

    # Build cash flows: initial outflow of -1, then coupons, final coupon+principal 1+coupon
    # Pre-allocate arrays for better performance
    n = length(coup_times)
    cfs = Vector{typeof(r)}(undef, n + 1)
    times = Vector{typeof(Δt)}(undef, n + 1)

    cfs[1] = -one(r)
    times[1] = zero(Δt)

    @inbounds for i in 1:n
        coup = stub ? c * (i == 1 ? first(coup_times) : 1 / frequency) : r
        cfs[i + 1] = i == n ? 1 + coup : coup
        times[i + 1] = coup_times[i]
    end

    r = FinanceCore.internal_rate_of_return(cfs, times)
    frequency_inner = 1 / Δt  # Simplified from min(1 / Δt, max(1 / Δt, frequency))
    if !isinteger(round(frequency_inner, digits = 8))
        throw(
            ArgumentError(
                "par(curve, $time; frequency=$frequency) implies a coupon period of $Δt and a compounding frequency of 1/Δt = $frequency_inner, which is not an integer and cannot be represented as a `Periodic` rate. Choose a maturity commensurate with the coupon frequency."
            )
        )
    end
    r = convert(Periodic(frequency_inner), r)
    return r
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
