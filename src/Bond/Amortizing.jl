# Amortizing loans and their payment schedules. Included in `module Bond`.

import ..__primal, ..__evalpoly_exact

"""
    Bond.Amortizing(rate, times, principal; accruals = diff([0; times]))

A loan or bond that repays `principal[i]` at `times[i]` and, with it, interest at the simple
annual `rate` on the balance outstanding since the previous payment: `rate * accruals[i] *
balance`. The initial balance is `sum(principal)`, so a balloon is a large last entry and the
balance is zero after the last payment. A negative entry adds to the balance (negative
amortization).

`accruals[i]` is the year fraction of payment `i`'s interest, by default the time since the
previous payment (from 0 for the first); pass it to apply a day count. Times are positive and
strictly increasing; the vectors are copied.

See also [`Bond.LevelPayment`](@ref), which builds the level-payment schedule of a mortgage, and
[`Bond.Pool`](@ref).

# Examples

```julia
# 100 repaid in four annual installments of 25, at 5%
loan = Bond.Amortizing(0.05, [1.0, 2.0, 3.0, 4.0], [25.0, 25.0, 25.0, 25.0])
collect(loan)   # 30.0, 28.75, 27.5, 26.25 at times 1, 2, 3, 4

# interest only, with a balloon at maturity
Bond.Amortizing(0.05, [1.0, 2.0, 3.0], [0.0, 0.0, 100.0])
```
"""
struct Amortizing{R <: Real, T <: Real, P <: Real, A <: Real} <: AbstractBond
    rate::R
    times::Vector{T}
    principal::Vector{P}
    accruals::Vector{A}
    function Amortizing(rate::R, times::AbstractVector{<:Real}, principal::AbstractVector{<:Real}, accruals::AbstractVector{<:Real}) where {R <: Real}
        n = length(times)
        n >= 1 || throw(ArgumentError("an amortizing schedule needs at least one payment"))
        (length(principal) == n && length(accruals) == n) || throw(
            DimensionMismatch("an amortizing schedule needs one principal amount and one accrual per payment time; got $n times, $(length(principal)) amounts, $(length(accruals)) accruals")
        )
        ts = collect(times)
        for (i, t) in enumerate(ts)
            p, prev = __primal(t), i == 1 ? zero(__primal(t)) : __primal(ts[i - 1])
            (isfinite(p) && p > prev) ||
                throw(ArgumentError("amortizing payment times must be finite, positive and strictly increasing; got $ts"))
        end
        ps, δs = collect(principal), collect(accruals)
        return new{R, eltype(ts), eltype(ps), eltype(δs)}(rate, ts, ps, δs)
    end
end
Amortizing(rate::Real, times, principal; accruals = diff([zero(eltype(times)); times])) = Amortizing(rate, times, principal, accruals)

maturity(b::Amortizing) = last(b.times)

"""
    Bond.LevelPayment(rate, frequency, maturity)

The [`Bond.Amortizing`](@ref) schedule of a unit loan repaid in equal payments `frequency` times a
year until `maturity`, at the nominal annual `rate` compounded at that frequency (a mortgage):
each payment is `i / (1 - (1 + i)^-n)` for `i = rate / frequency` and `n` payments. Interest
accrues `1 / frequency` of a year per payment, and the last payment repays the remaining balance,
so the schedule ends at exactly zero. A zero rate repays `1 / n` per payment, and rates near zero
keep their precision and their derivatives.

`frequency` is an integer or a `Periodic`, and `rate` a number or a `Periodic` rate of that
frequency. `maturity` must be a whole number of periods.

```julia
mortgage = Bond.LevelPayment(0.06, 12, 30)   # 360 monthly payments of about 0.0059955
```
"""
function LevelPayment(rate, frequency, maturity)
    f = __frequency(frequency).frequency
    r = __nominal(rate, frequency)
    __regular_schedule(maturity, f) ||
        throw(ArgumentError("Bond.LevelPayment needs a whole number of payment periods: maturity $maturity at frequency $f"))
    times = collect(coupon_times(maturity, f))
    n = length(times)
    δ = 1 / f
    payment = __level_payment(r * δ, n)
    principal = Vector{typeof(payment)}(undef, n)
    balance = one(payment)
    for k in 1:n
        # interest as `Amortizing` computes it, so that the payments are level
        p = k == n ? balance : payment - r * δ * balance
        principal[k] = p
        balance -= p
    end
    return Amortizing(r, times, principal, fill(δ, n))
end

# The level payment of a unit loan over `n` periods at the periodic rate `i`, i / (1 - (1+i)⁻ⁿ),
# written as (1/n)·Q(i)·B(n·log1p(i)) with Q(x) = x / log1p(x) and B(y) = y / (1 - e⁻ʸ). Both factors
# are 1 at 0, where each formula is 0/0; near 0 they are their series (Gregory coefficients and
# Bernoulli numbers), accurate to about 1e-20 below |x| = 1e-2, so a zero rate and its derivatives
# are exact limits.
const __GREGORY = (1 // 1, 1 // 2, -1 // 12, 1 // 24, -19 // 720, 3 // 160, -863 // 60480, 275 // 24192, -33953 // 3628800)
const __BERNOULLI_PLUS = (1 // 1, 1 // 2, 1 // 12, 0 // 1, -1 // 720, 0 // 1, 1 // 30240, 0 // 1, -1 // 1209600)
__x_over_log1p(x) = abs(x) < 1.0e-2 ? __evalpoly_exact(x, __GREGORY) : x / log1p(x)
__y_over_expm1(y) = abs(y) < 1.0e-2 ? __evalpoly_exact(y, __BERNOULLI_PLUS) : y / -expm1(-y)
__level_payment(i, n) = __x_over_log1p(i) * __y_over_expm1(n * log1p(i)) / n

"""
    Bond.__schedule(bond)

The payment schedule of a bond with principal, one named tuple per payment:
`(; time, accrual, opening, interest, principal, closing)`, with the balance before (`opening`)
and after (`closing`) the payment. Balances are accumulated from the last payment backward, so
the last closing balance is exactly zero.
"""
function __schedule(b::Amortizing)
    n = length(b.times)
    B = Vector{eltype(b.principal)}(undef, n + 1)
    B[n + 1] = zero(eltype(b.principal))
    for i in n:-1:1
        B[i] = B[i + 1] + b.principal[i]
    end
    return [
        (; time = b.times[i], accrual = b.accruals[i], opening = B[i], interest = b.rate * b.accruals[i] * B[i], principal = b.principal[i], closing = B[i + 1])
            for i in 1:n
    ]
end

# A fixed bond: unit principal until maturity, its coupons as its cashflow projection pays them
# (a short first stub accrues its actual length).
function __schedule(b::Fixed)
    ts = coupon_times(b)
    f = b.frequency.frequency
    regular = __regular_schedule(b.maturity, f)
    one_ = one(b.coupon_rate)
    return map(eachindex(ts)) do i
        t = ts[i]
        accrual = (i == 1 && !regular) ? oftype(float(t), t) : oftype(float(t), 1 / f)
        interest = (i == 1 && !regular) ? b.coupon_rate * t : b.coupon_rate / f
        last_payment = i == lastindex(ts)
        (; time = t, accrual, opening = one_, interest, principal = last_payment ? one_ : zero(one_), closing = last_payment ? zero(one_) : one_)
    end
end
