"""
Prescribed prepayment assumptions for loan pools ([`Bond.Pool`](@ref FinanceModels.Bond.Pool)).

A pool reads its assumption from a [`Models`](@ref FinanceModels.Models) store under its key. At
each payment, after the scheduled payment, the fraction `smm` of each surviving loan's remaining
balance prepays. Prepayment reduces the share of the pool that survives; it doesn't change the
schedule of the loans that remain.

- [`Prepayment.CPR`](@ref): annual conditional prepayment rates, converted to each payment
  interval `Δt` as `1 - (1 - CPR)^Δt`;
- [`Prepayment.SMM`](@ref): single monthly mortality rates, applied as given to pools that pay
  monthly.

Either takes one rate for every payment, or a vector with one rate per payment. The rates are
fixed: they don't depend on interest rates (no behavioral prepayment).
"""
module Prepayment
import ..AbstractModel
import ..__primal

export CPR, SMM

"""
    Prepayment.AbstractPrepayment

The supertype of prepayment assumptions. A type defines `Prepayment.__smm(p, i, Δt)`, the
fraction of the remaining balance that prepays at payment `i`, which ends an interval of
length `Δt`.
"""
abstract type AbstractPrepayment <: AbstractModel end

function __check_rates(x)
    for r in x
        p = __primal(r)
        0 <= p <= 1 || throw(DomainError(r, "a prepayment rate lies in [0, 1]"))
    end
    return x
end

"""
    Prepayment.CPR(rate)
    Prepayment.CPR(rates::AbstractVector)

An annual conditional prepayment rate: the fraction of a loan's remaining balance that
prepays over a year. Over a payment interval of length `Δt` (in years) the fraction
`1 - (1 - rate)^Δt` prepays at the payment. A vector gives one annual rate per payment.

```julia
pool = Bond.Pool(Bond.LevelPayment(0.06, 12, 30), "CPR")
present_value(Models(curve, Dict("CPR" => Prepayment.CPR(0.06))), pool)
```
"""
struct CPR{T} <: AbstractPrepayment
    rate::T
    CPR(rate::T) where {T <: Union{Real, AbstractVector{<:Real}}} = new{T}(__check_rates(rate))
end

"""
    Prepayment.SMM(rate)
    Prepayment.SMM(rates::AbstractVector)

A single monthly mortality rate: the fraction of a loan's remaining balance that prepays at a
monthly payment. A pool that doesn't pay monthly throws an `ArgumentError`; use
[`Prepayment.CPR`](@ref) for other payment intervals. A vector gives one rate per payment.
An `SMM` of `m` and a `CPR` of `1 - (1 - m)^12` prepay the same on a monthly schedule.
"""
struct SMM{T} <: AbstractPrepayment
    rate::T
    SMM(rate::T) where {T <: Union{Real, AbstractVector{<:Real}}} = new{T}(__check_rates(rate))
end

__rate(x::Real, i) = x
__rate(x::AbstractVector, i) = x[i]

# The fraction of the remaining balance that prepays at payment `i`, after an interval of `Δt`.
__smm(p::CPR, i, Δt) = -expm1(Δt * log1p(-__rate(p.rate, i)))
function __smm(p::SMM, i, Δt)
    abs(Δt - 1 / 12) <= 1.0e-8 || throw(
        ArgumentError("Prepayment.SMM applies to monthly payments; payment $i ends an interval of $Δt years. Use Prepayment.CPR.")
    )
    return __rate(p.rate, i)
end

# A vector of rates gives one per payment, so it must match the schedule.
__check_length(p::AbstractPrepayment, n) = __check_length(p.rate, n)
__check_length(::Real, n) = nothing
function __check_length(x::AbstractVector, n)
    length(x) == n ||
        throw(DimensionMismatch("a vector of prepayment rates needs one rate per payment: got $(length(x)) for $n payments"))
    return nothing
end

end
