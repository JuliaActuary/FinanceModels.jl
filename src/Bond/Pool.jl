# Homogeneous loan pools. Included in `module Bond`.

"""
    Bond.Pool(loan, key)

A homogeneous pool of loans that each follow `loan`'s schedule (a [`Bond.Amortizing`](@ref) or a
[`Bond.Fixed`](@ref)), per unit of the pool's initial balance, with prepayment read from a valuation
context under `key` (a [`Prepayment.CPR`](@ref FinanceModels.Prepayment.CPR) or
[`Prepayment.SMM`](@ref FinanceModels.Prepayment.SMM)).

At each payment the surviving loans pay their scheduled interest and principal, then the fraction
`smm` of their remaining balance prepays. With `F` the share of the pool that has not prepaid
before payment `i`, the pool pays `F · (interestᵢ + principalᵢ + smmᵢ · closingᵢ)` and `F` becomes
`F · (1 - smmᵢ)`. Prepayment shrinks the surviving share; it doesn't change the remaining loans'
schedule.

```julia
pool = Bond.Pool(Bond.LevelPayment(0.06, 12, 30), "CPR")
ctx = Models(Yield.Constant(0.05), Dict("CPR" => Prepayment.CPR(0.08)))
present_value(ctx, pool)
collect(Projection(pool, ctx))   # the pool's cashflows
```
"""
struct Pool{L <: AbstractBond, K} <: AbstractBond
    loan::L
    key::K
end

maturity(p::Pool) = maturity(p.loan)
