# Credit and Amortization

FinanceModels values loans and bonds whose principal amortizes, pools of loans that prepay, and
the default risk of their issuers. The assumptions a contract needs beyond the discount curve,
such as a prepayment speed or an issuer's hazard curve, live in a [`Models`](@ref
FinanceModels.Models) store under a key the contract carries, as a floating-rate bond reads its
index curve.

## Amortizing loans

[`Bond.Amortizing`](@ref FinanceModels.Bond.Amortizing) takes an annual interest rate, payment
times and the principal repaid at each. Interest accrues simply on the balance outstanding since
the previous payment, and the initial balance is the sum of the principal amounts:

```@example fixed_income
using FinanceModels

loan = Bond.Amortizing(0.05, [1.0, 2.0, 3.0, 4.0], [25.0, 25.0, 25.0, 25.0])
collect(loan)
```

A balloon is a large last principal amount, and `accruals` sets each payment's year fraction
when a day count applies. [`Bond.LevelPayment`](@ref FinanceModels.Bond.LevelPayment) builds the
schedule of a mortgage-style loan with equal payments:

```@example fixed_income
mortgage = Bond.LevelPayment(0.06, 12, 30)   # 30 years of monthly payments at 6%
first(collect(mortgage), 3)
```

## Loan pools and prepayment

[`Bond.Pool`](@ref FinanceModels.Bond.Pool) is a homogeneous pool of loans that share one
schedule, per unit of the pool's initial balance. At each payment the surviving loans pay as
scheduled, and then a fraction of their remaining balance prepays. The prepayment assumption is
read from the valuation context under the pool's key:

- [`Prepayment.CPR`](@ref FinanceModels.Prepayment.CPR): an annual conditional prepayment rate,
  `1 - (1 - CPR)^Δt` over a payment interval of `Δt` years;
- [`Prepayment.SMM`](@ref FinanceModels.Prepayment.SMM): a single monthly mortality rate, for
  monthly payments.

Either takes one rate for every payment or a vector with one rate per payment, which is how a
ramp such as a PSA-style speed is entered.

```@example fixed_income
pool = Bond.Pool(mortgage, "speed")
curve = Yield.Constant(0.05)
ctx = Models(curve, Dict("speed" => Prepayment.CPR(0.08)))
present_value(ctx, pool), present_value(curve, mortgage)
```

The 6% loans are worth more than par on a 5% curve, so prepayment at par lowers the pool's value.
The pool's cashflows are the projection under the same context:

```@example fixed_income
first(collect(Projection(pool, ctx)), 3)
```

Prepayment shrinks the share of the pool that survives; the remaining loans keep their schedule.
Prepayment rates are prescribed: they don't respond to interest rates.

## Hazard curves

A hazard curve gives the intensity of default of an issuer that has survived so far.
[`Credit.ConstantHazard`](@ref FinanceModels.Credit.ConstantHazard) has one intensity;
[`Credit.PiecewiseHazard`](@ref FinanceModels.Credit.PiecewiseHazard) is constant between tenors
and flat beyond the last:

```@example fixed_income
h = Credit.PiecewiseHazard([0.01, 0.02, 0.03], [1.0, 5.0, 10.0])
Credit.survival(h, 7.0), Credit.default_probability(h, 7.0, 8.0)
```

[`Credit.survival`](@ref FinanceModels.Credit.survival),
[`Credit.default_probability`](@ref FinanceModels.Credit.default_probability),
[`Credit.cumhazard`](@ref FinanceModels.Credit.cumhazard) and
[`Credit.hazard`](@ref FinanceModels.Credit.hazard) take one time (from 0) or an interval. Default
probabilities are computed as `-expm1(-cumhazard)`, so small ones keep their digits.

A [`Credit.Assumptions`](@ref FinanceModels.Credit.Assumptions) pairs a hazard curve with a
recovery fraction of the principal outstanding at default. It is the store entry that defaultable
contracts read by key.
