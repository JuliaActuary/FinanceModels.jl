# Credit, Amortization, and Optionality

FinanceModels values loans and bonds whose principal amortizes, pools of loans that prepay, the
default risk of their issuers, and the call and put rights on them. The assumptions a contract
needs beyond the discount curve, such as a prepayment speed or an issuer's hazard curve, live in a
[`Models`](@ref FinanceModels.Models) store under a key the contract carries, as a floating-rate
bond reads its index curve. [Worked Examples: a Callable Corporate Bond and a Mortgage Pool](@ref)
puts the pieces together.

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

## Defaultable contracts

[`Credit.Defaultable`](@ref FinanceModels.Credit.Defaultable) subjects an asset with principal (a
fixed bond, an amortizing loan or a pool) to its issuer's default, read from the context under
the issuer's key. On default the holder recovers the recovery fraction of the principal
outstanding immediately before it, without accrued interest, and receives nothing more:

```@example fixed_income
acme = Credit.Assumptions(Credit.ConstantHazard(0.02); recovery = 0.4)
ctx = Models(Yield.Constant(Continuous(0.04)), Dict("ACME" => acme))
bond = Bond.Fixed(0.05, Periodic(2), 10)
present_value(ctx, Credit.Defaultable(bond, "ACME")), present_value(Yield.Constant(Continuous(0.04)), bond)
```

The value is that of the expected cashflows: each payment weighted by the probability of
survival to it, and the expected recoveries between payments. Recovery accrues continuously, so
the projection lists it at the nodes of a fixed quadrature rule (4-point Gauss–Legendre on each
interval between payment dates and hazard tenors, in parts of at most a month), and discounting
those nodes on any curve is the integral:

```@example fixed_income
cfs = collect(Projection(Credit.Defaultable(bond, "ACME"), ctx))
length(cfs), sum(cf.amount * discount(ctx, cf.time) for cf in cfs)
```

The quadrature nodes depend only on the times, so derivatives through a value keep the value
exactly, and a portfolio of defaultable and default-free contracts values as the sum of its
parts. Hazard, recovery and the discount curve are separate inputs: a spread added to the
curve discounts the cashflows without changing the default probabilities.

## Calibrating a hazard curve

[`fit`](@ref FinanceModels.fit) bootstraps an issuer's piecewise hazard curve from the full
prices of its bonds, with one segment ending at each maturity, holding the discount curve and the
recovery assumption fixed:

```@example fixed_income
curve = Yield.Constant(Continuous(0.04))
bonds = [Credit.Defaultable(Bond.Fixed(0.05, Periodic(2), t), "ACME") for t in (2.0, 5.0, 10.0)]
quotes = Quote.([1.005, 0.985, 0.95], bonds)
template = Credit.Assumptions(Credit.PiecewiseHazard(); recovery = 0.4)
fitted = fit(template, quotes, Fit.Bootstrap(); discount = curve)
fitted.hazard.rates
```

Each segment is solved, shortest first, so that its bond reprices. A price above the value without
default risk over a segment needs a negative intensity and throws, as do two bonds with one
maturity. First and second derivatives with respect to the prices, the bonds' coupons, the
recovery and the discount curve's parameters are exact, by the implicit function theorem at each
segment:

```@example fixed_income
using ForwardDiff
target = Credit.Defaultable(Bond.Fixed(0.06, Periodic(1), 7.0), "ACME")
value(p) = present_value(Models(curve, Dict("ACME" => fit(template, Quote.(p, bonds), Fit.Bootstrap(); discount = curve))), target)
ForwardDiff.gradient(value, [q.price for q in quotes])
```

## Callable and putable assets

[`Bond.Callable`](@ref FinanceModels.Bond.Callable) gives the issuer the right to redeem an asset
on payment dates, immediately after the payment, at a price per unit of the remaining principal;
[`Bond.Putable`](@ref FinanceModels.Bond.Putable) gives the holder the right to sell it back. Their
value depends on future interest rates, so they are valued on a Hull–White lattice:
[`ShortRate.HullWhiteTree`](@ref FinanceModels.ShortRate.HullWhiteTree) wraps a
`ShortRate.HullWhite` model as the valuation model. A plain curve has no volatility, so valuing a
callable asset on one throws rather than assume none.

```@example fixed_income
curve = ZeroRateCurve([0.03, 0.035, 0.04, 0.042], [1.0, 3.0, 7.0, 12.0])
hw = ShortRate.HullWhite(0.1, 0.01, curve)
bond = Bond.Fixed(0.05, Periodic(2), 10)
callable = Bond.Callable(bond, 3.0:0.5:9.5, 1.0)   # callable at par after three years
present_value(ShortRate.HullWhiteTree(hw), callable), present_value(curve, bond)
```

The lattice reprices the curve at every grid time, so an exercise right that is never worth using
leaves the asset's value unchanged to rounding, and a callable bond with one call date converges to
the bond less a receiver swaption, Jamshidian's closed form. Credit and prepayment combine with
exercise: a defaultable callable bond rolls back survival and recovery at each node, and a callable
pool is redeemed after its scheduled payment and prepayment.

```@example fixed_income
acme = Credit.Assumptions(Credit.ConstantHazard(0.02); recovery = 0.4)
ctx = Models(ShortRate.HullWhiteTree(hw), Dict("ACME" => acme))
present_value(ctx, Credit.Defaultable(callable, "ACME"))
```

The lattice is part of the valuation model, so a value, its sensitivities and an option-adjusted
spread (a spread added to the curve the lattice is calibrated to) all use the same lattice
settings. First derivatives follow the exercise decisions the lattice makes. Second derivatives
through an exercise decision throw: differentiating along fixed decisions leaves out the effect of
the exercise boundary moving, which is where a callable bond's negative convexity comes from, so
measure convexity by revaluing at bumped curves. Derivatives also converge more slowly than values
as the steps shrink, so use a smaller `max_step` for sensitivities.

[`exercise_policy`](@ref FinanceModels.exercise_policy) returns the lattice's decisions as a
function of the exercise time and the short rate observed then, for projections along simulated
rate paths:

```@example fixed_income
policy = exercise_policy(Models(ShortRate.HullWhiteTree(hw)), callable)
policy(5.0, 0.02), policy(5.0, 0.08)
```

## Asset-liability projections

[`AssetProjection`](@ref FinanceModels.AssetProjection) projects an asset's balance movements under
projection assumptions, which need not be the pricing ones: one row per payment date with the
opening and closing principal, interest, scheduled principal, prepayment, the principal outstanding
at default and its recovery, and the principal called or put and the cash paid for it. Every row
reconciles, and each closing balance is the next opening one:

```@example fixed_income
loan = Bond.LevelPayment(0.06, 12, 10)
pool = Credit.Defaultable(Bond.Pool(loan, "speed"), "ACME")
store = Dict(
    "speed" => Prepayment.CPR(0.08),
    "ACME" => Credit.Assumptions(Credit.ConstantHazard(0.01); recovery = 0.5),
)
rows = collect(Projection(pool, Models(curve, store), AssetProjection()))
first(rows)
```

```@example fixed_income
maximum(abs(r.opening - r.principal - r.prepayment - r.defaulted - r.redeemed - r.closing) for r in rows)
```

By default, defaults are expected: survival weights the payments, and each row reports the expected
default and recovery over the interval it ends. With `default_times` or an `rng`, a bond's issuer
defaults at one sampled time, shared by every bond with the same key (pools keep expected
defaults). Supplying the default times lets a portfolio carry any dependence between issuers.

An exercisable asset follows the decisions of its pricing lattice along a simulated rate path:
`AssetProjection(; pricing)` computes [`exercise_policy`](@ref FinanceModels.exercise_policy) once
per asset, and each exercise date reads only the short rate of the path then. The simulation's
time grid must include the exercise dates. Along paths of the pricing model, the average discounted
projection converges to the lattice price:

```@example fixed_income
using Random
paths = simulate(hw; n_scenarios = 2000, timestep = 1 / 12, horizon = 10.0, rng = Xoshiro(1))
kind = AssetProjection(; pricing = Models(ShortRate.HullWhiteTree(hw)))
cash(r) = r.interest + r.principal + r.recovery + r.redemption
pv_path(path) = sum(cash(r) * discount(path, r.time) for r in collect(Projection(callable, Models(path), kind)))
sum(pv_path, paths) / length(paths), present_value(ShortRate.HullWhiteTree(hw), callable)
```
