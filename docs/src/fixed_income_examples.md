# Worked Examples: a Callable Corporate Bond and a Mortgage Pool

Two complete examples of [Credit, Amortization, and Optionality](@ref): a ten-year callable
corporate bond valued on a fitted hazard curve and a Hull–White lattice, and a pool of level-payment
loans projected for asset-liability management. Option-adjusted spreads and rate sensitivities of
these assets are computed by ActuaryUtilities (`zspread`, `duration`, `convexity`), which accepts the
same valuation contexts.

## A callable corporate bond

The issuer has three straight bonds outstanding. Their full prices, the risk-free curve and a
recovery assumption determine its hazard curve, one segment per maturity:

```@example callable
using FinanceModels

curve = ZeroRateCurve([0.030, 0.034, 0.038, 0.041, 0.043], [1.0, 2.0, 5.0, 10.0, 20.0])
straight = [Credit.Defaultable(Bond.Fixed(c, Periodic(2), t), "ACME") for (c, t) in ((0.045, 2.0), (0.05, 5.0), (0.055, 10.0))]
prices = [1.0045, 1.0010, 0.9980]
acme = fit(Credit.Assumptions(Credit.PiecewiseHazard(); recovery = 0.4), Quote.(prices, straight), Fit.Bootstrap(); discount = curve)
acme.hazard.rates
```

The fitted assumptions reprice the bonds they were fitted to:

```@example callable
pricing = Models(curve, Dict("ACME" => acme))
[present_value(pricing, b) for b in straight] .- prices
```

The callable bond pays a 5.5% semiannual coupon for ten years and is callable at par on every
coupon date from year four. Its value needs interest-rate volatility, so the valuation model is a
Hull–White lattice calibrated to the same curve:

```@example callable
callable = Credit.Defaultable(Bond.Callable(Bond.Fixed(0.055, Periodic(2), 10), 4.0:0.5:9.5, 1.0), "ACME")
hw = ShortRate.HullWhite(0.05, 0.01, curve)
lattice = Models(ShortRate.HullWhiteTree(hw; max_step = 1 / 24), Dict("ACME" => acme))
value = present_value(lattice, callable)
noncall = present_value(pricing, Credit.Defaultable(Bond.Fixed(0.055, Periodic(2), 10), "ACME"))
value, noncall, noncall - value   # the call option's value to the issuer
```

The effective duration along the lattice's decisions is a first derivative, so ForwardDiff
computes it through the lattice, and through the hazard calibration too. Holding the fitted hazard
fixed measures rate risk at a fixed credit assumption, as ActuaryUtilities'
`duration(lattice, callable)` does. Refitting the hazard at each shifted curve, with the straight
bonds' prices held, lets the issuer's credit spread absorb the rate move on its own bonds, which
leaves little rate risk:

```@example callable
using ForwardDiff
tree_at(shifted) = ShortRate.HullWhiteTree(ShortRate.HullWhite(0.05, 0.01, shifted); max_step = 1 / 24)
fixed_credit(shift) = present_value(Models(tree_at(curve + Continuous(shift)), Dict("ACME" => acme)), callable)
function refit_credit(shift)
    shifted = curve + Continuous(shift)
    credit = fit(Credit.Assumptions(Credit.PiecewiseHazard(); recovery = 0.4), Quote.(prices, straight), Fit.Bootstrap(); discount = shifted)
    present_value(Models(tree_at(shifted), Dict("ACME" => credit)), callable)
end
-ForwardDiff.derivative(fixed_credit, 0.0) / value, -ForwardDiff.derivative(refit_credit, 0.0) / value
```

Second derivatives through the call decisions throw; measure convexity by revaluing at bumped
curves.

The decisions themselves are a function of the call date and the short rate then:

```@example callable
policy = exercise_policy(lattice, callable)
[policy(5.0, r) for r in (0.01, 0.03, 0.05, 0.07)]
```

## A mortgage pool

A pool of 30-year, 6% level-payment loans prepays at 8% CPR when priced; for the
asset-liability projection the company assumes 12% CPR and a 0.5% annual default intensity with
40% recovery:

```@example pool
using FinanceModels

loans = Bond.LevelPayment(0.06, 12, 30)
pool = Credit.Defaultable(Bond.Pool(loans, "speed"), "pool")
curve = Yield.Constant(0.05)
pricing = Models(curve, Dict("speed" => Prepayment.CPR(0.08), "pool" => Credit.Assumptions(Credit.ConstantHazard(0.005); recovery = 0.4)))
projection = Models(curve, Dict("speed" => Prepayment.CPR(0.12), "pool" => Credit.Assumptions(Credit.ConstantHazard(0.005); recovery = 0.4)))
present_value(pricing, pool)
```

The projection gives each month's balance movements per unit of initial principal:

```@example pool
rows = collect(Projection(pool, projection, AssetProjection()))
rows[1:3]
```

Every row reconciles, the balances chain from month to month, and the principal that leaves the
pool, scheduled, prepaid, defaulted or redeemed, adds up to its initial balance:

```@example pool
reconciled = maximum(abs(r.opening - r.principal - r.prepayment - r.defaulted - r.redeemed - r.closing) for r in rows)
chained = all(rows[i].closing == rows[i + 1].opening for i in 1:(length(rows) - 1))
runoff = sum(r.principal + r.prepayment + r.defaulted + r.redeemed for r in rows)
reconciled, chained, runoff
```

Annual cashflows for an asset-liability model are sums of the monthly rows:

```@example pool
year(r) = ceil(Int, r.time - 1e-9)
cash(r) = r.interest + r.principal + r.prepayment + r.recovery
annual = [sum(cash(r) for r in rows if year(r) == y) for y in 1:30]
annual[1:5]
```

At the pricing assumptions, the discounted projected cash differs from the pool's price only by
when recoveries are paid: the projection reports each month's recovery at the month's end, the
valuation when defaults occur.

```@example pool
priced = collect(Projection(pool, pricing, AssetProjection()))
sum(cash(r) * discount(curve, r.time) for r in priced), present_value(pricing, pool)
```
