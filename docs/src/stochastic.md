# Stochastic Interest Rate Models

## Introduction

FinanceModels.jl includes stochastic short-rate models that are first-class yield models. Because Vasicek, Cox-Ingersoll-Ross (CIR), and Hull-White all have closed-form zero-coupon bond prices, they implement `discount(model, t)` analytically. This means the entire existing valuation infrastructure -- `zero`, `forward`, `par`, `present_value`, and `fit` -- works unchanged with these models.

For stochastic-cashflow analysis (e.g. Monte Carlo valuation), `simulate()` generates scenario yield curves that also plug into the existing `present_value`.

## Available Models

| Model | Dynamics | Parameters |
|:------|:---------|:-----------|
| [`ShortRate.Vasicek`](@ref FinanceModels.ShortRate.Vasicek) | `dr = a(b - r)dt + σ dW` | `a`, `b`, `σ`, `initial` |
| [`ShortRate.CoxIngersollRoss`](@ref FinanceModels.ShortRate.CoxIngersollRoss) | `dr = a(b - r)dt + σ√r dW` | `a`, `b`, `σ`, `initial` |
| [`ShortRate.HullWhite`](@ref FinanceModels.ShortRate.HullWhite) | `dr = (θ(t) - ar)dt + σ dW` | `a`, `σ`, `curve` |

Where:
- `a` is the speed of mean reversion
- `b` is the long-term mean rate
- `σ` is the volatility
- `initial` is the initial short rate `r₀` (a `Rate` or scalar)
- `curve` is an existing yield model (for Hull-White, which calibrates to an initial term structure)

## Constructing Models

### Vasicek

The Vasicek model is the simplest mean-reverting short-rate model. The short rate `r(t)` reverts to a long-term level `b` at speed `a`, with constant volatility `σ`.

```julia
using FinanceModels

v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))
```

The `initial` rate can be passed as a scalar (interpreted as continuous) or as an explicit `Rate`:

```julia
# These are equivalent:
ShortRate.Vasicek(0.136, 0.0168, 0.0119, 0.01)
ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))
```

### Cox-Ingersoll-Ross

The CIR model is similar to Vasicek but the volatility is proportional to `√r`, which prevents negative rates when the Feller condition `2ab > σ²` is satisfied.

```julia
cir = ShortRate.CoxIngersollRoss(0.3, 0.05, 0.1, Continuous(0.03))
```

### Hull-White

The Hull-White model takes an existing yield curve and adds stochastic dynamics. The drift is calibrated so that the model exactly reproduces the initial term structure.

```julia
curve = fit(Spline.Cubic(), CMTYield.([0.04, 0.05, 0.055, 0.06], [1, 5, 10, 30]), Fit.Bootstrap())
hw = ShortRate.HullWhite(0.1, 0.01, curve)

# discount factors match the initial curve exactly:
discount(hw, 10) == discount(curve, 10) # true
```

## Using Stochastic Models as Yield Curves

Since these models implement `discount()`, all standard yield curve operations work:

```julia
v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))

# Discount factors (closed-form zero-coupon bond prices)
discount(v, 5)        # P(0, 5)
discount(v, 2, 10)    # P(2, 10)

# Zero rates
zero(v, 5)            # continuous zero rate at t=5

# Forward rates
forward(v, 2, 3)      # forward rate from t=2 to t=3

# Par yields
par(v, 10)            # par yield at t=10
```

### Conditional Discount Factors

Stochastic models also support the **conditional** zero-coupon bond price
`P(t, T | r(t) = r)` — the price at time `t` of a bond maturing at `T`,
given that the short rate at `t` is `r`:

```julia
v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))

# P(2, 10 | r(2) = 0.04) under Vasicek
discount(v, 2.0, 10.0, 0.04)
```

This is distinct from the 2-argument form `discount(v, 10)` which gives the
**unconditional** `P(0, T)` from the initial term structure.

The conditional form is used internally for derivative pricing (swaption
Jamshidian decomposition) and can be useful for scenario analysis where
the short rate at a future time is known.

### Extracting the Short Rate from a Simulated Path

After simulation, you can extract the instantaneous short rate `r(t)` from a
`RatePath` using `short_rate`:

```julia
scenarios = simulate(v; n_scenarios=10, timestep=1/12, horizon=10.0)
short_rate(scenarios[1], 5.0)  # r(5) for the first scenario
```

### Valuing Fixed-Income Contracts

Since stochastic models are yield models, `present_value` works directly for deterministic-cashflow instruments:

```julia
v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))
bond = Bond.Fixed(0.05, Periodic(2), 10)

# Analytical present value using closed-form discount factors
present_value(v, bond)
```

## Calibrating Models with `fit`

Stochastic models support calibration via `fit`, using the same API as all other models. The optimizer uses ForwardDiff automatic differentiation to find parameters that best match observed market quotes.

```julia
# Observed market zero-coupon yields
quotes = ZCBYield.([0.02, 0.025, 0.03], [1, 5, 10])

# Initial guess
v0 = ShortRate.Vasicek(0.1, 0.02, 0.01, Continuous(0.01))

# Fit to market data
v_fitted = fit(v0, quotes)

# Verify: the fitted model reprices the quotes
map(q -> present_value(v_fitted, q.instrument), quotes)
```

The default parameter bounds for fitting are:

| Model | `a` | `b` | `σ` | `initial` |
|:------|:----|:----|:----|:----------|
| Vasicek | `0.0 .. 5.0` | `-0.1 .. 0.5` | `0.0 .. 1.0` | `-0.05 .. 0.2` |
| CIR | `0.0 .. 5.0` | `0.0 .. 0.5` | `0.0 .. 1.0` | `0.0 .. 0.2` |
| Hull-White | `0.0 .. 5.0` | -- | `0.0 .. 1.0` | -- |

Custom bounds can be passed via the `variables` keyword argument to `fit`.

## Monte Carlo Simulation

### Generating Scenarios with `simulate`

`simulate` uses Euler-Maruyama discretisation to generate interest-rate paths. Each path is returned as a [`RatePath`](@ref FinanceModels.RatePath), which is itself an `AbstractYieldModel` -- so `discount`, `zero`, `forward`, `par`, and `present_value` all work on individual scenarios.

```julia
using Random

v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))

scenarios = simulate(v;
    n_scenarios = 1000,    # number of paths
    timestep    = 1/12,    # monthly steps
    horizon     = 30.0,    # 30-year horizon
    rng         = MersenneTwister(42),  # reproducible
)

length(scenarios)            # 1000
scenarios[1] isa RatePath    # true

# Each scenario is a full yield model:
discount(scenarios[1], 5)
zero(scenarios[1], 10)
present_value(scenarios[1], Bond.Fixed(0.05, Periodic(2), 10))
```

### Monte Carlo Present Value with `pv_mc`

`pv_mc` is a convenience function that simulates scenarios and averages `present_value` across them:

```julia
v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))
bond = Bond.Fixed(0.05, Periodic(2), 10)

# Monte Carlo expected PV
mc = pv_mc(v, bond; n_scenarios=5000, timestep=1/12)

# Compare to analytical (closed-form) PV
analytical = present_value(v, bond)

# These should be close (within ~1-2% for 5000 scenarios)
```

The signature is:

```julia
pv_mc(model, contract;
    n_scenarios = 1000,
    timestep    = 1/12,
    horizon     = nothing,   # defaults to maturity + 1
    rng         = Random.default_rng(),
)
```

### Working with Individual Scenarios

Since each scenario is a yield model, you can do per-scenario analysis:

```julia
v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))
bond = Bond.Fixed(0.05, Periodic(2), 10)

scenarios = simulate(v; n_scenarios=1000, timestep=1/12, horizon=11.0)

# Distribution of present values
pvs = [present_value(sc, bond) for sc in scenarios]

# Percentiles
sort!(pvs)
p95 = pvs[950]   # 95th percentile PV
p05 = pvs[50]    # 5th percentile PV
mean_pv = sum(pvs) / length(pvs)
```

### Projecting Cashflows Across Scenarios

For fixed-coupon bonds, the cashflows themselves don't change across scenarios -- only the discount factors (and thus the present value) change. You can get the cashflows directly with `collect`:

```julia
bond = Bond.Fixed(0.05, Periodic(2), 3)
collect(bond)
# 6-element Vector{Cashflow}:
#  Cashflow(0.025, 0.5)
#  Cashflow(0.025, 1.0)
#  ...
#  Cashflow(1.025, 3.0)
```

For **floating-rate bonds**, the cashflows depend on forward rates, which differ across scenarios. Use `Projection` with a `Dict` mapping the reference rate key to the scenario's yield model:

```julia
using FinanceModels, Random

v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))
scenarios = simulate(v; n_scenarios=3, timestep=1/12, horizon=4.0,
                     rng=MersenneTwister(42))

# Floating bond: 2% spread over "SOFR", semiannual, 3-year
bond_float = Bond.Floating(0.02, Periodic(2), 3.0, "SOFR")

for (i, sc) in enumerate(scenarios)
    proj = Projection(bond_float, Dict("SOFR" => sc), CashflowProjection())
    cfs = collect(proj)
    println("Scenario $i cashflows:")
    for cf in cfs
        println("  t=$(cf.time): $(round(cf.amount; digits=4))")
    end
end
# Each scenario produces different coupon amounts because the
# forward rates from each RatePath are different.
```

This also works with composite contracts like interest rate swaps, where the floating leg references a scenario:

```julia
curve = Yield.Constant(0.05)
swap = InterestRateSwap(curve, 5)

sc = scenarios[1]
proj = Projection(swap, Dict("OIS" => sc), CashflowProjection())
cfs = collect(proj)  # net cashflows (received fixed, paid floating)
```

## Mean Reversion Behaviour

A key feature of these models is mean reversion. With strong mean reversion (`a` large), the long-term zero rate converges to `b`:

```julia
# Strong mean reversion: long-term rate approaches b=0.05
v_strong = ShortRate.Vasicek(2.0, 0.05, 0.01, Continuous(0.10))
zero(v_strong, 30)  # close to Continuous(0.05)

# Weak mean reversion: initial rate persists longer
v_weak = ShortRate.Vasicek(0.01, 0.05, 0.01, Continuous(0.10))
zero(v_weak, 30)    # still far from 0.05
```

## Gaussian Model Derivative Pricing

Both the Vasicek and Hull-White models are Gaussian short-rate models and share the same closed-form ZCB option formula (Black's formula). This means all derivative pricing — ZCB options, caps, floors, and swaptions — works identically for both models.

For Hull-White, since `discount(hw, t)` simply returns the initial curve's discount factor (the model is calibrated to match the curve exactly), bond prices alone cannot identify `a` and `σ`. Instead, Hull-White is typically calibrated to **derivative prices** — caps, floors, swaptions, or zero-coupon bond options.

### Zero-Coupon Bond Options

```julia
curve = Yield.Constant(Continuous(0.05))
hw = ShortRate.HullWhite(0.1, 0.015, curve)

# Call on a ZCB: right to buy at time T=1 a ZCB maturing at S=5 for strike K
call = present_value(hw, Option.ZCBCall(1.0, 5.0, 0.75))
put  = present_value(hw, Option.ZCBPut(1.0, 5.0, 0.75))
```

### Caps and Floors

A cap is a portfolio of caplets, each paying `max(L - K, 0) · τ` where `L` is the simply-compounded forward rate. Under Hull-White, each caplet is equivalent to a scaled put on a zero-coupon bond.

```julia
hw = ShortRate.HullWhite(0.03, 0.02, Yield.Constant(Continuous(0.01)))

# 3% strike, quarterly resets, 2-year maturity
cap = present_value(hw, Option.Cap(0.03, 4, 2.0))
flr = present_value(hw, Option.Floor(0.03, 4, 2.0))
```

Cap-floor parity holds: `Cap(K) - Floor(K) = forward swap value`.

### Swaptions

A European swaption gives the right to enter a swap at expiry. Pricing uses the Jamshidian (1989) decomposition into zero-coupon bond options.

```julia
hw = ShortRate.HullWhite(0.03, 0.02, Yield.Constant(Continuous(0.01)))

# 1y into 4y payer swaption, 1.1% strike, quarterly
payer = present_value(hw, Option.Swaption(1.0, 5.0, 0.011, 4; payer=true))

# Receiver swaption
receiver = present_value(hw, Option.Swaption(1.0, 5.0, 0.011, 4; payer=false))
```

### Calibrating Hull-White to Derivatives

With derivative pricing, `fit` can calibrate Hull-White's `a` and `σ` to market swaption or cap prices:

```julia
curve = Yield.Constant(Continuous(0.03))

# Initial guess
hw0 = ShortRate.HullWhite(0.05, 0.01, curve)

# Market swaption prices (here generated from a "true" model)
hw_true = ShortRate.HullWhite(0.1, 0.015, curve)
instruments = [
    Option.Swaption(1.0, 6.0, 0.03, 2),
    Option.Swaption(2.0, 7.0, 0.03, 2),
    Option.Swaption(3.0, 8.0, 0.03, 2),
]
quotes = [Quote(present_value(hw_true, inst), inst) for inst in instruments]

# Calibrate
hw_fit = fit(hw0, quotes)
```

## Important Notes

### Model Behaviour

- **CIR Feller condition:** The CIR model requires `2ab > σ²` for the short rate to stay strictly positive. When this condition is violated, the rate can reach zero. The simulation uses the full truncation scheme (absorption at zero) in that case.
- **Vasicek negative rates:** The Vasicek model allows negative rates by design. For very negative rates or long horizons, discount factors may exceed 1 (equivalently, zero rates are negative). This is consistent with observed negative-rate environments.

### Performance Guidance

- `simulate` allocates one `LinearInterpolation` per scenario. For 100k scenarios with a 30-year monthly horizon, expect execution times on the order of seconds.
- `pv_mc` is embarrassingly parallel across scenarios but is not parallelised internally. For large-scale simulations, users can parallelise with `Threads.@threads` or `Distributed`:

```julia
using Base.Threads

scenarios = simulate(v; n_scenarios=100_000, timestep=1/12, horizon=31.0)
pvs = Vector{Float64}(undef, length(scenarios))
@threads for i in eachindex(scenarios)
    pvs[i] = present_value(scenarios[i], bond)
end
mean_pv = sum(pvs) / length(pvs)
```

## Summary

| Function | Description |
|:---------|:------------|
| `discount(model, t)` | Closed-form ZCB price `P(0,t)` |
| `discount(model, t, T, r_t)` | Conditional ZCB price `P(t,T \| r(t)=r)` |
| `zero(model, t)` | Continuous zero rate at `t` |
| `forward(model, t1, t2)` | Forward rate from `t1` to `t2` |
| `par(model, t)` | Par yield at `t` |
| `present_value(model, contract)` | Analytical present value |
| `present_value(m, Option.ZCBCall(...))` | Vasicek/Hull-White ZCB option price |
| `present_value(m, Option.Cap(...))` | Vasicek/Hull-White cap price |
| `present_value(m, Option.Swaption(...))` | Vasicek/Hull-White swaption price |
| `fit(model, quotes)` | Calibrate to market data |
| `simulate(model; ...)` | Generate `Vector{RatePath}` scenarios |
| `pv_mc(model, contract; ...)` | Monte Carlo expected present value |
| `short_rate(path, t)` | Extract `r(t)` from a simulated `RatePath` |
