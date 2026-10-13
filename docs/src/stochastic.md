# Stochastic Interest Rate Models

## Introduction

FinanceModels.jl includes stochastic short-rate models that are first-class yield models. Because Vasicek, Cox-Ingersoll-Ross (CIR), and Hull-White all have closed-form zero-coupon bond prices, they implement `discount(model, t)` analytically. This means the entire existing valuation infrastructure -- `zero`, `forward`, `par`, `present_value`, and `fit` -- works unchanged with these models.

For stochastic-cashflow analysis (e.g. Monte Carlo valuation), `simulate()` generates scenario yield curves that also plug into the existing `present_value`.

!!! note "Relationship to EconomicScenarioGenerators.jl"
    Short-rate and equity scenario generators were previously provided by [EconomicScenarioGenerators.jl](https://github.com/JuliaActuary/EconomicScenarioGenerators.jl). The `FinanceModels.ShortRate` models and `Equity.BlackScholesMerton`, together with `simulate` and `pv_mc`, are the maintained successors of its interest-rate and equity generators. Its copula-correlated, multi-model simulation (the `Correlated` generator) remains in EconomicScenarioGenerators.jl.

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

The CIR model has volatility proportional to `√r`. With positive parameters its continuous-time rate is nonnegative; the Feller condition `2ab ≥ σ²` makes zero inaccessible from a positive initial rate. The discretized simulation uses full truncation, described below.

```julia
cir = ShortRate.CoxIngersollRoss(0.3, 0.05, 0.1, Continuous(0.03))
```

### Hull-White

The Hull-White model takes an existing yield curve and adds stochastic dynamics. The drift is calibrated so that the model exactly reproduces the initial term structure.

```julia
curve = fit(Spline.Linear(), CMTYield.([0.04, 0.05, 0.055, 0.06], [1, 5, 10, 30]))
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
discount(v, 2, 10)    # P(0, 10) / P(0, 2): interval factor on the initial curve

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

This uses the future state `r(t)`. In contrast, `discount(v, T)` gives the time-zero price `P(0, T)`, and `discount(v, t, T)` gives `P(0, T) / P(0, t)`. The latter is an interval factor on the initial curve, not a conditional bond price given a simulated future rate.

The conditional form is used internally for derivative pricing (swaption
Jamshidian decomposition) and can be useful for scenario analysis where
the short rate at a future time is known.

### The Simulated Short Rate on a Path

Each simulated path records the short rate at its grid times. `short_rate` reads it as a
`Continuous` rate, and it is the state to pass to the conditional discount factor:

```julia
scenarios = simulate(v; n_scenarios=10, timestep=1/12, horizon=10.0)
path = scenarios[1]
r5 = short_rate(path, 5.0)    # the simulated r(5) on this path
discount(v, 5.0, 10.0, r5)    # P(5, 10 | r(5)) on this path
```

`short_rate` answers at the path's grid times, allowing floating-point rounding around a grid time. `FinanceModels.simulation_times(path)` returns a copy of that grid. Other times throw a `DomainError`, because the simulation has no recorded rate there. A hand-built `RatePath(interp)` has no recorded rates, so `short_rate` throws an `ArgumentError` on it. `Yield.instantaneous_forward(path, t)` is the slope of the path's discount integral over the
step that starts at `t`, the average of that step's two rates. It uses the rate at the end of the
step, so it is not the state at `t`.

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

# Compare fitted and observed prices; a parametric fit need not reprice exactly
[(fitted=present_value(v_fitted, q.instrument), observed=q.price) for q in quotes]
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

`simulate` generates interest-rate paths using the exact Gaussian transition density for Vasicek and Hull-White (no discretisation bias in the short rate at any timestep) and the full truncation scheme (Lord, Koekkoek & Van Dijk, 2010) for Cox-Ingersoll-Ross. Each path is returned as a [`RatePath`](@ref FinanceModels.RatePath), which is itself an `AbstractYieldModel` -- so `discount`, `zero`, `forward`, `par`, and `present_value` all work on individual scenarios. Path discount factors use trapezoidal integration of the short rate, so they still have time-discretization error even when rate transitions are exact. A path is defined only through its last grid time: evaluating it later throws. An aligned horizon is included exactly; otherwise the final full step ends just beyond the requested horizon. Choose a `horizon` that covers every cashflow you value and inspect `FinanceModels.simulation_times(path)` when grid alignment matters.

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

# Differences reflect both sampling error and time-discretization error.
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

For floating coupons fixed in advance, use the **conditional term structure at each reset**. A `RatePath` stores realized discount integrals: `forward(path, reset, payment)` uses rates realized after the reset and therefore cannot serve as a fixing known at that reset. Passing a path as the index curve to `Projection` does not perform conditional valuation automatically.

Here is a single-curve Vasicek example with a 2% margin, semiannual coupons, and unit principal. At each reset, convert the conditional bond price into the next period's simple interest accrual:

```@example floating_scenarios
using FinanceModels, Random

v = ShortRate.Vasicek(0.136, 0.0168, 0.0119, Continuous(0.01))
scenarios = simulate(v; n_scenarios=3, timestep=1/12, horizon=3.0,
                     rng=MersenneTwister(42))

function floating_cashflows(model, path; maturity=3.0, accrual=0.5, margin=0.02)
    [begin
         reset = t - accrual
         conditional_df = discount(model, reset, t, short_rate(path, reset))
         coupon = inv(conditional_df) - 1 + margin * accrual
         Cashflow(coupon + (t == maturity ? 1 : 0), t)
     end for t in accrual:accrual:maturity]
end

cashflows = [floating_cashflows(v, path) for path in scenarios]
[present_value(path, cfs) for (path, cfs) in zip(scenarios, cashflows)]
```

The reset dates must lie on the simulation grid. The example assumes full coupon periods and one curve for projection and discounting. Use a valuation callback with the appropriate fixing logic for more complex contracts; `pv_mc` does not supply that logic for a floating contract automatically.

## Equity Price Paths

`simulate` also takes an `Equity.BlackScholesMerton` model and returns [`PricePath`](@ref FinanceModels.PricePath)s on the
same grid as rate paths:

```julia
m = Equity.BlackScholesMerton(0.03, 0.01, 0.2)   # r, q, σ
paths = simulate(m; n_scenarios=1000, timestep=1/12, horizon=10.0)
p = paths[1]
p(5.0)                                            # the price at year 5, relative to the start
100 .* p.(FinanceModels.simulation_times(p))      # prices in currency from a spot of 100
```

Each step is the exact lognormal transition, so the prices at grid times have no discretisation
bias. A path starts at 1 and is the ex-dividend price relative to its start; one share with its
dividends reinvested is worth `exp(q * t) * p(t)`. Like [`short_rate`](@ref), a path answers only
at its grid times: the simulation has no price between them, so a payoff at a time needs a grid
that contains it.

The drift is the model's constant `r - q`: risk-neutral when `r` is the risk-free rate, real-world
when `r` is an expected return. It is never a simulated short rate. `σ` is a number or a
`Volatility.Constant`; `Volatility.Constant()` is zero volatility, so its paths are deterministic.

## Mean Reversion Behaviour

The expected short rate in Vasicek and CIR tends toward `b`. The long-maturity zero rate is different because bond pricing averages discounted outcomes. For Vasicek with `a > 0`, its limit is `b - σ² / (2a²)`, close to `b` when mean reversion is strong relative to volatility:

```julia
# Strong mean reversion: long-maturity zero-rate limit is 0.0499875
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

- **CIR Feller condition:** For the continuous-time CIR model started above zero, `2ab ≥ σ²` makes zero inaccessible. When this condition fails, the rate can reach zero. The simulation uses the full truncation scheme: an auxiliary process may go negative while the observed rate is floored at zero. Check convergence with a finer `timestep`, especially when the Feller condition is strongly violated. ForwardDiff derivatives through clipped states remain finite when the clipped state has zero derivatives; the clipping boundary itself is a kink.
- **Vasicek negative rates:** The Vasicek model allows negative rates by design. For very negative rates or long horizons, discount factors may exceed 1 (equivalently, zero rates are negative). This is consistent with observed negative-rate environments.

### Performance Guidance

- `simulate` allocates one `LinearInterpolation` and one vector of recorded short rates per scenario. Memory therefore grows with both path count and grid size; retaining the simulated rates adds one value per grid time per path.
- `pv_mc` is embarrassingly parallel across scenarios but is not parallelised internally. For large-scale simulations, users can parallelise with `Threads.@threads` or `Distributed`:

```julia
using Base.Threads

scenarios = simulate(v; n_scenarios=10_000, timestep=1/12, horizon=31.0)
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
| `short_rate(path, t)` | The simulated `r(t)` at a grid time of a `RatePath`, a `Continuous` rate |
| `simulate(Equity.BlackScholesMerton(...); ...)` | Generate `Vector{PricePath}` equity price paths |
| `p(t)` | The simulated price at a grid time of a `PricePath`, relative to its start |

`FinanceModels.simulation_times(path)` lists the grid times of a `RatePath` or `PricePath`; `FinanceModels.simulation_steps(horizon, timestep)` reports the step count and whether the requested horizon is aligned.
