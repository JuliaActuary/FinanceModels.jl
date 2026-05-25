# Models, Valuation, Projections, and Fitting

## Introduction

Conceptually, we have an iterative process:

1. We use models to value contracts
2. We use observed (or assumed) prices to calibrate models

Thus the discussion of model calibration and valuation of contracts is inextricably linked together.

## Yield (Interest Rate) models

### Rates

We should first discuss `Rate`s, which are reexported from [`FinanceCore.jl`](https://github.com/JuliaActuary/FinanceCore.jl)

Rates are types that wrap scalar values to provide information about how to determine `discount` and `accumulation` factors. These allow for explicit handling of rate compounding conventions which, if not explicit, is often a source of errors in practice.

There are two `Frequency` types:

- `Periodic(m)` for rates that compound `m` times per period (e.g. `m` times per year if working with annual rates).
- `Continuous()` for continuously compounding rates.

#### Examples

```julia
Continuous(0.05)       # 5% continuously compounded
Periodic(0.05,2)       # 5% compounded twice per period
```

These are both subtypes of the parent `Rate` type and are instantiated as:

```julia
Rate(0.05,Continuous())       # 5% continuously compounded
Rate(0.05,Periodic(2))        # 5% compounded twice per period
```

Broadcast over a vector to create `Rates` with the given compounding:

```julia
Periodic.([0.02,0.03,0.04],2) 
Continuous.([0.02,0.03,0.04]) 
```

Rates can also be constructed by specifying the `CompoundingFrequency` and then passing a scalar rate:

```julia
Periodic(1)(0.05)
Continuous()(0.05)
```

#### Conversion

Convert rates between different types with `convert`. E.g.:

```julia
r = Rate(0.01,Periodic(12))             # rate that compounds 12 times per rate period (ie monthly)

convert(Periodic(1),r)                  # convert monthly rate to annual effective
convert(Continuous(),r)          # convert monthly rate to continuous
```

To get the scalar value out of the `Rate`, use `FinanceModels.rate(r)`:

```julia-rel
julia> r = Rate(0.01,Periodic(12));   
julia> rate(r)
0.01

```

### Available Models - Yields

- [`FinanceModels.Yield.Constant`](@ref)
- Bootstrapped [`Spline`](@ref FinanceModels.Spline)s
- [`FinanceModels.Yield.SmithWilson`](@ref)
- [`FinanceModels.Yield.NelsonSiegel`](@ref)
- [`FinanceModels.Yield.NelsonSiegelSvensson`](@ref)
- [`FinanceModels.Yield.CairnsPritchard`](@ref)
- [`FinanceModels.Yield.CairnsPritchardExtended`](@ref)

#### Curve Transformations

- `curve + curve` — additive composition of zero rates ([`FinanceModels.Yield.CompositeYield`](@ref))
- `curve * scalar` / `curve / scalar` — scale zero rates ([`FinanceModels.Yield.ScaledYield`](@ref))
- [`FinanceModels.Yield.TenorShift`](@ref) — lazy zero-rate shift depending on tenor: `(rate, tenor) -> Rate` (formerly `TransformedYield`, retained as alias)
- [`FinanceModels.Yield.ProjectedShift`](@ref) — lazy zero-rate shift depending on tenor *and* a projection time `τ`: `(τ, rate, tenor) -> Rate`
- [`FinanceModels.Yield.ForwardStarting`](@ref) — rebase a curve to a new time-zero

### Available Models - Stochastic Short Rates

Stochastic short-rate models with closed-form zero-coupon bond prices. These are full yield models that also support Monte Carlo simulation. See the [Stochastic Models](@ref) guide for details and examples.

- [`FinanceModels.ShortRate.Vasicek`](@ref)
- [`FinanceModels.ShortRate.CoxIngersollRoss`](@ref)
- [`FinanceModels.ShortRate.HullWhite`](@ref)

#### Arithmetic

Adding, subtracting, multiplying, dividing, and comparing rates is supported.

Yield models can also be composed. Here is an example of fitting rates and spread separately and then adding the two models together:

```julia-repl
julia> q_rate = ZCBYield([0.01,0.02,0.03]);

julia> q_spread = ZCBYield([0.01,0.01,0.01]);

julia> model_rate = fit(Spline.Linear(),q_rate,Fit.Bootstrap());⠀

julia> model_spread = fit(Spline.Linear(),q_spread,Fit.Bootstrap());

julia> forward(model_spread + model_rate,0,1)
Rate{Float64, Continuous}(0.01980262729617973, Continuous())

julia> forward(model_spread + model_rate,0,1) |> Periodic(1)
Rate{Float64, Periodic}(0.020000000000000018, Periodic(1))

julia> discount(model_spread + model_rate,0,3)
0.8889963586709149

julia> discount(0.04,3)
0.8889963586709148
```

!!! warning "Caution with Spreads"

    It is fairly common to see spreads and rates provided separately where both are quoted in par convention. For example, US Treasury par rates with the associated par risk spreads. Because par rates are dependent on the amount and path of rates preceeding the given tenor, **it is not valid to construct a "spread curve" with par rates and then use it in composition with a "rate curve"**.
    
    That is, while the zero rates and spreads in the preceeding example allow for additive or subtractive composition, it is not the case for par rates and spreads. Note the different discount factors produced:

    ```julia
    q_rate = ParYield([0.01,0.02,0.03]);
    q_spread = ParYield([0.01,0.01,0.01]);
    q_yield = ParYield([0.02,0.03,0.04]);

    model_rate = fit(Spline.Linear(),q_rate,Fit.Bootstrap());
    model_spread = fit(Spline.Linear(),q_spread,Fit.Bootstrap());
    model_yield = fit(Spline.Linear(),q_yield,Fit.Bootstrap());

    # The curves are different!
    discount(model_spread + model_rate,3)
    # 0.8889963586709149

    discount(model_yield,3)
    # 0.8864366955434709
    ```

#### Yield Shifts: `TenorShift` and `ProjectedShift`

Two concrete subtypes of [`Yield.AbstractYieldShift`](@ref FinanceModels.Yield.AbstractYieldShift) apply lazy transformations to a base curve's zero rates without discretizing or refitting. They differ in whether the shift depends on one or two time axes.

##### `TenorShift` — shift depends on tenor only

[`Yield.TenorShift`](@ref FinanceModels.Yield.TenorShift) (formerly `TransformedYield`, retained as a deprecated alias) applies a shift that may vary with the tenor `t`. The rule function receives the base curve's `Continuous` zero rate and the tenor, and returns a new rate:

```julia-repl
julia> base = Yield.Constant(0.05);

julia> # Parallel shift (+100 bp) using Rate arithmetic
       shifted = base + (z, t) -> z + Periodic(0.01, 1);

julia> zero(shifted, 10)  # ≈ 0.05 + log(1.01), not 0.06 — Rate arithmetic converts Periodic → Continuous
Rate{Float64, Continuous}(0.05995033085316808, Continuous())

julia> # Continuous shift (simpler when convention is known)
       shifted2 = base + (z, t) -> z + Continuous(0.01);

julia> zero(shifted2, 10)
Rate{Float64, Continuous}(0.06, Continuous())

julia> # Tenor-dependent twist (steepener that fades at 30y)
       twist = base + (z, t) -> z + Continuous(0.02 * max(0.0, 1.0 - t/30.0));

julia> zero(twist, 1)   # ≈ 5% + 1.93%
Rate{Float64, Continuous}(0.06933333333333334, Continuous())

julia> zero(twist, 30)  # shift is zero at 30y
Rate{Float64, Continuous}(0.05, Continuous())
```

The rule function has the signature `(z::Rate, t) -> Rate`. The return value is type-asserted as `Rate`, so rules must carry compounding convention explicitly — returning a plain `Real` raises a `TypeError`. Because the rule receives the zero rate itself, `Rate` arithmetic like `z + Periodic(0.01, 1)` handles compounding conversion correctly — no manual convention juggling needed.

The `+` operator provides ergonomic construction: `curve + f` or `f + curve` both create a `TenorShift`. This is distinct from `curve + scalar` (which creates a `CompositeYield`).

##### `ProjectedShift` — shift depends on tenor *and* a projection time

[`Yield.ProjectedShift`](@ref FinanceModels.Yield.ProjectedShift) extends `TenorShift` with a second time axis `τ`: the **projection time** (as-of / valuation-date offset) at which the curve is being evaluated. The rule signature is `(τ, z::Rate, t) -> Rate`, separating the projection time from the tenor.

This is the natural shape for scenario frameworks where a shift's magnitude evolves across a multi-year horizon — BMA SBA phase-ins, IFRS17 macroeconomic scenarios, embedded-value runoffs:

```julia-repl
julia> base = Yield.Constant(0.05);

julia> # -150 bp parallel shift, phased in linearly over 10 projection years.
       phase_in = (τ, z, _) -> z + Continuous(-0.015 * min(τ, 10) / 10);

julia> # Curve as seen at projection year 3 (30% phased in → -45 bp).
       c3 = Yield.ProjectedShift(base, phase_in, 3.0);

julia> zero(c3, 5)
Rate{Float64, Continuous}(0.04550000000000001, Continuous())

julia> # Curve as seen at projection year 10 (fully phased in → -150 bp).
       c10 = Yield.ProjectedShift(base, phase_in, 10.0);

julia> zero(c10, 5)
Rate{Float64, Continuous}(0.035, Continuous())
```

The intended pattern: store `phase_in` once as a first-class, year-independent value, then call `ProjectedShift(base, phase_in, τ)` at each projection time `τ` in a projection loop.

There is no `+` operator sugar for `ProjectedShift` — fixing `τ` at composition time defeats the purpose of storing the rule as a year-independent value. Always use the explicit constructor.

### Creating New Yield Models

See the [FinanceModels.jl Guide](@ref) for an example of creating a model from scratch. Some additional aspects to note:

- The only method that must be defined to calculate the [`FinanceCore.present_value`](@ref) of something is [`FinanceCore.discount`](@ref). Other methods will be inferred.
- Other methods that are imputed by default, but can be extended include: [`FinanceCore.accumulation`](@ref), [`FinanceModels.forward`](@ref), [`FinanceModels.par`](@ref), [`FinanceModels.zero`](@ref), and [`FinanceModels.rate`](@ref).

## Equity and Volatility Models

### Available Models - Option Valuation

- [`FinanceModels.Equity.BlackScholesMerton`](@ref)

### Available Models - Volatility

- [`FinanceModels.Volatility.Constant`](@ref)

#### Creating new Volatility Models

A volatility model must extend `volatility(vol::Volatility.MyNewModel, strike_ratio, time_to_maturity)`.
