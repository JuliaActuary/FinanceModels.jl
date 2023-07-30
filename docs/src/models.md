# Models, Valuation, Projections, and Fitting

## Page Contents

```@contents
Pages = ["models.md"]
Depth = 4
```

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

- `Yields.Periodic(m)` for rates that compound `m` times per period (e.g. `m` times per year if working with annual rates).
- `Yields.Continuous()` for continuously compounding rates.

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

```julia-repl
r = Rate(Yields.Periodic(12),0.01)             # rate that compounds 12 times per rate period (ie monthly)

convert(Yields.Periodic(1),r)                  # convert monthly rate to annual effective
convert(Yields.Continuous(),r)          # convert monthly rate to continuous
```

#### Arithmetic

Adding, subtracting, multiplying, dividing, and comparing rates is supported.

### Available Models - Yields

- [`Yield.Constant`](@ref)
- Bootstrapped [`Spline`](@ref)
- [`Yield.SmithWilson`](@ref)
- [`Yield.NelsonSiegel`](@ref)
- [`Yield.NelsonSiegelSvensson`](@ref)

Yield models can also be composed. Here is an example of fitting rates and spread separately and then adding the two models together:

```julia-repl
julia> q_rate = ZCBYield([0.01,0.02,0.03]);

julia> q_spread = ZCBYield([0.01,0.01,0.01]);

julia> m_rate = fit(Spline.Linear(),q_rate,Fit.Bootstrap());⠀           

julia> m_spread = fit(Spline.Linear(),q_spread,Fit.Bootstrap());

julia> forward(m_spread + m_rate,0,1)
Rate{Float64, Continuous}(0.01980262729617973, Continuous())

julia> forward(m_spread + m_rate,0,1) |> Periodic(1)
Rate{Float64, Periodic}(0.020000000000000018, Periodic(1))

julia> discount(m_spread + m_rate,0,3)
0.8889963586709149

julia> discount(0.04,3)
0.8889963586709148
```

### Creating New Yield Models

See the [Guide](guide.html) for an example of creating a model from scratch. Some additional aspects to note:

- The only method that must be defined to calculate the [`present_value`](@ref) of something is [`discount`](@ref). Other methods will be inferred.
- Other methods that are imputed by default, but can be extended include: [`accumulation`](@ref), [`forward`](@ref), [`par`](@ref), [`zero`](@ref), and [`rate`](@ref).

## Equity and Volatility Models

### Available Models - Option Valuation

- [`Option.BlackScholesMerton`](@ref)

### Available Models - Volatility

- [`Volatility.Constant`](@ref)

#### Creating new Volatility Models

A volatility model must extend `volatility(vol::Volatility.MyNewModel, strike_ratio, time_to_maturity)`.
