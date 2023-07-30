# Models, Valuation, and Fitting

## Page Contents

```@contents
Pages = ["contracts.md"]
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

Adding, substracting, multiplying, dividing, and comparing rates is supported.

## Available Models

TODO

