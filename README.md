# Yields

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaActuary.github.io/Yields.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaActuary.github.io/Yields.jl/dev)
[![Build Status](https://github.com/JuliaActuary/Yields.jl/workflows/CI/badge.svg)](https://github.com/JuliaActuary/Yields.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaActuary/Yields.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaActuary/Yields.jl)
[![lifecycle](https://img.shields.io/badge/LifeCycle-Developing-blue)](https://www.tidyverse.org/lifecycle/)


**Yields** provides a simple interface for constructing, manipulating, and using yield curves for modeling purposes.

It's intended to provide common functionality around modeling interest rates, spreads, and miscellaneous yields across the JuliaActuary ecosystem (though not limited to use in JuliaActuary packages).

## QuickStart

```julia
using Yields

riskfree_maturities = [0.5, 1.0, 1.5, 2.0]
riskfree    = [5.0, 5.8, 6.4, 6.8] ./ 100     #spot rates, annual effective if unspecfied

spread_maturities = [0.5, 1.0, 1.5, 3.0]      # different maturities
spread    = [1.0, 1.8, 1.4, 1.8] ./ 100       # spot spreads

rf_curve = Yields.Zero(riskfree,riskfree_maturities)
spread_curve = Yields.Zero(spread,spread_maturities)


yield = rf_curve + spread_curve               # additive combination of the two curves

discount(yield,1.5) # 1 / (1 + 0.064 + 0.014) ^ 1.5
```

## Usage

### Rates

Rates are types that wrap scalar values to provide information about how to determine `discount` and `accumulation` factors.

There are two `CompoundingFrequency` types:

- `Yields.Periodic(m)` for rates that compound `m` times per period (e.g. `m` times per year if working with annual rates).
- `Yields.Continuous()` for continuously compounding rates.

#### Examples

```julia
Rate(0.05,Yields.Continuous())       # 5% continuously compounded
Yields.Continuous(0.05)              # alternate constructor

Rate(0.05, Yields.Periodic(2))       # 5% compounded twice per period
Yields.Periodic(0.05, 2)             # alternate constructor

# construct a vector of rates with the given compounding
Rate.(0.02,0.03,0.04,Yields.Periodic(2)) 
```

### Yields

There are a several ways to construct a yield curve object. `rates` can be a vector of `Rate`s described above, or will assume `Yields.Periodic(1)` if the functions are given `Real` number values

- `Yields.Zero(rates,maturities)`  using a vector of zero, or spot, rates
- `Yields.Forward(rates,maturities)` using a vector of one-period (or `periods`-long) forward rates
- `Yields.Constant(rate)` takes a single constant rate for all times
- `Yields.Step(rates,maturities)` doesn't interpolate - the rate is flat up to the corresponding time in `times`
- `Yields.Par(rates,maturities)` takes a series of yields for securities priced at par.Assumes that maturities <= 1 year do not pay coupons and that after one year, pays coupons with frequency equal to the CompoundingFrequency of the corresponding rate.
- `Yields.CMT(rates,maturities)` takes the most commonly presented rate data (e.g. [Treasury.gov](https://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=yield)) and bootstraps the curve given the combination of bills and bonds.
- `Yields.CMT(rates,maturities)` takes the most commonly presented rate data (e.g. [Treasury.gov](https://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=yield)) and bootstraps the curve given the combination of bills and bonds.

### Functions

Most of the above yields have the following defined (goal is to have them all):

- `discount(curve,from,to)` or `discount(curve,to)` gives the discount factor
- `accumulation(curve,from,to)` or `accumulation(curve,to)` gives the accumulation factor
- `forward(curve,from,to)` gives the average rate between the two given times
- `zero(curve,time)` or `zero(curve,time,CompoundingFrequency)` gives the zero-coupon spot rate for the given time.

### Combinations

Different yield objects can be combined with addition or subtraction. See the [Quickstart](#quickstart) for an example. 

When adding a `Yields.AbstractYield` with a scalar or vector, that scalar or vector will be promoted to a yield type via [`Yield()`](#yield). For example:

```julia
y1 = Yields.Constant(0.05)
y2 = y1 + 0.01                # y2 is a yield of 0.06
```

### Conversion

Convert rates between different types with `convert`. E.g.:

```julia-repl
r = Rate(Yields.Periodic(12),0.01)             # rate that compounds 12 times per rate period (ie monthly)

convert(Yields.Periodic(1),r)                  # convert monthly rate to annual effective
convert(Yields.Continuous(),r)          # convert monthly rate to continuous
```

## Internals

For time-variant yields (ie yield *curves*), the inputs are converted to spot rates and linearly interpolated (using [`Interpolations.jl`](https://github.com/JuliaMath/Interpolations.jl)). 

If you want more precise curvature (e.g. cubic spline interpolation) you can pre-process your rates into a greater number of input points before creating the `Yields` representation. `Yields.jl` uses `Interpolations.jl` as it is a pure-Julia interpolations package and enables auto-differentiation (AD) in `Yields.jl` usage. For example, [`ActuaryUtilities.jl`](https://github.com/JuliaActuary/ActuaryUtilities.jl) uses AD for `duration` and `convexity`.

### Combination Implementation

[Combinations](#combinations) track two different curve objects and are not combined into a single underlying data structure. This means that you may achieve better performance if you combine the rates before constructing a `Yields` representation. The exception to this is `Constant` curves, which *do* get combined into a single structure that is as performant as pre-combined rate structure.

## Related Packages 

- [**`InterestRates.jl`**](https://github.com/felipenoris/InterestRates.jl) specializes in fast rate calculations aimed at valuing fixed income contracts, with business-day-level accuracy. 
  - Comparative comments: **`Yields.jl`** does not try to provide as precise controls over the timing, structure, and interpolation of the curve. Instead, **`Yields.jl`** provides a minimal interface for common modeling needs.
