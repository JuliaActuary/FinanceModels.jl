# Yields

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaActuary.github.io/Yields.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaActuary.github.io/Yields.jl/dev)
[![Build Status](https://github.com/JuliaActuary/Yields.jl/workflows/CI/badge.svg)](https://github.com/JuliaActuary/Yields.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaActuary/Yields.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaActuary/Yields.jl)
[![lifecycle](https://img.shields.io/badge/LifeCycle-Experimental-orange)](https://www.tidyverse.org/lifecycle/)


**Yields** provides a simple interface for constructing, manipulating, and using yield curves for modeling purposes.

It's intended to provide common functionality around modeling interest rates, spreads, and miscellaneous yields across the JuliaActuary ecosystem (though not limited to use in JuliaActuary packages).

## QuickStart

```julia
using Yields

riskfree_maturities = [0.5, 1.0, 1.5, 2.0]
riskfree    = [5.0, 5.8, 6.4, 6.8] ./ 100     #spot rates

spread_maturities = [0.5, 1.0, 1.5, 3.0]      # different maturities
spread    = [1.0, 1.8, 1.4, 1.8] ./ 100       # spot spreads

rf_curve = Yields.Zero(riskfree,riskfree_maturities)
spread_curve = Yields.Zero(spread,spread_maturities)


yield = rf_curve + spread_curve               # additive combination of the two curves

discount(yield,1.0) # 1 / (1 + 0.058 + 0.018)
```

## Usage

### Yields

There are a few ways to construct a yield curve object:

- `Zero(rates)` or `Zero(rates,maturities)` using a vector of zero, or spot, rates
- `Forward(rates)` or `Forward(rates,periods) using a vector of one-period (or `periods`-long) forward rates
- `Constant(rate)` takes a single constant rate for all times
- `USTreasury(rates)` takes the most commonly presented rate data (e.g. [Treasury.gov](https://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=yield)) and bootstraps the curve given the combination of bills and bonds.

### Functions

Most of the above yields have the following defined (goal is to have them all):
- `rate(curve,time)` gives the rate at `time` 
- `discount(curve,time)` gives the discount factor through `time`
- `accumulate(curve,time)` gives the accumulation factor through `time`
- `forward(curve,time_from,time_to)` gives the average rate between the two given times

### Combinations

Different yield objects can be combined with addition or subtraction. See the [Quickstart](#quickstart) for an example.


## Related Packages 

- [**`InterestRates.jl`**](https://github.com/felipenoris/InterestRates.jl) specializes in fast rate calculations aimed at valuing fixed income contracts, with business-day-level accuracy. 
  - Comparative comments: **`Yields.jl`** does not try to provide as precise controls over the timing, structure, and interpolation of the curve. Instead, **`Yields.jl`** provides a minimal interface for common modeling needs.
