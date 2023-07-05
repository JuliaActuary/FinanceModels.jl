# FinanceModels.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaActuary.github.io/FinanceModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaActuary.github.io/FinanceModels.jl/dev)
[![Build Status](https://github.com/JuliaActuary/FinanceModels.jl/workflows/CI/badge.svg)](https://github.com/JuliaActuary/FinanceModels.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaActuary/FinanceModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaActuary/FinanceModels.jl)

**FinanceModels.jl** provides a simple interface for constructing, manipulating, and using yield curves for modeling purposes.

It's intended to provide common functionality around modeling interest rates, spreads, and miscellaneous FinanceModels across the JuliaActuary ecosystem (though not limited to use in JuliaActuary packages).

![anim_fps2](https://user-images.githubusercontent.com/711879/174458687-860c5d7f-e125-46a9-a706-7d113f1e243b.gif)


## QuickStart

```julia
using FinanceModels

riskfree_maturities = [0.5, 1.0, 1.5, 2.0]
riskfree    = [5.0, 5.8, 6.4, 6.8] ./ 100     #spot rates, annual effective if unspecified

spread_maturities = [0.5, 1.0, 1.5, 3.0]      # different maturities
spread    = [1.0, 1.8, 1.4, 1.8] ./ 100       # spot spreads

rf_curve = FinanceModels.Zero(riskfree,riskfree_maturities)
spread_curve = FinanceModels.Zero(spread,spread_maturities)


yield = rf_curve + spread_curve               # additive combination of the two curves

discount(yield,1.5)                           # 1 / (1 + 0.064 + 0.014) ^ 1.5
```

## Usage

### Rates

Rates are types that wrap scalar values to provide information about how to determine `discount` and `accumulation` factors.

There are two `Frequency` types:

- `FinanceModels.Periodic(m)` for rates that compound `m` times per period (e.g. `m` times per year if working with annual rates).
- `FinanceModels.Continuous()` for continuously compounding rates.

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

Rates can also be constructed by specifying the `Frequency` and then passing a scalar rate:

```julia
Periodic(1)(0.05)
Continuous()(0.05)
```

#### Conversion

Convert rates between different types with `convert`. E.g.:

```julia-repl
r = Rate(FinanceModels.Periodic(12),0.01)             # rate that compounds 12 times per rate period (ie monthly)

convert(FinanceModels.Periodic(1),r)                  # convert monthly rate to annual effective
convert(FinanceModels.Continuous(),r)          # convert monthly rate to continuous
```

#### Arithmetic

Adding, substracting, multiplying, dividing, and comparing rates is supported.

### Curves

There are a several ways to construct a yield curve object. If `maturities` is omitted, the method will assume that the timepoints corresponding to each rate are the indices of the `rates` (e.g. generally one to the length of the array for standard, non-offset arrays). 

#### Fitting Curves to Rates

There is a set of constructor methods which will return a yield curve calibrated to the given inputs. 

- `FinanceModels.Zero(rates,maturities)`  using a vector of zero rates (sometimes referred to as "spot" rates)
- `FinanceModels.Forward(rates,maturities)` using a vector of forward rates
- `FinanceModels.Par(rates,maturities)` takes a series of FinanceModels for securities priced at par. Assumes that maturities <= 1 year do not pay coupons and that after one year, pays coupons with frequency equal to the Frequency of the corresponding rate (2 by default).
- `FinanceModels.CMT(rates,maturities)` takes the most commonly presented rate data (e.g. [Treasury.gov](https://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=yield)) and bootstraps the curve given the combination of bills and bonds.
- `FinanceModels.OIS(rates,maturities)` takes the most commonly presented rate data for overnight swaps and bootstraps the curve. Rates assume a single settlement for <1 year and quarterly settlements for 1 year and above.

##### Fitting techniques

There are multiple curve fitting methods available:

- `Boostrap(interpolation_method)` (the default method)
  - where `interpolation` can be one of the built-in `QuadraticSpline()` (the default) or `LinearSpline()`, or a user-supplied function.
- Two methods from the Nelson-Siegel-Svensson family, where τ_initial is the starting τ point for the fitting optimization routine: 
  - `NelsonSiegel(τ_initial=1.0)`
  - `NelsonSiegelSvensson(τ_initial=[1.0,1.0])`

To specify which fitting method to use, pass the object to as the first parameter to the above set of constructors, for example: `FinanceModels.Par(NelsonSiegel(),rates,maturities)`.

#### Kernel Methods

- `FinanceModels.SmithWilson` curve (used for [discounting in the EU Solvency II framework](https://www.eiopa.europa.eu/sites/default/files/risk_free_interest_rate/12092019-technical_documentation.pdf)) can be constructed either directly by specifying its inner representation or by calibrating to a set of cashflows with known prices.
  - These cashflows can conveniently be constructed with a Vector of `FinanceModels.ZeroCouponQuote`s, `FinanceModels.SwapQuote`s, or `FinanceModels.BulletBondQuote`s.

#### Other Curves

- `FinanceModels.Constant(rate)` takes a single constant rate for all times
- `FinanceModels.Step(rates,maturities)` doesn't interpolate - the rate is flat up to the corresponding time in `times`

### Functions

Most of the above FinanceModels have the following defined (goal is to have them all):

- `discount(curve,from,to)` or `discount(curve,to)` gives the discount factor
- `accumulation(curve,from,to)` or `accumulation(curve,to)` gives the accumulation factor
- `zero(curve,time)` or `zero(curve,time,Frequency)` gives the zero-coupon spot rate for the given time.
- `forward(curve,from,to)` gives the zero rate between the two given times
- `par(curve,time)` gives the coupon-paying par equivalent rate for the given time.

### Combinations

Different yield objects can be combined with addition or subtraction. See the [Quickstart](#quickstart) for an example.

When adding a `FinanceModels.AbstractYield` with a scalar or vector, that scalar or vector will be promoted to a yield type via [`Yield()`](#yield). For example:

```julia
y1 = FinanceModels.Constant(0.05)
y2 = y1 + 0.01                # y2 is a yield of 0.06
```

### Forward Starting Curves

Constructed curves can be shifted so that a future timepoint becomes the effective time-zero for a said curve.

```julia-repl
julia> zero = [5.0, 5.8, 6.4, 6.8] ./ 100
julia> maturity = [0.5, 1.0, 1.5, 2.0]
julia> curve = FinanceModels.Zero(zero, maturity)
julia> fwd = FinanceModels.ForwardStarting(curve, 1.0)

julia> discount(curve,1,2)
0.9275624570410582

julia> discount(fwd,1) # `curve` has effectively been reindexed to `1.0`
0.9275624570410582
```

## Exported vs Un-exported Functions

Generally, CamelCase methods which construct a datatype are exported as they are unlikely to conflict with other parts of code that may be written. For example, `rate` is un-exported (it must be called with `FinanceModels.rate(...)`) because `rate` is likely a very commonly defined variable within actuarial and financial contexts and there is a high risk of conflicting with defined variables.

Consider using `import FinanceModels` which would require qualifying all methods, but alleviates any namespace conflicts and has the benefit of being explicit about the calls (internally we prefer this in the package design to keep dependencies and their usage clear). 

## Internals

For time-variant FinanceModels (ie yield *curves*), the inputs are converted to spot rates and interpolated using quadratic B-splines by default (see documentation for alternatives, such as linear interpolations).

### Combination Implementation

[Combinations](#combinations) track two different curve objects and are not combined into a single underlying data structure. This means that you may achieve better performance if you combine the rates before constructing a `FinanceModels` representation. The exception to this is `Constant` curves, which *do* get combined into a single structure that is as performant as pre-combined rate structure.

## Related Packages

- [**`InterestRates.jl`**](https://github.com/felipenoris/InterestRates.jl) specializes in fast rate calculations aimed at valuing fixed income contracts, with business-day-level accuracy.
  - Comparative comments: **`FinanceModels.jl`** does not try to provide as precise controls over the timing, structure, and interpolation of the curve. Instead, **`FinanceModels.jl`** provides a minimal, but flexible and intuitive interface for common modeling needs.
