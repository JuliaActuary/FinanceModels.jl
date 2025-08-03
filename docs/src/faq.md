# Frequently Asked Questions

## How can I handle `Date`s instead of real timepoints?

Currently, you must convert the `Date` into a real-valued timepoint for use within the models and contracts. Future releases may contemplate more explicit, built-in handling of dates. In the meantime, you may find these packages helpful if you need precise date-level accuracy:

- [Miletus.jl](https://github.com/JuliaComputing/Miletus.jl)
- [InterestRates.jl](https://github.com/felipenoris/InterestRates.jl)
- [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl)
- [DayCounts.jl](https://github.com/JuliaFinance/DayCounts.jl)
- [QuantLib.jl](https://github.com/pazzo83/QuantLib.jl)

## Why does the package rely on using Transducers?

Transducers are a way of defining logic to be applied to a reducible collection. They can compose together efficiently and the compiler can optimize them well. In rewriting the package from v3 to v4, Transducers vastly simplified the iteration and state handling needed when projecting the contracts. The performance remains excellent and made a lot of the internals much simpler.

Transducers are a [rich and powerful way to express programs](https://www.youtube.com/watch?v=6mTbuzafcII) and can seem somewhat unfamiliar at first encounter. For users of FinanceModels, very of transducers are needed/exposed:

- To regular end-users who just use what is given to them here, the transducers internals are effectively completely hidden
- To moderately advanced users who want to extend the functionality, as the examples show the only real exposure here is a weird function name ( `__foldl__`) with for loop with a `return` signature that has some extra information.

A number of examples of extending the package are given on the [FinanceModels.jl Guide](@ref) page and the of course the source code itself offers examples of existing `Projection`s and `Contract`s.

## Composite Yield/Discount Curves

Curves can be added or subtracted together, but note that this is not always the same thing as adding or subtracting spreads with rates. If spreads and base rates are expressed as zero rates, then the curve addition/subtraction has the same effect as re-fitting the yield model with the rate+spread inputs added together first. Non-zero rates (e.g. par rates) do not have this same property. Zero-coupon rates have a direct, linear relationship with the underlying discount factors. Par-coupon rates have a complex, non-linear relationship with the underlying discount factors and so the curve addition/subtraction does not work the same way.

Example:

```julia
using FinanceModels
using Test


rates = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
spreads = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
mats = [1 / 12, 2 / 12, 3 / 12, 6 / 12, 1, 2, 3, 5, 7, 10, 20, 30]


### Zero coupon rates/spreads

q_rf_z = ZCBYield.(rates,mats)
q_s_z = ZCBYield.(spreads,mats)
q_y_z = ZCBYield.(rates + spreads,mats)

c_rf_z = fit(Spline.Linear(),q_rf_z,Fit.Bootstrap())
c_s_z = fit(Spline.Linear(),q_s_z,Fit.Bootstrap())
c_y_z = fit(Spline.Linear(),q_y_z,Fit.Bootstrap())

# adding curves when the spreads were zero spreads DOES works
discount(c_rf_z+c_s_z,20) ≈ discount(c_y_z,20) #true


### Par coupon rates/spreads

q_rf = CMTYield.(rates,mats)
q_s = CMTYield.(spreads,mats)
q_y = CMTYield.(rates + spreads,mats)

c_rf = fit(Spline.Linear(),q_rf,Fit.Bootstrap())
c_s = fit(Spline.Linear(),q_s,Fit.Bootstrap())
c_y = fit(Spline.Linear(),q_y,Fit.Bootstrap())

# adding curves when the spreads were par spreads does NOT work
discount(c_rf+c_s,20) ≈ discount(c_y,20) # false



```



## I have another question

Ask on the discussion forum here: [https://github.com/JuliaActuary/FinanceModels.jl/discussions](https://github.com/JuliaActuary/FinanceModels.jl/discussions)

