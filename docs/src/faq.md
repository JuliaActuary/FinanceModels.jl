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

Transducers are a [rich and powerful way to express programs](https://www.youtube.com/watch?v=6mTbuzafcII) and can seem somewhat unfamiliar at first encounter. For users of FinanceModels, few details of Transducers are exposed:

- To regular end-users who just use what is given to them here, the transducers internals are effectively completely hidden
- To moderately advanced users who want to extend the functionality, as the examples show the only real exposure here is a weird function name ( `__foldl__`) with for loop with a `return` signature that has some extra information.

A number of examples of extending the package are given on the [FinanceModels.jl Guide](@ref) page and, of course, the source code itself offers examples of existing `Projection`s and `Contract`s.

## Composite Yield/Discount Curves

Curve addition multiplies discount factors, which adds **continuously compounded zero rates**. Discount factors depend exponentially on zero rates, not linearly. Adding annual effective rates or par yields before fitting generally gives a different result.

For linear interpolation on a common knot grid, interpolating the sum of continuous zero rates is equivalent to adding the separately interpolated curves. This identity does not generally extend to nonlinear interpolation methods such as PCHIP or MonotoneConvex.

```@example composition
using FinanceModels

rates = [0.01, 0.02, 0.03]
spreads = [0.005, 0.01, 0.015]
tenors = [1.0, 5.0, 10.0]

base = ZeroRateCurve(Continuous.(rates), tenors, Spline.Linear())
spread = ZeroRateCurve(Continuous.(spreads), tenors, Spline.Linear())
combined = ZeroRateCurve(Continuous.(rates + spreads), tenors, Spline.Linear())
@assert discount(base + spread, 7.0) ≈ discount(combined, 7.0) # hide
discount(base + spread, 7.0) ≈ discount(combined, 7.0)
```

Par yields depend on all coupon-date discount factors. If the desired quotes are base par yields plus par spreads, combine those quotes first and fit a single curve:

```@example composition
base_par = fit(Spline.Linear(), ParYield.(rates, tenors))
spread_par = fit(Spline.Linear(), ParYield.(spreads, tenors))
combined_par = fit(Spline.Linear(), ParYield.(rates + spreads, tenors))
@assert !isapprox(discount(base_par + spread_par, 7.0), discount(combined_par, 7.0)) # hide
discount(base_par + spread_par, 7.0) ≈ discount(combined_par, 7.0) # false
```

To apply par spreads over an existing curve without refitting it, use [`quote_shift`](@ref FinanceModels.Yield.quote_shift); it meets the shifted quotes at the quoted tenors only.

## I have another question

Ask on the discussion forum here: [https://github.com/JuliaActuary/FinanceModels.jl/discussions](https://github.com/JuliaActuary/FinanceModels.jl/discussions)
