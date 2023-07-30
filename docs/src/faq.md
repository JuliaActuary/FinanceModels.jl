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

A number of examples of extending the package are given on the [Overview](/overview) page and the of course the source code itself offers examples of existing `Projection`s and `Contract`s.


## I have another question

Ask on the discussion forum here: https://github.com/JuliaActuary/Yields.jl/discussions

