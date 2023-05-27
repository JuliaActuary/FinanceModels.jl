```@meta
CurrentModule = FinanceModels
```

# Developer Notes

## Custom Curve Types

Types that subtype `FinanceModels.AbstractYieldCurve` should implement a few key methods:

- `discount(curve,to)` should return the discount factor for the given curve through time `to`

For example:

```julia
struct MyYield <: FinanceModels.AbstractYieldCurve
    rate
end

FinanceModels.discount(c::MyYield,to) = exp(-c.rate * to)
```


By defining the `discount` method as above and subtyping `FinanceModels.AbstractYieldCurve`, FinanceModels.jl has generic functions that will work:

- `zero(curve,to)` returns the zero rate at time `to`
- `discount(curve,from,to)` is the discount factor between the two timepoints
- `forward(curve,to)` and `forward(curve,from,to)` is the forward zero rate
- `accumulation(curve,to)` and `accumulation(curve,from,to)` is the inverse of the discount factor.

If creating a new type of curve, you may find that it's most natural to define one of the functions versus the other, or that you may define specialized functions which are more performant than the generic implementations.

### `__ratetype`

In some contexts, such as creating performant iteration of curves in [EconomicScenarioGenerators.jl](https://github.com/JuliaActuary/EconomicScenarioGenerators.jl), Julia wants to know what type should be expected given an object type. For this reason, we define an internal, un-exported function which returns the `Rate` type expected given a Yield curve.

Sometimes it is most natural or convenient to expect a certain kind of `Rate` from a given curve. In many advanced use-cases (differentiation, stochastic rates), `Continuous` rates are most natural. For this reason, the `DEFAULT_COMPOUNDING` constant within FinanceModels.jl is $(FinanceModels.DEFAULT_COMPOUNDING). Two comments on this:

1. Becuase FinanceModels.jl returns `Rate` types (e.g. `Rate(0.05,Continuous()`) instead of single scalars (e.g. `0.05`) functions within the `JuliaActuary` universe (e.g. `ActuaryUtilities.present_value) know how to treat rates differently and in general users should not ever need to worry about converting between different compounding conventions.
2. Developers implementing new `AbstractYieldCurve` types can define their own default. For example, using the `MyYield` example above:

  - `__ratetype(::Type{MyYield}) = FinanceModels.Rate{Float64, Continuous}`

If the `CompoundingFrequency` is `Continuous`, then it's currently not necessary to define `__ratetype`, as it will fall back onto the generic method defined for `AbstractYieldCurve`s.

If the preferred compounding frequency is `Periodic`, then you must either define the methods (`zero`, `forward`,...) for your type or to use the generic methods then you must define `FinanceModels.CompoundingFrequency(curve::MyCurve)` to return the `Periodic` compounding datatype of the rates to return. 

For example, if we wanted `MyCurve` to return `Periodic(1)` rates, then we would define:

`FinanceModels.CompoundingFrequency(curve::MyCurve) = Periodic(1)`

This latter step is necessary and distinct from `__ratetype`. This is due to `__ratetype` relying on type-only information. The `Periodic` type contains as a datafield the compounding frequency. Therefore, the frequency is not known to the type system and is available only at runtime. 
