```@meta
CurrentModule = Yields
```

# Developer Notes

## Custom Curve Types

Types that subtype `Yields.AbstractYield` should implement a few key methods:

- `discount(curve,to)` should return the discount factor for the given curve through time `to`

For example:

```julia
struct MyYield <: Yields.AbstractYield
    rate
end

Yields.discount(c::MyYield,to) = exp(-c.rate * to)
```


By defining the `discount` method as above and subtyping `Yields.AbstractYield`, Yields.jl has generic functions that will work:

- `zero(curve,to)` returns the zero rate at time `to`
- `discount(curve,from,to)` is the discount factor between the two timepoints
- `forward(curve,to)` and `forward(curve,from,to)` is the forward zero rate
- `accumulation(curve,to)` and `accumulation(curve,from,to)` is the inverse of the discount factor.

If creating a new type of curve, you may find that it's most natural to define one of the functions versus the other, or that you may define specialized functions which are more performant than the generic implementations.

### `__ratetype`

In some contexts, such as creating performant iteration of curves in [EconomicScenarioGenerators.jl](https://github.com/JuliaActuary/EconomicScenarioGenerators.jl), Julia wants to know what type should be expected given an object type. For this reason, we define an internal, un-exported function which returns the `Rate` type expected given a Yield curve.

Sometimes it is most natural or convenient to expect a certain kind of `Rate` from a given curve. In many advanced use-cases (differentiation, stochastic rates), `Continuous` rates are most natural. For this reason, the `DEFAULT_COMPOUNDING` constant within Yields.jl is $(Yields.DEFAULT_COMPOUNDING). Two comments on this:

1. Becuase Yields.jl returns `Rate` types (e.g. `Rate(0.05,Continuous()`) instead of single scalars (e.g. `0.05`) functions within the `JuliaActuary` universe (e.g. `ActuaryUtilities.present_value) know how to treat rates differently and in general users should not ever need to worry about converting between different compounding conventions.
2. Developers implementing new `AbstractYield` types can define their own default. For example, using the `MyYield` example above:

  - `__ratetype(::Type{MyYield}) = Yields.Rate{Float64, Continuous}`
