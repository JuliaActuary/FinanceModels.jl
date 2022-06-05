```@meta
CurrentModule = Yields
```

# Developer Notes

Types that subtype `Yields.AbstractYield` should implement a few key methods:

- `discount(curve,to)` should return the discount factor for the given curve through time `to`

By defining the `discount` method as above and subtyping `Yields.AbstractYield`, Yields.jl has generic functions that will work:

- `zero(curve,to)` returns the zero rate at time `to`
- `discount(curve,from,to)` is the discount factor between the two timepoints
- `forward(curve,to)` and `forward(curve,from,to)` is the forward zero rate
- `accumulation(curve,to)` and `accumulation(curve,from,to)` is the inverse of the discount factor.

If creating a new type of curve, you may find that it's most natural to define one of the functions versus the other, or that you may define specialized functions which are more performant than the generic implementations.
