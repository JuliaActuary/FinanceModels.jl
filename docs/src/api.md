# Yields API Reference

Please [open an issue](https://github.com/JuliaActuary/Yields.jl/issues) if you encounter any issues or confusion with the package.

## Rate Types

When JuliaActuary packages return a rate, they will be of a `Rate` type, such as `Rate(0.05,Periodic(2))` for a 5% rate compounded twice per period. It is recommended to keep rates typed and use them throughout the ecosystem without modifying it. 

For example, if we construct a curve like this:

```julia
# 2021-03-31 rates from Treasury.gov
rates =[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
mats = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]
  
curve = Yields.CMT(rates,mats)
```

Then rates from this curve will be typed. For example:

```julia
z = zero(c,10)
```

Now, `z` will be: `Yields.Rate{Float64, Continuous}(0.01779624378877313, Continuous())`

This `Rate` has both the rate an the compounding convention embedded in the datatype.

You can now use that rate throughout the JuliaActuary ecosystem, such as with ActuaryUtilities.jl:

```julia
using ActuaryUtilities
present_values(z,cashflows)
```

If you need to extract the rate for some reason, you can get the rate by calling `Yields.rate(...)`. Using the above example, `Yields.rate(z)` will return `0.01779624378877313`. 

```@index
```

```@autodocs
Modules = [Yields]
```
