```@meta
CurrentModule = Yields
```


# Upgrading from Prior Versions

## v3 to v4

### tl;dr

#### Use the `curve` method

There's now just a single primary constructor method named [`curve`](@ref)

`Yields.Par(rates,maturities)` -> `curve(ParYield.(rates,maturities))`
`Yields.CMT(rates,maturities)` -> `curve(CMTYield.(rates,maturities))`
`Yields.Zero(rates,maturities)` -> `curve(ZCBYield.(rates,maturities))`
`Yields.Forward(rates,maturities)` -> `curve(ForwardYield.(rates,maturities))`

The resulting curves can be used just as before.

#### Picking different `CurveMethod`]s

Pass [`CurveMethod`](@ref) to the first argument in [`curve`](@ref). The default is still `Bootstrap(QuadraticSpline())`. Examples:

```julia
curve(ZCBYield.(rates,maturities))
curve(Bootstrap(),ZCBYield.(rates,maturities))
curve(Bootstrap(LinearSpline()),ZCBYield.(rates,maturities))
curve(SmithWilson(ufr, α),ZCBYield.(rates,maturities))
curve(NelsonSiegel(τ_initial),ZCBYield.(rates,maturities))
curve(NelsonSiegel(τ_initial),ZCBYield.(rates,maturities))
```

### More details on the changes

The package has been enhanced in a couple of ways:

1. A more explicit and generic way to represent financial instruments [`Instruments`](@ref) that are used to generate curves of different kinds.
2. A simplified interface utilizing multiple dispatch: the [`curve`](@ref) method.

#### Instruments

These provide an explicit way to specify what a quoted price or yield represents. This allows simplified curve construction logic because it no longer needs to explicitly handle all types of rates (i.e. the old `Par`, `CMT`, `OIS`, `Forward`, and `Zero`).

For more details, see the [Instruments and Quotes](@ref) page.

#### `curve` method

A single primary method allows:

1. a simpler, unified interface
2. simplified logic by utilizing multiple dispatch on the new [`Instruments`](@ref)
3. Easier to extend package for a user-defined `CurveMethod`s