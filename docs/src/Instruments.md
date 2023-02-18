```@meta
CurrentModule = Yields
```

# Instruments and Quotes

## Instruments

Instruments are a way to represent basic units of traded instruments and cashflows.

Instrument Types:

```plaintext
Instrument  #an  abstract type
 - Cashflow(amount,time)
 - Bond(coupon_rate,frequency,maturity)
 - Forward(time,Instrument)
 - Composite(Instrument1,Instrument2)
```

These immutable `struct`s provide a generic and composable way to represent a wide range of financial instruments. By design, they are not intended to fully represent complete attributes of traded securities. If you have modeled more complex securities, it should be simple to translate them into the data types used in JuliaActuary packages.


### `Cashflow`

A data type indicating an amount and a time for a cashflow.

### `Bond`

A data type represent a bullet or zero coupon bond. A `Bond` can be decomposed into its constituent `Cashflow`s by calling `collect(Bond(coupon_rate,frequency,maturity))`

## Quotes

The above instruments are combined with `Quote`s to represent an observed price or yield.



## Convenience Functions

Utility functions are provided to make it easy to define representative instruments without any boilerplate:

```@docs
ZCBPrice
ZCBYield
ParYield
ParSwapYield
CMTYield
OISYield
ForwardYield
```

## Time

All of the times represented in the datatype represent a time relative to time zero. E.g. if `time=1.0` in `Cashflow`, that means 1 unit period from time zero. If you need to manage forward instruments relative to time zero, that's what `Forward(time,Instrument)` does. If you need to represent instruments held at future points in time (i.e. "future time zeros") then you can use a variety of Julia datastructures to handle them. For example, create an array or dictionary containing the sets of instruments you need. For example:


```julia
[
    # outer loop time zero
    ZCBYield.([5.0, 5.8, 6.4, 6.8] ./ 100, [0.5, 1.0, 1.5, 2.0]),

    # outer loop time one
    ZCBYield.([5.2, 5.9, 6.5, 6.9] ./ 100, [0.5, 1.0, 1.5, 2.0]),
]
```

Also note that [`ForwardStarting`](@ref) is a way to create curves implied by a given curve at a forward point in time.