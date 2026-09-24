# A curve view that evaluates without dual numbers, so internal derivatives never
# mix a caller's ForwardDiff tag with an internal one.
struct __PrimalCurve{C} <: AbstractYieldModel
    curve::C
end
FinanceCore.discount(c::__PrimalCurve, t) = __primal(FinanceCore.discount(c.curve, t))
function Base.zero(c::__PrimalCurve, t)
    z = convert(Continuous(), Base.zero(c.curve, t))
    return Continuous(__primal(FinanceCore.rate(z)))
end

"""
    implied_quote(curve, family, maturity; guess = 0.0, bracket = (-0.5, 1.0))

Return the quote `x` for which `family(x, maturity)` reprices on `curve`, so that
`present_value(curve, q.instrument) == q.price` for `q = family(x, maturity)`.

`family` is a quote constructor taking `(quote, maturity)`, such as
[`CMTYield`](@ref FinanceModels.Bond.CMTYield), [`OISYield`](@ref FinanceModels.Bond.OISYield),
[`ZCBYield`](@ref FinanceModels.Bond.ZCBYield), [`ZCBPrice`](@ref FinanceModels.Bond.ZCBPrice),
or a closure like `(r, t) -> ParYield(r, t; frequency = 1)`. The result is expressed in
the family's own convention: for example an annual-effective rate for `ZCBYield`,
a semiannual par yield for `CMTYield` beyond one year, or a price for `ZCBPrice`.

The solve starts from `guess` and falls back to a bracketed search on `bracket`.
First-order ForwardDiff derivatives with respect to curve parameters are exact:
they come from the implicit function theorem at the solution, not from solver
iterations. Nested dual numbers, or a `family` that closes over dual numbers, throw
an `ArgumentError`.

# Examples

```julia-repl
julia> curve = Yield.Constant(0.04);

julia> implied_quote(curve, ZCBYield, 5.0) ≈ 0.04
true

julia> implied_quote(curve, (r, t) -> ParYield(r, t; frequency = 1), 5.0) ≈ 0.04
true
```

See also [`par`](@ref).
"""
function implied_quote(curve, family::F, maturity; guess = 0.0, bracket = (-0.5, 1.0)) where {F}
    primal = __PrimalCurve(curve)
    residual(c, x) = (q = family(x, maturity); FinanceCore.present_value(c, q.instrument) - q.price)
    g(x) = residual(curve, x)
    g_primal(x) = residual(primal, x)
    __ad_depth(g_primal(float(guess))) == 0 || throw(
        ArgumentError(
            "implied_quote differentiates through the curve only; the quote family must not close over dual numbers"
        )
    )
    # the size of the quote's price and value at the solution, whatever its notional
    scale(x) = (q = family(x, maturity); max(abs(FinanceCore.present_value(primal, q.instrument)), abs(q.price)))
    return __implicit_root(g, g_primal, guess; bracket, who = "implied_quote", scale)
end
