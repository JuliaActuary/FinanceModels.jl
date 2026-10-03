# A curve view that evaluates without dual numbers, so internal derivatives never mix a caller's
# ForwardDiff tag with an internal one. Like `HullWhite`, it forwards every capability in
# `__FORWARDED_CAPABILITIES`.
struct __PrimalCurve{C} <: AbstractYieldModel
    curve::C
end
FinanceCore.discount(c::__PrimalCurve, t) = __primal(FinanceCore.discount(c.curve, t))
FinanceCore.discount(c::__PrimalCurve, from, to) = __primal(FinanceCore.discount(c.curve, from, to))
__log_discount(c::__PrimalCurve, t) = __primal(__log_discount(c.curve, t))
__log_interval(c::__PrimalCurve, from, to) = __primal(__log_interval(c.curve, from, to))
function __log_tail(c::__PrimalCurve)
    t = __log_tail(c.curve)
    return __LogTail(__primal(t.a2), __primal(t.a1), __primal(t.a0))
end
Base.zero(c::__PrimalCurve, t) = Continuous(__primal(__continuous(Base.zero(c.curve, t))))
instantaneous_forward(c::__PrimalCurve, t) = __primal(instantaneous_forward(c.curve, t))

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
    # the size of the quote's price and value, whatever its notional
    magnitude(x) = (q = family(x, maturity); max(abs(FinanceCore.present_value(primal, q.instrument)), abs(q.price)))
    # Solve relative to the quote's size at the guess, so the solvers' absolute tolerances do
    # not depend on units such as a notional. A zero size means the price and the value both
    # vanish at the guess, which is then an exact root.
    s0 = magnitude(float(guess))
    σ = iszero(s0) ? one(s0) : s0
    g(x) = residual(curve, x) / σ
    g_primal(x) = residual(primal, x) / σ
    scale(x) = magnitude(x) / σ
    return __implicit_root(
        g, g_primal, guess; bracket, scale, who = "implied_quote",
        hint = "implied_quote differentiates through the curve only, not through `family` or `guess`."
    )
end

"""
    par(curve,time;frequency=2)

Calculate the par yield for maturity `time` for the given `curve` and `frequency`. Returns a `Rate` object with periodicity corresponding to the `frequency`.

If `time` is shorter than one regular coupon period (e.g. `time=0.5` with `frequency=1`), the single stub payment implies a compounding frequency of `1/time`: the result is quoted as `Periodic(1/time)` when `1/time` is a (near-)integer, and otherwise an `ArgumentError` is thrown because the implied frequency cannot be represented as a `Periodic` rate.

If `time` is longer than one coupon period but not a whole number of periods, the schedule has a short first stub which accrues its actual length (see `Bond.coupon_times`); the result is the internal rate of return of the par-priced true-accrual schedule, quoted as `Periodic(frequency)` — so `par` of a flat curve recovers the curve's rate at any maturity. On such stub schedules this yield quote differs from the annualized par *coupon* `c` solving `c·Σᵢ δᵢ·DF(tᵢ) + DF(T) = 1`, which is what `InterestRateSwap` uses for its fixed leg.

# Examples

```julia-repl
julia> c = Yield.Constant(0.04);

julia> par(c,4)
Periodic(0.03960780543711406, 2)

julia> par(c,4;frequency=1)
Periodic(0.040000000000000036, 1)

julia> par(c,0.6;frequency=4)
Periodic(0.039413626195875295, 4)

julia> par(c,0.2;frequency=4)
Periodic(0.039374942589460726, 5)

julia> par(c,2.5)
Periodic(0.03960780543711406, 2)
```
"""
function par(curve, time; frequency = 2)
    coup_times = coupon_times(time, frequency)
    mat_disc = discount(curve, time)
    coupon_pv = sum(discount(curve, t) for t in coup_times)
    Δt = step(coup_times)
    r = (1 - mat_disc) / coupon_pv

    # A maturity that is not a whole number of periods leaves a short first stub
    # which accrues its actual length: such schedules pay `c·δᵢ` with the
    # annualized coupon `c` solving the true-accrual par condition, so the IRR
    # below is that of the correctly-accrued par-priced bond. Whole-period
    # schedules pay `r` every period (the historical construction, unchanged).
    stub = length(coup_times) > 1 && !__regular_schedule(time, frequency)
    c = stub ? __par_coupon(curve, time, frequency) : r

    # Build cash flows: initial outflow of -1, then coupons, final coupon+principal 1+coupon
    n = length(coup_times)
    cfs = Vector{typeof(r)}(undef, n + 1)
    times = Vector{typeof(Δt)}(undef, n + 1)

    cfs[1] = -one(r)
    times[1] = zero(Δt)

    @inbounds for i in 1:n
        coup = stub ? c * (i == 1 ? first(coup_times) : 1 / frequency) : r
        cfs[i + 1] = i == n ? 1 + coup : coup
        times[i + 1] = coup_times[i]
    end

    r = FinanceCore.internal_rate_of_return(cfs, times)
    frequency_inner = 1 / Δt
    if !isinteger(round(frequency_inner, digits = 8))
        throw(
            ArgumentError(
                "par(curve, $time; frequency=$frequency) implies a coupon period of $Δt and a compounding frequency of 1/Δt = $frequency_inner, which is not an integer and cannot be represented as a `Periodic` rate. Choose a maturity commensurate with the coupon frequency."
            )
        )
    end
    r = convert(Periodic(frequency_inner), r)
    return r
end
