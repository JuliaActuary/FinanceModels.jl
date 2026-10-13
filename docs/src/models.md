# Models, Valuation, Projections, and Fitting

## Introduction

Conceptually, we have an iterative process:

1. We use models to value contracts
2. We use observed (or assumed) prices to calibrate models

Thus the discussion of model calibration and valuation of contracts is inextricably linked together.

## Yield (Interest Rate) models

### Rates

We should first discuss `Rate`s, which are reexported from [`FinanceCore.jl`](https://github.com/JuliaActuary/FinanceCore.jl)

Rates are types that wrap scalar values to provide information about how to determine `discount` and `accumulation` factors. These allow for explicit handling of rate compounding conventions which, if not explicit, is often a source of errors in practice.

There are two `Frequency` types:

- `Periodic(m)` for rates that compound `m` times per period (e.g. `m` times per year if working with annual rates).
- `Continuous()` for continuously compounding rates.

#### Examples

```julia
Continuous(0.05)       # 5% continuously compounded
Periodic(0.05,2)       # 5% compounded twice per period
```

These construct `Rate` values with the chosen compounding convention, equivalently:

```julia
Rate(0.05,Continuous())       # 5% continuously compounded
Rate(0.05,Periodic(2))        # 5% compounded twice per period
```

Broadcast over a vector to create `Rates` with the given compounding:

```julia
Periodic.([0.02,0.03,0.04],2) 
Continuous.([0.02,0.03,0.04]) 
```

Rates can also be constructed by specifying the `Frequency` and then passing a scalar rate:

```julia
Periodic(1)(0.05)
Continuous()(0.05)
```

#### Conversion

Convert rates between different types with `convert`. E.g.:

```julia
r = Rate(0.01,Periodic(12))             # rate that compounds 12 times per rate period (ie monthly)

convert(Periodic(1),r)                  # convert monthly rate to annual effective
convert(Continuous(),r)          # convert monthly rate to continuous
```

To get the scalar value out of the `Rate`, use `FinanceModels.rate(r)`:

```julia-repl
julia> r = Rate(0.01,Periodic(12));

julia> rate(r)
0.009999999999998899
```

(`Rate`s internally store the equivalent continuously compounded rate, so the nominal value returned by `rate` may show floating point artifacts of the round-trip conversion.)

#### [Rate conventions](@id rate-conventions)

A bare number has no compounding convention of its own, so each input reads it by its role:

| Quantity | Convention of a number | Examples |
|---|---|---|
| Curve coordinates: knot rates, zero rates, instantaneous forwards | Continuous | `ZeroRateCurve([0.03, 0.04], [1, 5])`; results `knot_rates(curve)`, `zero(curve, t)`, `Yield.instantaneous_forward(curve, t)` |
| Scalar rates for FinanceCore and flat curves | Annual effective, `Periodic(1)` | `Yield.Constant(0.05)`, `curve + 0.01`, `discount(0.05, t)`, `ZCBYield(0.05, t)` |
| Coupons, margins and interest-rate strikes | Nominal, paid at the contract's frequency | `Bond.Fixed(0.05, Periodic(2), 10)`, `Bond.Floating(0.002, Periodic(4), 5, :SOFR)`, `Option.Cap(0.03, 4, 5)` |

Rate-valued results are `Rate`s, which carry their convention; curve coordinates are `Continuous`.
Each of these inputs also takes a `Rate`:

- a curve coordinate or other continuous input (knot rates, Smith–Wilson's `ufr`, a conditional
  short rate) converts the `Rate` to its continuously compounded value, so
  `ZeroRateCurve(Periodic.([0.03, 0.04], 1), [1, 5])` builds the curve of those annual rates;
- a coupon, margin or strike must be a `Periodic` rate of the contract's frequency, and its
  nominal rate is used: `Bond.Fixed(Periodic(0.05, 2), Periodic(2), 10)` is
  `Bond.Fixed(0.05, Periodic(2), 10)`. Another frequency throws an `ArgumentError`, because its
  nominal rate and its converted rate give different coupons; convert it explicitly. A typed
  coupon can differ from the number in the last bits, since the `Rate` stores its continuously
  compounded equivalent.

A few quantities stay numbers: raw model fields (volatilities, mean-reversion speeds, curve
coefficients, the stored `ufr` and long-run means), quote coordinates, and sensitivities. Vasicek's
and CIR's `initial` short rate is a `Rate`: a number is stored as `Continuous`.
[`implied_quote`](@ref FinanceModels.Yield.implied_quote) returns a number in its quote family's
convention; wrap it when that convention is known, for example
`Periodic(implied_quote(curve, OISYield, 7.0), 1)`.

`rate` gives the number of a `Rate`, for example as a ForwardDiff input:
`rate.(knot_rates(curve))`. A `Rate` never equals a number, so `knot_rates(curve) == [0.03, 0.04]`
is `false`; compare with `Continuous.([0.03, 0.04])`.

### Available Models - Yields

- [`FinanceModels.Yield.Constant`](@ref)
- Fitted and bootstrapped [`Spline`](@ref FinanceModels.Spline)s
- [`FinanceModels.Yield.SmithWilson`](@ref)
- [`FinanceModels.Yield.NelsonSiegel`](@ref)
- [`FinanceModels.Yield.NelsonSiegelSvensson`](@ref)
- [`FinanceModels.Yield.CairnsPritchard`](@ref)
- [`FinanceModels.Yield.CairnsPritchardExtended`](@ref)
- [`FinanceModels.Yield.MonotoneConvex`](@ref) (Hagan-West)
- [`FinanceModels.Yield.ZeroRateCurve`](@ref) (direct construction from zero rates and tenors)

#### Curve Transformations

- `curve + curve` — additive composition of zero rates ([`FinanceModels.Yield.CompositeYield`](@ref))
- `curve * scalar` / `curve / scalar` — scale zero rates ([`FinanceModels.Yield.ScaledYield`](@ref))
- [`FinanceModels.Yield.TenorShift`](@ref) — lazy zero-rate shift depending on tenor: `(rate, tenor) -> Rate` (formerly `TransformedYield`)
- [`FinanceModels.Yield.ProjectedShift`](@ref) — lazy zero-rate shift depending on tenor *and* a projection time `τ`: `(τ, rate, tenor) -> Rate`
- [`FinanceModels.Yield.ForwardStarting`](@ref) — rebase a curve to a new time-zero
- [`FinanceModels.Yield.Blend`](@ref) — weighted blend of curves in discount-factor, zero-rate or forward-rate space, with constant weights or weights that vary with tenor (see [Blending curves](@ref))

### Available Models - Stochastic Short Rates

Stochastic short-rate models with closed-form zero-coupon bond prices. These are full yield models that also support Monte Carlo simulation. See the [Stochastic Interest Rate Models](@ref) guide for details and examples.

- [`FinanceModels.ShortRate.Vasicek`](@ref)
- [`FinanceModels.ShortRate.CoxIngersollRoss`](@ref)
- [`FinanceModels.ShortRate.HullWhite`](@ref)

#### Arithmetic

Adding, subtracting, multiplying, dividing, and comparing rates is supported.

Yield models can also be composed. Here is an example of fitting rates and spread separately and then adding the two models together:

```julia-repl
julia> q_rate = ZCBYield([0.01,0.02,0.03]);

julia> q_spread = ZCBYield([0.01,0.01,0.01]);

julia> model_rate = fit(Spline.Linear(),q_rate,Fit.Bootstrap());

julia> model_spread = fit(Spline.Linear(),q_spread,Fit.Bootstrap());

julia> forward(model_spread + model_rate,0,1)
Continuous(0.019900661706336083)

julia> forward(model_spread + model_rate,0,1) |> Periodic(1)
Periodic(0.020100000000000007, 1)

julia> discount(model_spread + model_rate,0,3)
0.8882274785263332

julia> discount(model_spread,0,3) * discount(model_rate,0,3)
0.8882274785263332
```

Note that combining two 1% (annual effective) curves gives a first-year rate of 2.01% effective, not 2%: curve addition compounds the two deflators, which retains the cross-term that simple rate addition drops. See [Yield Curve Arithmetic](@ref) for details.

!!! warning "Caution with Spreads"

    It is fairly common to see spreads and rates provided separately where both are quoted in par convention. For example, US Treasury par rates with the associated par risk spreads. Because par rates are dependent on the amount and path of rates preceeding the given tenor, **it is not valid to construct a "spread curve" with par rates and then use it in composition with a "rate curve"**.
    
    That is, while the zero rates and spreads in the preceeding example allow for additive or subtractive composition, it is not the case for par rates and spreads. Note the different discount factors produced:

    ```julia
    q_rate = ParYield([0.01,0.02,0.03]);
    q_spread = ParYield([0.01,0.01,0.01]);
    q_yield = ParYield([0.02,0.03,0.04]);

    model_rate = fit(Spline.Linear(),q_rate,Fit.Bootstrap());
    model_spread = fit(Spline.Linear(),q_spread,Fit.Bootstrap());
    model_yield = fit(Spline.Linear(),q_yield,Fit.Bootstrap());

    # The curves are different!
    discount(model_spread + model_rate,3)
    # 0.8864304594826338

    discount(model_yield,3)
    # 0.8864384275063881
    ```

    To apply par spreads over an existing curve without refitting it, use [`quote_shift`](@ref FinanceModels.Yield.quote_shift); it meets the shifted quotes at the quoted tenors only.

#### Yield Shifts: `TenorShift` and `ProjectedShift`

Two concrete subtypes of [`Yield.AbstractYieldShift`](@ref FinanceModels.Yield.AbstractYieldShift) apply lazy transformations to a base curve's zero rates without discretizing or refitting. They differ in whether the shift depends on one or two time axes.

##### `TenorShift` — shift depends on tenor only

[`Yield.TenorShift`](@ref FinanceModels.Yield.TenorShift) (formerly `TransformedYield`) applies a shift that may vary with the tenor `t`. The rule function receives the base curve's `Continuous` zero rate and the tenor, and returns a new rate:

```julia-repl
julia> base = Yield.Constant(0.05);  # 5% annual effective: continuous zero = log(1.05) ≈ 0.0488

julia> # Parallel shift of +100 bp annual effective: convert to annual effective, add, convert back
       shifted = base + (z, t) -> Continuous(Periodic(1)(z) + 0.01);

julia> zero(shifted, 10)  # = log(1.06)
Continuous(0.05826890812397578)

julia> # `z + Periodic(0.01, 1)` adds in `z`'s (continuous) convention, so it adds log(1.01):
       # the annual effective rate becomes 6.05%, not 6%
       zero(base + (z, t) -> z + Periodic(0.01, 1), 10)
Continuous(0.05874049502260009)

julia> # Continuous shift (simpler when convention is known)
       shifted2 = base + (z, t) -> z + Continuous(0.01);

julia> zero(shifted2, 10)  # = log(1.05) + 0.01
Continuous(0.05879016416943201)

julia> # Tenor-dependent twist (steepener that fades at 30y)
       twist = base + (z, t) -> z + Continuous(0.02 * max(0.0, 1.0 - t/30.0));

julia> zero(twist, 1)   # ≈ 4.88% + 1.93%
Continuous(0.06812349750276533)

julia> zero(twist, 30)  # shift is zero at 30y
Continuous(0.04879016416943201)
```

The rule function has the signature `(z::Rate, t) -> Rate`. The return value is type-asserted as `Rate`, so rules must carry compounding convention explicitly — returning a plain `Real` raises a `TypeError`. Because the rule receives the zero rate as a `Rate`, it can convert explicitly: `Continuous(Periodic(1)(z) + 0.01)` adds 100 bp annual effective. `Rate + Rate` adds in the left operand's convention.

The `+` operator provides ergonomic construction: `curve + f` or `f + curve` both create a `TenorShift`. This is distinct from `curve + scalar` (which creates a `CompositeYield`).

##### `ProjectedShift` — shift depends on tenor *and* a projection time

[`Yield.ProjectedShift`](@ref FinanceModels.Yield.ProjectedShift) extends `TenorShift` with a second time axis `τ`: the **projection time** (as-of / valuation-date offset) at which the curve is being evaluated. The rule signature is `(τ, z::Rate, t) -> Rate`, separating the projection time from the tenor.

This is the natural shape for scenario frameworks where a shift's magnitude evolves across a multi-year horizon — BMA SBA phase-ins, IFRS17 macroeconomic scenarios, embedded-value runoffs:

```julia-repl
julia> base = Yield.Constant(0.05);  # continuous zero = log(1.05) ≈ 0.0488

julia> # -150 bp parallel shift, phased in linearly over 10 projection years.
       phase_in = (τ, z, _) -> z + Continuous(-0.015 * min(τ, 10) / 10);

julia> # Curve as seen at projection year 3 (30% phased in → -45 bp).
       c3 = Yield.ProjectedShift(base, phase_in, 3.0);

julia> zero(c3, 5)
Continuous(0.04429016416943201)

julia> # Curve as seen at projection year 10 (fully phased in → -150 bp).
       c10 = Yield.ProjectedShift(base, phase_in, 10.0);

julia> zero(c10, 5)
Continuous(0.03379016416943201)
```

The intended pattern: store `phase_in` once as a first-class, year-independent value, then call `ProjectedShift(base, phase_in, τ)` at each projection time `τ` in a projection loop.

There is no `+` operator sugar for `ProjectedShift` — fixing `τ` at composition time defeats the purpose of storing the rule as a year-independent value. Always use the explicit constructor.

#### Quote shifts: `quote_shift`

`TenorShift` moves zero rates. To move *quoted* rates instead (par yields, CMT or OIS rates, swap rates), use [`quote_shift`](@ref FinanceModels.Yield.quote_shift). It takes a quote family and its tenors, applies a rule to each tenor's quote, and returns the curve plus continuous zero-rate bumps at those tenors that reprice every shifted quote. The rule receives the quote `q` that [`implied_quote`](@ref FinanceModels.Yield.implied_quote) gives at tenor `t`, and returns the new quote in the family's own units:

```julia-repl
julia> tenors = [1, 2, 3, 5, 7, 10, 20, 30];

julia> par2(y, t) = ParYield(y, t; frequency = 2);

julia> curve = fit(Spline.MonotoneConvex(), par2.([0.045, 0.042, 0.04, 0.039, 0.04, 0.042, 0.046, 0.045], tenors));

julia> up = quote_shift(curve, par2, tenors) do q, t
           q + 0.01   # every quoted par yield +100 bp
       end;

julia> implied_quote(curve, par2, 10.0), implied_quote(up, par2, 10.0)
(0.04199999997286451, 0.05199999997286452)

julia> implied_quote(curve, par2, 15.0), implied_quote(up, par2, 15.0)  # between tenors: about +100 bp
(0.04464519250910856, 0.054628543967210313)

julia> key_rate = quote_shift(curve, CMTYield, tenors) do q, t
           t == 10 ? q + 0.0001 : q   # 10-year CMT +1 bp
       end;

julia> implied_quote(key_rate, CMTYield, 7.0) - implied_quote(curve, CMTYield, 7.0)
0.0
```

The same call applies par spreads over an existing curve without refitting it, for example `quote_shift((q, t) -> q + spread[t], curve, par2, tenors)`.

The result is a [`CompositeYield`](@ref FinanceModels.Yield.CompositeYield): `curve + ZeroRateCurve(bumps, tenors, Spline.Linear(); extrapolation = :flat_zero)`. Each bump is a continuous zero-rate hat, linear from the previous tenor to the next and flat beyond the first and last tenors. What holds:

- **Exact at the tenors only.** Each shifted quote reprices within solver tolerance. Between tenors the quotes move approximately (above, +99.83 bp at 15 years), and forwards jump at the tenors.
- **No shape guarantees.** The bumps do not keep the base curve's properties, such as `MonotoneConvex`'s positive forwards.
- **Exact first-order derivatives**, through the base curve and the new quotes, as described in [Quote shifts and quote risk](@ref quote-shift-risk).
- **Quotes it accepts.** Each quote may pay only up to its own tenor and must pay something after the previous tenor. There, its amounts may change sign at most once, from negative to positive, and must end positive. Its price must exceed the value of its payments that the bump cannot move. Every standard family qualifies, including at negative rates.

To recalibrate instead, shifting the quotes that built the curve and fitting again with the same construction, call `fit` on the shifted quotes. That gives risk to the quotes passed to `fit`, and moves the curve between tenors as the construction dictates.

### Creating New Yield Models

See the [FinanceModels.jl Guide](@ref) for an example of creating a model from scratch. Some additional aspects to note:

- The only method that must be defined to calculate the [`FinanceCore.present_value`](@ref) of something is [`FinanceCore.discount`](@ref). Other methods will be inferred.
- For a curve that defines only `discount`, the interval factor `discount(curve, from, to)` is the ratio `discount(curve, to) / discount(curve, from)`. The built-in curves instead compute interval factors from their cumulative log-discount (for example z(t)·t from a zero rate), so their far-tail intervals stay finite where both discount factors underflow.
- Other methods that are imputed by default, but can be extended include: [`FinanceCore.accumulation`](@ref), [`FinanceModels.forward`](@ref), [`FinanceModels.par`](@ref), [`FinanceModels.zero`](@ref), and [`FinanceModels.rate`](@ref).

## Equity and Volatility Models

### Available Models - Option Valuation

- [`FinanceModels.Equity.BlackScholesMerton`](@ref)

### Available Models - Volatility

- [`FinanceModels.Volatility.Constant`](@ref)

#### Creating new Volatility Models

A volatility model must extend `volatility(vol::Volatility.MyNewModel, strike_ratio, time_to_maturity)`.

## Foreign Exchange Models

- [`FinanceModels.FX.Forwards`](@ref) — covered-interest-parity outright forward rates from a spot rate and a discount curve per currency; fit to [`FX.Outright`](@ref FinanceModels.FX.Outright) quotes.

See the [Foreign Exchange](@ref) guide page for pair conventions, curve construction, and cross-currency basis handling.
