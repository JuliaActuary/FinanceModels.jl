# Foreign Exchange

The [`FX`](@ref FinanceModels.FX) module extends the contracts/models/`fit` architecture
to currencies. It rests on one observation: under covered interest parity (CIP), an FX
forward curve is a *ratio of discount factors* scaled by spot,

```math
F(t) = S_0 \cdot \frac{DF_{\text{foreign}}(t)}{DF_{\text{domestic}}(t)}
```

so FX requires no new curve mathematics — every spline, parametric model, bootstrap, and
piece of curve arithmetic in FinanceModels applies to FX curve construction unchanged.
What the module adds is the *semantic* layer where FX errors actually happen: explicit
pair direction, explicit points scale, explicit currency of every cashflow, and loud
failures on mismatches.

## Conventions: pairs, direction, and points

A currency pair is created with [`FX.Pair`](@ref). The first currency is the **base**
(the thing being priced — one unit of it), the second is the **quote** (also called
*domestic*, *terms*, or *counter*) currency (the units the price is expressed in). This
matches the market's concatenated naming:

| Pair                  | Market name | A rate of | means                |
|:----------------------|:------------|:----------|:---------------------|
| `FX.Pair(:EUR, :USD)` | EURUSD      | `1.10`    | \$1.10 per €1        |
| `FX.Pair(:USD, :JPY)` | USDJPY      | `150.0`   | ¥150 per \$1         |

`inv` flips the direction — `inv(FX.Pair(:EUR, :USD)) == FX.Pair(:USD, :EUR)` — and
correspondingly inverts the rate.

Currencies are usually `Symbol`s, but any values work — strings (plain or the
[InlineStrings](https://github.com/JuliaStrings/InlineStrings.jl) codes CSV readers
produce) or ISO 4217 numeric codes (`FX.Pair(978, 840)`). Pairs compare by *content*:
string-typed currencies match whenever their characters do, regardless of storage type
(`FX.Pair(String3("EUR"), String3("USD")) == FX.Pair("EUR", "USD")`, since
`AbstractString` equality and hashing are content-based). A `Symbol` pair and a string
pair are *not* equal, however — pick one convention per system, normalizing at the data
boundary (e.g. `Symbol.(codes)`) if sources disagree.

Direction confusion is the classic FX bug: an inverted rate is not obviously wrong on
sight (`0.91` and `1.10` are both plausible EUR/USD numbers), and a crossed pair can pass
through an entire calculation silently. Because every FX contract and model carries its
pair, pricing a contract against a model quoted in a different pair — including the
inverted one — throws an immediate `ArgumentError` naming both pairs rather than
returning a wrong number.

Forward *points* are the second convention trap. The interbank market quotes forwards as
pips over spot, and the pip size differs by pair: `0.0001` for most pairs but `0.01` for
JPY-quoted pairs. There is deliberately no points constructor: a function that guessed
the scale would be the classic silent JPY error, and one that required it would only
rename the arithmetic. Convert points explicitly at the data boundary, where the
pair-dependent pip factor stays visible at the call site:
`FX.Outright.(pair, spot .+ points ./ 10_000, times)`.

## Covered interest parity

Why must the forward rate equal ``S_0 \cdot DF_f(t)/DF_d(t)``? Compare two ways of
arranging to hold one unit of the base currency at time ``t``:

1. **Buy now and deposit**: buy ``DF_f(t)`` units of base currency today (cost
   ``S_0 \cdot DF_f(t)`` in quote currency) and invest them at the base-currency rate;
   they accumulate to exactly 1 unit at ``t``.
2. **Deposit and buy forward**: invest ``F(t) \cdot DF_d(t)`` of quote currency at the
   quote-currency rate (worth ``F(t)`` at ``t``) and contract today, at zero cost, to
   exchange it for 1 base unit at the forward rate ``F(t)``.

Both strategies deliver the same thing, so they must cost the same today:
``S_0 \, DF_f(t) = F(t) \, DF_d(t)``, which rearranges to the CIP formula. If the traded
forward deviates, the difference is a riskless profit ("covered interest arbitrage").

[`FX.Forwards`](@ref) is this equation as a model — a pair, a spot rate, and one discount
curve per currency, where each curve can be *any* yield model (constant, spline,
Nelson–Siegel, a `CompositeYield` sum, …):

```julia
using FinanceModels

eurusd = FX.Pair(:EUR, :USD)
usd = Yield.Constant(Continuous(0.05))
eur = Yield.Constant(Continuous(0.03))

m = FX.Forwards(eurusd, 1.10, usd, eur)

forward(m, 0.0) # 1.10 — the spot rate
forward(m, 1.0) # 1.10 * exp(0.05 - 0.03) ≈ 1.1222
m(1.0)          # callable shorthand for forward(m, 1.0)

mi = inv(m)     # the USDEUR model: forward(mi, t) == 1 / forward(m, t)
```

Note the direction of the effect: the *higher*-rate currency's forwards depreciate — here
USD rates exceed EUR rates, so EURUSD forwards sit *above* spot (EUR at a forward
premium). This is exactly Hull's introductory example (*Options, Futures, and Other
Derivatives*, §5.10): AUD/USD spot 0.6200 with 5% AUD and 7% USD rates gives a 2-year
forward of ``0.62\,e^{(0.07-0.05)\cdot 2} = 0.6453``, which is reproduced in this
package's test suite.

## Contracts and quotes

[`FX.Forward`](@ref) is the outright forward contract: at `time`, receive one unit of
base currency and pay `strike` units of quote currency. Its value, in quote currency, is
the discounted difference between the model forward and the strike:

```julia
c = FX.Forward(eurusd, 1.1222, 1.0)  # receive €1, pay $1.1222, at t = 1
present_value(m, c)                  # (forward(m, 1) - 1.1222) * discount(usd, 1)
```

A forward struck *at* the market rate costs nothing to enter. That is precisely how the
market quotes forwards, so FX quotes are zero-price `Quote`s — the FX analog of a par
bond quoting at 1.0:

```julia
quotes = FX.Outright.(eurusd, [1.1055, 1.1113, 1.1225, 1.1459], [0.25, 0.5, 1.0, 2.0])

# equivalently, from screen points over spot (pips: 1/10_000 here, 1/100 for JPY quotes):
points = [55.0, 113.0, 225.0, 459.0]
quotes = FX.Outright.(eurusd, 1.10 .+ points ./ 10_000, [0.25, 0.5, 1.0, 2.0])
```

## Constructing a curve from market quotes

Set the pricing equation of a quoted forward to its (zero) price and solve: with
``p = (F(t) - K)\,DF_d(t)`` and ``F(t) = S_0\,DF_f(t)/DF_d(t)``,

```math
DF_f(t) = \frac{K \cdot DF_d(t) + p}{S_0}
```

Given spot and the domestic curve, *each forward quote pins the implied foreign discount
factor in closed form* — no optimizer, no iteration. Curve construction therefore reduces
to a quote transformation followed by the ordinary fitting machinery, which is what the
one-step `fit` methods do:

```julia
m = fit(FX.Forwards(eurusd, 1.10, usd, Spline.Cubic()), quotes, Fit.Bootstrap())

forward(m, 1.0)                        # ≈ 1.1225: every quote reprices exactly
present_value(m, quotes[3].instrument) # ≈ 0.0
```

The intermediate transform is available directly as [`FX.implied_zcb_quotes`](@ref),
returning zero-coupon `Quote`s you can feed to *any* model or fitting method. Parametric
foreign curves calibrate through the generic optimizer path (the foreign curve's
parameters are the free variables; spot and the domestic curve are inputs):

```julia
m = fit(FX.Forwards(eurusd, 1.10, usd, Yield.Constant()), quotes)
```

## Cross-currency basis and hedge pricing

Textbook CIP stopped holding exactly in 2008: regulatory balance-sheet costs and demand
for dollar funding mean market forwards deviate from the CIP value computed off each
currency's own OIS curve. The deviation, expressed as a spread on one leg, is the
**cross-currency basis** (persistently negative vs. USD for EUR and especially JPY in the
post-crisis era).

The model handles this structurally, because the `foreign` field is *defined by what
reprices the FX market* rather than as any particular reference curve:

- **Absorbed basis.** When you `fit` an `FX.Forwards` to market forward quotes, the
  fitted `foreign` curve is the basis-adjusted (CSA / collateralized) discount curve —
  ``DF_f(t) = F_{\text{mkt}}(t) \cdot DF_d(t) / S_0`` by construction. Everything priced
  with the fitted model — forwards, hedge rolls, converted cashflows — is
  basis-consistent automatically: the cost or benefit of the currency hedge is embedded
  in the curve.

- **Explicit basis.** Because curve arithmetic composes in zero-rate space, the basis is
  a first-class curve when you want to see it:

  ```julia
  # from an OIS curve and a fitted basis spread curve:
  m = FX.Forwards(eurusd, 1.10, usd_ois, eur_ois + basis)

  # or extract the basis implied by a market-fit model:
  basis = m_fit.foreign - eur_ois
  ```

  A ``-20``bp EUR basis (`Yield.Constant(Continuous(-0.002))`) lowers effective EUR
  discounting and raises the CIP forwards, exactly as quoted markets do.

One practical consequence for hedged portfolios: the domestic-equivalent yield curve of a
fully-hedged foreign asset is `asset_curve - m.foreign + m.domestic` (all zero-rate
arithmetic) — the asset's spread over the basis-adjusted foreign curve, re-anchored to
domestic rates. The basis shows up as exactly the hedge drag it is in practice.

## Multi-currency valuation and cross-currency swaps

Valuing a foreign-currency contract in domestic terms is a *projection* concern: project
the contract's cashflows, convert each into the domestic currency, then discount as
usual. [`FX.Converted`](@ref) is that conversion as a wrapper contract. It multiplies
each projected amount by the CIP forward for that cashflow's time, looking the FX model
up from the projection's model store by key — the same key/value pattern
[`Bond.Floating`](@ref FinanceModels.Bond.Floating) uses for its reference rate:

```julia
store = Dict("EURUSD" => m, "ESTR" => eur)
eur_bond = Bond.Fixed(0.04, Periodic(1), 2.0)  # a EUR-denominated bond

p = Projection(FX.Converted(eur_bond, "EURUSD"), store, CashflowProjection())
collect(p)  # USD cashflows: each EUR amount × forward(m, t)
```

Converting at forwards and discounting domestically is *identical* to discounting on the
(basis-adjusted) foreign curve and converting at spot:

```math
\sum_t c_t \, F(t) \, DF_d(t) \;=\; S_0 \sum_t c_t \, DF_f(t)
```

so hedged cross-currency valuation is consistent whichever way it is sliced — this
identity is asserted directly in the test suite, for fixed and floating legs alike.

### A worked example: Hull's fixed-for-fixed currency swap

The classic currency-swap valuation (Hull, *Options, Futures, and Other Derivatives*,
§7.9): term structures are flat at 9% in the US and 4% in Japan (continuously
compounded), spot is ¥110 = \$1, and an institution receives 5% annually on ¥1,200M and
pays 8% annually on \$10M for three more years, with principals exchanged at maturity.
A currency swap is just two bonds — one converted:

```julia
jpyusd = FX.Pair(:JPY, :USD)
usd9 = Yield.Constant(Continuous(0.09))
jpy4 = Yield.Constant(Continuous(0.04))
fx = FX.Forwards(jpyusd, 1 / 110, usd9, jpy4)

swap = Composite(
    FX.Converted(Bond.Fixed(0.05, Periodic(1), 3) |> Map(cf -> cf * 1200.0), "JPYUSD"), # receive ¥ leg
    Bond.Fixed(0.08, Periodic(1), 3) |> Map(cf -> cf * -10.0),                          # pay $ leg
)

p = Projection(swap, Dict("JPYUSD" => fx), CashflowProjection())
pv(usd9, p) # ≈ 1.5430 ($ millions)
```

Hull values this swap both as a portfolio of forward contracts (his Table 7.9: forwards
0.009557, 0.010047, 0.010562 — compare `forward(fx, 1.0)` etc.) and as two bonds
converted at spot (``1{,}230.55/110 - 9.6439 = 1.543``). The projection above *is* the
forward-portfolio route; the two-bond route is
`1200 * pv(jpy4, yen_bond) / 110 - 10 * pv(usd9, usd_bond)`; both give \$1.543M, and the
test suite asserts they agree to numerical precision.

### Floating legs and basis swaps

Because a converted contract still sees the full model store, a converted
[`Bond.Floating`](@ref FinanceModels.Bond.Floating) resolves its own reference curve
while its cashflows are converted — so the market-standard constant-notional
floating-for-floating cross-currency basis swap is also just a `Composite`:

```julia
store = Dict("SOFR" => usd_ois, "ESTR" => eur_ois, "EURUSD" => m)

basis_swap = Composite(
    FX.Converted(Bond.Floating(-0.0020, Periodic(4), 5.0, "ESTR"), "EURUSD"), # receive €STR − 20bp, in USD terms
    Bond.Floating(0.0, Periodic(4), 5.0, "SOFR") |> Map(-),                   # pay SOFR flat
)
```

(The `Bond.Floating` legs are bond-style — they include the final principal — so the
principal exchange at maturity is represented; a par basis swap's initial exchange nets
against a unit `Cashflow` at time zero if you need it explicitly.)

### Calibrating the long-dated basis from par basis-swap quotes

Forward points cover the liquid short end (out to a year or two); beyond that the
cross-currency basis trades as **par basis-swap spreads** — "5y EUR/USD basis −18bp"
means the five-year swap of €STR − 18bp against SOFR flat, with notional exchanges,
prices to par. [`FX.ParBasisSwap`](@ref) turns those quotes into `fit` inputs.

The quoting assumption, standard for collateralized (OIS-discounted) swaps, is that the
domestic leg pays the forwards of the same curve used for domestic discounting, so it is
worth par and drops out. What remains is a condition on the foreign side alone:

```math
\sum_i (f_i + b)\,\delta\,DF_f(t_i) + DF_f(T) = 1
```

where ``f_i`` are forwards from the foreign *projection* curve (e.g. a fitted €STR
curve — a known input passed as the `reference` keyword), ``b`` is the quoted spread,
``\delta`` the accrual fraction, and ``DF_f`` the basis-adjusted discount curve being
calibrated. Because the coupon amounts are fully determined by the projection curve,
each quote materializes into a deterministic base-currency cashflow strip
([`FX.BasisSwapLeg`](@ref)) that must price to par on the curve being fit — and
deterministic strips flow through the same fitting machinery as the implied zero-coupon
quotes. One `fit` call takes both quote types:

```julia
estr = ...   # EUR projection curve, fitted from EUR OIS quotes
sofr = ...   # USD discount curve

quotes = [
    # short end: outrights from screen points over spot (pips = 1/10_000)
    FX.Outright.(eurusd, 1.10 .+ [25.0, 55.0, 120.0] ./ 10_000, [0.25, 0.5, 1.0]);
    # long end: par basis-swap spreads
    FX.ParBasisSwap.(eurusd, [-0.0012, -0.0015, -0.0018], [2.0, 5.0, 10.0]; reference = estr)
]

m = fit(FX.Forwards(eurusd, 1.10, sofr, Spline.Linear()), quotes, Fit.Bootstrap())
```

Every quote — outright or swap — reprices under the fitted model, and the composite
basis swap from the previous section values to zero against it (the test suite asserts
both). Two practical notes:

- Prefer *local* interpolation (e.g. `Spline.Linear`) when bootstrapping coupon-bearing
  quotes: a global spline reshapes earlier segments as later knots are added, drifting
  already-solved swaps off par.
- Interbank basis swaps are quoted on the mark-to-market (resetting-notional)
  structure. Under deterministic curves the constant-notional and MTM par spreads agree
  to first order, so quoted spreads calibrate the constant-notional representation
  directly; the difference only becomes material alongside FX–rates correlation, i.e.
  with a stochastic FX model.

Deliberately *not* built yet: the mark-to-market (resetting-notional) basis-swap
*contract* itself, pending convention decisions (which leg resets, settlement timing,
accrual on the adjustment). Open an issue if your use case needs it sooner.

Also deliberately deferred: **multi-pair consistency**. Each `FX.Forwards` is a
self-contained model of one pair, and independently fitted pair models are not forced
to be mutually consistent — nothing makes a directly fitted EURJPY model agree with the
EURJPY forwards implied by triangulating EURUSD × USDJPY, and each fit implies its own
version of any shared currency's discount curve. The scaling design is a per-currency
(*collateral-consistent*) market object: pick one numéraire/collateral currency, store
one spot rate and one collateral-consistent discount curve per currency, and *derive*
every pair's forwards from those. For ``N`` currencies that is ``N`` curves instead of
up to ``N(N-1)/2`` independent pair models, and triangle consistency holds by
construction rather than by calibration. The tradeoff is exact repricing: where the
market quotes a direct cross basis inconsistent with triangulation, such a model
reprices the cross at its derived arbitrage-consistent forward, not at the screen
quote. `FX.Forwards` is the two-currency case of that design and would become the
pairwise view a market object returns, so nothing in the current API forecloses it.

## Relationship to other models

- **FX options**: the Garman–Kohlhagen formula is `Equity.BlackScholesMerton(r_domestic,
  r_foreign, σ)` applied to one unit of base currency — the foreign rate plays the role
  of the dividend yield. Direct pricing methods for `FX` types and delta-conventioned
  volatility quotes (ATM/risk-reversal/butterfly) are future work.
- **Stochastic FX** (simulated spot paths over correlated short-rate models, for
  real-world/unhedged projection) would slot into the existing
  `simulate`/[`RatePath`](@ref FinanceModels.RatePath) machinery and is likewise future
  work.
