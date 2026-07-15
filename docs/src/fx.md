# Foreign Exchange

The [`FX`](@ref FinanceModels.FX) module extends the contracts/models/`fit` architecture to
currencies. The design observation is that under covered interest parity (CIP), an FX
forward curve is a *ratio of discount factors* scaled by spot:

```math
F(t) = S_0 \cdot \frac{DF_{\text{foreign}}(t)}{DF_{\text{domestic}}(t)}
```

so FX needs no new curve mathematics — every spline, parametric model, bootstrap, and
piece of curve arithmetic in FinanceModels applies to FX curve construction unchanged.
What the module adds is the *semantic* layer where FX errors actually happen: explicit
pair direction, explicit points scale, and loud failures on mismatched currencies.

## Currency pairs

A currency pair is a type-level object. `FX.Pair(:EUR, :USD)` is the market's "EURUSD":
the price of one unit of the *base* currency (EUR) in units of the *quote* (domestic)
currency (USD). `inv` flips the direction:

```julia
eurusd = FX.Pair(:EUR, :USD)
inv(eurusd) # FX.Pair(:USD, :EUR)
```

Because the pair travels in the type, pricing a contract against a model quoted in a
different pair (including the inverted one) is an immediate `ArgumentError` naming both
pairs — not a silently crossed rate.

## The forward curve model

[`FX.Forwards`](@ref) holds the spot rate and the two discount curves. Any yield model
works for either curve:

```julia
using FinanceModels

usd = Yield.Constant(Continuous(0.05))
eur = Yield.Constant(Continuous(0.03))

m = FX.Forwards(eurusd, 1.10, usd, eur)

forward(m, 0.0) # 1.10 — the spot
forward(m, 1.0) # 1.10 * exp(0.05 - 0.03) ≈ 1.1222 — CIP forward
m(1.0)          # callable shorthand for forward(m, 1.0)

mi = inv(m)     # the USDEUR model: forward(mi, t) == 1 / forward(m, t)
```

An outright forward contract and its valuation (in the quote currency):

```julia
c = FX.Forward(eurusd, 1.1222, 1.0)  # receive 1 EUR, pay 1.1222 USD, at t = 1
present_value(m, c)                  # (forward(m, 1) - 1.1222) * discount(usd, 1)
```

## Building a curve from market quotes

FX forwards are quoted either as outrights or as points over spot. Both return zero-price
`Quote`s (entering a forward at the market rate costs nothing, just as a par bond quotes
at 1.0):

```julia
quotes = FX.Outright.(eurusd, [1.1055, 1.1113, 1.1225, 1.1459], [0.25, 0.5, 1.0, 2.0])

# or, from points (note the explicit pip scale — 10_000 for most pairs, 100 for JPY quotes):
quotes = FX.ForwardPoints.(eurusd, [55.0, 113.0, 225.0, 459.0], [0.25, 0.5, 1.0, 2.0]; spot = 1.10)
```

Given spot and the domestic curve, each forward quote pins the implied foreign discount
factor *in closed form*: ``DF_f(t) = F(t) \cdot DF_d(t) / S_0``. Curve construction is
therefore just a quote transformation followed by the ordinary fitting machinery, and the
one-step `fit` methods do exactly that:

```julia
m = fit(FX.Forwards(eurusd, 1.10, usd, Spline.Cubic()), quotes, Fit.Bootstrap())

forward(m, 1.0)                       # ≈ 1.1225: every quote reprices exactly
present_value(m, quotes[3].instrument) # ≈ 0.0
```

The intermediate transform is also available directly via
[`FX.implied_zcb_quotes`](@ref), which returns zero-coupon `Quote`s you can feed to any
model — a parametric curve, a different spline, `Fit.Loss` least squares, etc.

Parametric foreign curves fit through the generic optimizer path (the foreign curve's
parameters are the free variables; spot and the domestic curve are inputs):

```julia
m = fit(FX.Forwards(eurusd, 1.10, usd, Yield.Constant()), quotes)
```

## Cross-currency basis and hedge pricing

Since 2008, covered interest parity does not hold against "textbook" curves: market
forwards embed a cross-currency basis. The model handles this structurally, because the
`foreign` field is *defined by what reprices the FX market*, not as any particular
reference curve:

- **Absorbed basis.** When you `fit` an `FX.Forwards` to market forward quotes, the
  fitted `foreign` curve is the basis-adjusted (CSA / collateralized) discount curve —
  ``DF_f(t) = F_{\text{mkt}}(t) \cdot DF_d(t) / S_0`` by construction. Anything you then
  price with the model (forwards, hedged cashflow conversion) is basis-consistent
  automatically: the cost or benefit of the hedge is embedded in the curve.

- **Explicit basis.** Because curve arithmetic composes in zero-rate space, the basis is
  a first-class object when you want it to be:

  ```julia
  # from an OIS curve and a fitted basis spread curve:
  m = FX.Forwards(eurusd, 1.10, usd_ois, eur_ois + basis)

  # or extract the basis implied by the market-fit model:
  basis = m_fit.foreign - eur_ois
  ```

  A ``-20``bp EUR basis (`basis = Yield.Constant(Continuous(-0.002))`) lowers the
  effective EUR discounting rate and raises the CIP forwards, exactly as quoted markets
  do.

Long-dated basis calibration from cross-currency *swap* quotes (rather than forward
points) requires projecting floating legs in two currencies and is planned alongside the
multi-currency projection wrapper — see the roadmap note below.

## Relationship to other models

- **FX options**: the Garman–Kohlhagen model is `Equity.BlackScholesMerton(r_domestic,
  r_foreign, σ)` applied to one unit of base currency — the foreign rate plays the role
  of the dividend yield. Direct pricing methods for `FX` types and delta-conventioned
  volatility quotes are future work.
- **Multi-currency projection**: a wrapper contract that converts a foreign contract's
  projected cashflows at CIP forwards (enabling cross-currency swaps and hedged-portfolio
  valuation through the `Projection` model store) is the planned next step; `FX.Forwards`
  is the primitive it consumes.
