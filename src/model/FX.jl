"""
The `FX` module provides foreign exchange models, contracts, and quote conventions:

- [`FX.Pair`](@ref) — a currency pair as a type-level object, e.g. `FX.Pair(:EUR, :USD)`
- [`FX.Forwards`](@ref) — a covered-interest-parity model of outright FX forward rates, built from a spot rate and two discount curves
- [`FX.Forward`](@ref) — an outright FX forward contract
- [`FX.Outright`](@ref) and [`FX.ForwardPoints`](@ref) — market quote conventions, returning `Quote`s
- [`FX.implied_zcb_quotes`](@ref) — transform FX forward quotes into the implied base-currency zero-coupon quotes used for curve construction

The design philosophy mirrors the rest of FinanceModels: market conventions that are a
common source of silent errors (pair direction, points scale, which curve discounts what)
are encoded explicitly in types and keyword arguments, and the curves inside an FX model
are ordinary yield models, so splines, parametric curves, curve arithmetic
(e.g. `foreign_ois + basis`), and `fit` all apply unchanged.

# Example

```julia
using FinanceModels

eurusd = FX.Pair(:EUR, :USD)
usd = Yield.Constant(Continuous(0.05))

# market outright forwards (or use FX.ForwardPoints for points quotes)
quotes = FX.Outright.(eurusd, [1.1055, 1.1113, 1.1225, 1.1459], [0.25, 0.5, 1.0, 2.0])

# closed-form implied EUR discount factors, bootstrapped through a spline
m = fit(FX.Forwards(eurusd, 1.10, usd, Spline.Cubic()), quotes, Fit.Bootstrap())

forward(m, 1.0)  # ≈ 1.1225, and every quote reprices to zero PV
```
"""
module FX

import ..AbstractModel
import ..FinanceCore: Timepoint
using ..FinanceCore

abstract type AbstractFXModel <: AbstractModel end

"""
    FX.Pair(base::Symbol, quote::Symbol)

A currency pair as a singleton type. `FX.Pair(:EUR, :USD)` denotes the price of one unit
of the *base* currency (`:EUR`) expressed in units of the *quote* — also called domestic
or terms — currency (`:USD`), matching the market's "EURUSD" naming.

Carrying the pair at the type level means that direction errors — the classic FX bug —
surface as immediate, descriptive errors rather than silently inverted prices: an
[`FX.Forwards`](@ref) model only prices [`FX.Forward`](@ref) contracts denominated in
the *same* pair.

`inv` flips the pair: `inv(FX.Pair(:EUR, :USD)) == FX.Pair(:USD, :EUR)`.

# Examples

```julia-repl
julia> FX.Pair(:EUR, :USD)
FinanceModels.FX.Pair{:EUR, :USD}()

julia> inv(FX.Pair(:EUR, :USD))
FinanceModels.FX.Pair{:USD, :EUR}()
```
"""
struct Pair{B, Q} end

Pair(base::Symbol, quote_currency::Symbol) = Pair{base, quote_currency}()

Base.inv(::Pair{B, Q}) where {B, Q} = Pair{Q, B}()
Base.Broadcast.broadcastable(p::Pair) = Ref(p)

"""
    FX.Forwards(pair, spot, domestic, foreign)

A covered-interest-parity (CIP) model of outright FX forward rates for the given
[`FX.Pair`](@ref):

    forward(m, t) = spot * discount(foreign, t) / discount(domestic, t)

where `forward(m, t)` is the arbitrage-free forward exchange rate (quote currency per one
unit of base currency) for delivery at time `t`, and `forward(m, 0) == spot`. The model is
also callable: `m(t) == forward(m, t)`.

# Arguments
- `pair`: the [`FX.Pair`](@ref) this model prices, e.g. `FX.Pair(:EUR, :USD)`.
- `spot`: the current exchange rate (quote currency per unit of base currency).
- `domestic`: the *quote*-currency discount curve (any yield model).
- `foreign`: the *base*-currency discount curve (any yield model). See the note on
  cross-currency basis below.

`inv(m)` returns the model for the flipped pair: spot is inverted and the curves swap
roles, so `forward(inv(m), t) == 1 / forward(m, t)`.

# Cross-currency basis

`foreign` is defined by what reprices the FX market, *not* necessarily the base currency's
own OIS or government curve. Since the collapse of covered interest parity against
textbook curves (post-2008), market forwards embed a cross-currency basis. Two equivalent
ways to handle it:

- **Absorbed**: `fit` the model to market forward quotes (see [`FX.implied_zcb_quotes`](@ref)
  and the `fit` methods for `FX.Forwards`). The fitted `foreign` curve *is* the
  basis-adjusted ("CSA" / collateralized) discount curve, and the basis is embedded —
  hedge pricing with the fitted model is basis-consistent by construction.
- **Explicit**: build `foreign` compositionally from a base-currency OIS curve and a
  fitted basis spread curve using ordinary curve arithmetic: `foreign_ois + basis`
  (zero-rate addition via `CompositeYield`). The basis curve is then a first-class,
  inspectable object: `basis = implied_foreign - foreign_ois`.

# Examples

```julia
eurusd = FX.Pair(:EUR, :USD)
usd = Yield.Constant(Continuous(0.05))
eur = Yield.Constant(Continuous(0.03))

m = FX.Forwards(eurusd, 1.10, usd, eur)

forward(m, 0.0)  # 1.10 (spot)
forward(m, 1.0)  # 1.10 * exp(0.05 - 0.03) ≈ 1.12222
m(1.0)           # same as forward(m, 1.0)

# explicit −20bp EUR/USD cross-currency basis on the EUR leg:
m_basis = FX.Forwards(eurusd, 1.10, usd, eur + Yield.Constant(Continuous(-0.002)))
```
"""
struct Forwards{P <: Pair, S, D, F} <: AbstractFXModel
    pair::P
    spot::S
    domestic::D
    foreign::F
end

FinanceCore.forward(m::Forwards, to) = m.spot * discount(m.foreign, to) / discount(m.domestic, to)
(m::Forwards)(t) = FinanceCore.forward(m, t)

Base.inv(m::Forwards) = Forwards(inv(m.pair), inv(m.spot), m.foreign, m.domestic)

"""
    FX.Forward(pair, strike, time)

An outright FX forward contract: at `time`, receive one unit of the base currency of
`pair` and pay `strike` units of the quote currency. Its value is expressed in the quote
currency.

Under an [`FX.Forwards`](@ref) model `m` for the same pair:

    present_value(m, c) = (forward(m, c.time) - c.strike) * discount(m.domestic, c.time)

A forward struck at the market outright has zero present value, which is how market
quotes are represented — see [`FX.Outright`](@ref).

Pricing a contract whose pair differs from the model's pair throws an `ArgumentError`
naming both pairs, rather than producing a silently inverted or crossed price.

# Examples

```julia
eurusd = FX.Pair(:EUR, :USD)
m = FX.Forwards(eurusd, 1.10, Yield.Constant(Continuous(0.05)), Yield.Constant(Continuous(0.03)))

atm = FX.Forward(eurusd, forward(m, 1.0), 1.0)
present_value(m, atm)  # 0.0

off = FX.Forward(eurusd, 1.10, 1.0)  # struck below the forward
present_value(m, off)  # (forward(m, 1.0) - 1.10) * discount at 1.0 > 0
```
"""
struct Forward{P <: Pair, K, T <: Timepoint} <: FinanceCore.AbstractContract
    pair::P
    strike::K
    time::T
end

FinanceCore.maturity(c::Forward) = c.time

function FinanceCore.present_value(m::Forwards{P}, c::Forward{P}) where {P <: Pair}
    return (FinanceCore.forward(m, c.time) - c.strike) * discount(m.domestic, c.time)
end

# a mismatched pair would otherwise fall through to the generic projection-based
# `present_value` and fail with an unrelated-looking iteration error; name the actual
# problem instead. The diagonal `{P,P}` method above is more specific, so matched pairs
# never reach this.
function FinanceCore.present_value(m::Forwards, c::Forward)
    throw(ArgumentError("cannot price an FX.Forward on $(c.pair) with an FX.Forwards model for $(m.pair)"))
end

"""
    FX.Outright(pair, forward_rate, time)

A `Quote` for an outright FX forward: entering a forward at the market rate costs
nothing, so this returns a zero-price `Quote` on an [`FX.Forward`](@ref) struck at
`forward_rate`. This is the analog of par-bond quoting for FX and is the input to the
`fit` methods for [`FX.Forwards`](@ref).

Use broadcasting to create a set of quotes: `FX.Outright.(pair, rates, times)`.

See also [`FX.ForwardPoints`](@ref) for the points-over-spot convention.

# Examples

```julia-repl
julia> FX.Outright(FX.Pair(:EUR, :USD), 1.1225, 1.0)
Quote{Float64, FinanceModels.FX.Forward{FinanceModels.FX.Pair{:EUR, :USD}, Float64, Float64}}(0.0, FinanceModels.FX.Forward{FinanceModels.FX.Pair{:EUR, :USD}, Float64, Float64}(FinanceModels.FX.Pair{:EUR, :USD}(), 1.1225, 1.0))
```
"""
Outright(pair::Pair, forward_rate, time) = Quote(zero(forward_rate), Forward(pair, forward_rate, time))

"""
    FX.ForwardPoints(pair, points, time; spot, scale=10_000)

A `Quote` for an FX forward quoted as *forward points* over spot, the interbank
convention: the outright forward rate is `spot + points / scale`.

`scale` is the pip factor and is deliberately explicit: it is `10_000` for most pairs
(one pip = 0.0001), but `100` for JPY-quoted pairs (one pip = 0.01). Returns the same
zero-price `Quote` as [`FX.Outright`](@ref).

# Examples

```julia
eurusd = FX.Pair(:EUR, :USD)
FX.ForwardPoints(eurusd, 25.0, 0.5; spot = 1.10)                      # outright 1.1025
FX.ForwardPoints(FX.Pair(:USD, :JPY), -30.0, 1.0; spot = 150.0, scale = 100)  # outright 149.70
```
"""
ForwardPoints(pair::Pair, points, time; spot, scale = 10_000) = Outright(pair, spot + points / scale, time)

"""
    FX.implied_zcb_quotes(quotes, spot, domestic)
    FX.implied_zcb_quotes(m::FX.Forwards, quotes)

Transform FX forward `quotes` (as produced by [`FX.Outright`](@ref) /
[`FX.ForwardPoints`](@ref)) into zero-coupon-bond `Quote`s for the *implied
base-currency discount curve*, given the `spot` rate and the `domestic`
(quote-currency) discount curve.

Because a forward's price is `(F(t) - K) * DF_d(t)`, the implied base-currency discount
factor at each quote's maturity is closed-form:

    DF_f(t) = (K * DF_d(t) + price) / spot

(for the usual zero-price outright quotes this is `F(t) * DF_d(t) / spot`). The returned
quotes can be fed to any curve-fitting method — this is how the `fit` methods for
[`FX.Forwards`](@ref) construct the foreign curve, and the resulting curve is the
basis-adjusted (CSA) discount curve described in [`FX.Forwards`](@ref).

All quotes must share a single currency pair (and match the model's pair in the second
form); mixed pairs throw an `ArgumentError`, since they would otherwise blend into a
meaningless curve.

# Examples

```julia
eurusd = FX.Pair(:EUR, :USD)
usd = Yield.Constant(Continuous(0.05))
quotes = FX.Outright.(eurusd, [1.1055, 1.1113, 1.1225], [0.25, 0.5, 1.0])

implied = FX.implied_zcb_quotes(quotes, 1.10, usd)   # Vector of ZCB-price Quotes
eur = fit(Spline.Cubic(), implied, Fit.Bootstrap())  # the implied EUR discount curve
```
"""
function implied_zcb_quotes(quotes, spot, domestic)
    pair = first(quotes).instrument.pair
    return map(quotes) do q
        c = q.instrument
        if c.pair != pair
            throw(ArgumentError("FX quotes must share a single currency pair; got both $(pair) and $(c.pair)"))
        end
        df = (c.strike * discount(domestic, c.time) + q.price) / spot
        Quote(df, Cashflow(one(df), c.time))
    end
end

function implied_zcb_quotes(m::Forwards, quotes)
    qpair = first(quotes).instrument.pair
    if qpair != m.pair
        throw(ArgumentError("quotes are for $(qpair) but the model prices $(m.pair)"))
    end
    return implied_zcb_quotes(quotes, m.spot, m.domestic)
end

end
