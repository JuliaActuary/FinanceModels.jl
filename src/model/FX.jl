"""
The `FX` module provides foreign exchange models, contracts, and quote conventions:

- [`FX.Pair`](@ref) — a currency pair, e.g. `FX.Pair(:EUR, :USD)`
- [`FX.Forwards`](@ref) — a covered-interest-parity model of outright FX forward rates, built from a spot rate and two discount curves
- [`FX.Forward`](@ref) — an outright FX forward contract
- [`FX.Converted`](@ref) — a wrapper that converts a contract's projected cashflows into the quote currency at CIP forward rates, enabling cross-currency swaps and hedged multi-currency valuation
- [`FX.Outright`](@ref) and [`FX.ForwardPoints`](@ref) — market quote conventions, returning `Quote`s
- [`FX.ParBasisSwap`](@ref) — par cross-currency basis-swap spread quotes for calibrating the long-dated basis (with [`FX.BasisSwapLeg`](@ref) as the underlying instrument)
- [`FX.implied_zcb_quotes`](@ref) — transform FX forward quotes into the implied base-currency zero-coupon quotes used for curve construction

The design philosophy mirrors the rest of FinanceModels: market conventions that are a
common source of silent errors (pair direction, points scale, which curve discounts what)
are encoded explicitly in types and keyword arguments, and the curves inside an FX model
are ordinary yield models, so splines, parametric curves, curve arithmetic
(e.g. `foreign_ois + basis`), and `fit` all apply unchanged.

Present values follow one rule: `present_value` is denominated in the currency of the
contract's own cashflows — the *quote* currency for an [`FX.Forward`](@ref) (it settles
in quote-currency units), the *base* currency for an [`FX.BasisSwapLeg`](@ref) (a strip
of base-currency cashflows). Convert before combining values across denominations: a
base-currency present value times spot is the quote-currency present value.

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
import ..Bond
import ..FinanceCore: Timepoint
using ..FinanceCore

abstract type AbstractFXModel <: AbstractModel end

"""
    FX.Pair(base, quote)

A currency pair. `FX.Pair(:EUR, :USD)` denotes the price of one unit of the *base*
currency (`:EUR`) expressed in units of the *quote* — also called domestic or terms —
currency (`:USD`), matching the market's "EURUSD" naming.

Every FX contract and model carries its pair, so direction errors — the classic FX
bug — surface as immediate, descriptive `ArgumentError`s rather than silently inverted
prices: an [`FX.Forwards`](@ref) model only prices contracts denominated in an equal
pair.

`inv` flips the pair: `inv(FX.Pair(:EUR, :USD)) == FX.Pair(:USD, :EUR)`.

Currencies are usually `Symbol`s, but any values work — strings (including the
[InlineStrings](https://github.com/JuliaStrings/InlineStrings.jl) codes CSV readers
produce) or ISO 4217 numeric codes (`FX.Pair(978, 840)`). Pairs compare by *content*:
string-typed currencies are equal whenever their characters match, regardless of
storage type — `FX.Pair(String3("EUR"), String3("USD")) == FX.Pair("EUR", "USD")` —
because `AbstractString` equality and hashing are content-based. A `Symbol` pair and a
string pair are *not* equal, however; pick one convention for a given system,
normalizing at the data boundary (e.g. `Symbol.(codes)`) if sources disagree.

A degenerate pair such as `FX.Pair(:USD, :USD)` is deliberately permitted: with unit
spot and identical curves it is the *identity* pair (`forward ≡ 1`), which lets generic
multi-currency code route domestic contracts through the same [`FX.Converted`](@ref)
machinery as foreign ones.

# Examples

```julia-repl
julia> FX.Pair(:EUR, :USD)
FX.Pair(:EUR, :USD)

julia> inv(FX.Pair(:EUR, :USD))
FX.Pair(:USD, :EUR)
```
"""
struct Pair{B, Q}
    base::B
    quote_currency::Q # `quote` is a reserved word in Julia
end

Base.inv(p::Pair) = Pair(p.quote_currency, p.base)
# content-based equality with a matching hash: pairs match whenever their currencies
# compare equal, even across storage types — e.g. `String3("EUR")` vs `String7("EUR")`,
# which CSV readers legitimately assign to different columns of the same codes
Base.:(==)(a::Pair, b::Pair) = a.base == b.base && a.quote_currency == b.quote_currency
Base.hash(p::Pair, h::UInt) = hash(p.quote_currency, hash(p.base, hash(:FXPair, h)))
Base.Broadcast.broadcastable(p::Pair) = Ref(p)
Base.show(io::IO, p::Pair) = print(io, "FX.Pair(", repr(p.base), ", ", repr(p.quote_currency), ")")

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
roles, so `forward(inv(m), t) == 1 / forward(m, t)`. Inverting is calibration-exact —
a model fit to market quotes reprices the inverted quotes at `1 / F` identically, with
no refit.

There is deliberately no two-time `forward(m, from, to)` (a `MethodError`): under CIP
the fair forward for delivery at `to` is `forward(m, to)` no matter when the exchange
is contracted, so the yield-curve forward-forward form has no FX analog.

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

    present_value(m, c, cur_time = 0.0) = (forward(m, c.time) - c.strike) * discount(m.domestic, cur_time, c.time)

A forward struck at the market outright has zero present value, which is how market
quotes are represented — see [`FX.Outright`](@ref). `cur_time` moves only the
discounting date — the model is *not* rolled: `forward(m, c.time)` is still today's
forward, the package's usual static-curve valuation convention. A forward that settled
before `cur_time` is worth zero, consistent with how cashflows before `cur_time` are
treated everywhere else in the package.

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

function FinanceCore.present_value(m::Forwards, c::Forward, cur_time = 0.0)
    # a mismatched pair would otherwise price a silently crossed/inverted rate — the
    # classic FX bug this module exists to prevent
    if c.pair != m.pair
        throw(ArgumentError("cannot price an FX.Forward on $(c.pair) with an FX.Forwards model for $(m.pair)"))
    end
    v = (FinanceCore.forward(m, c.time) - c.strike) * discount(m.domestic, cur_time, c.time)
    # a forward settled before the valuation time contributes nothing, matching the
    # projection-path convention (`Filter(cf -> cf.time >= cur_time)` in Projection.jl)
    return c.time < cur_time ? zero(v) : v
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
julia> q = FX.Outright(FX.Pair(:EUR, :USD), 1.1225, 1.0);

julia> q.price, q.instrument.strike, q.instrument.time
(0.0, 1.1225, 1.0)
```
"""
Outright(pair::Pair, forward_rate, time) = Quote(zero(forward_rate), Forward(pair, forward_rate, time))

"""
    FX.ForwardPoints(pair, points, time; spot, scale)

A `Quote` for an FX forward quoted as *forward points* over spot, the interbank
convention: the outright forward rate is `spot + points / scale`.

`scale` is the pip factor and is required rather than defaulted: it is `10_000` for
most pairs (one pip = 0.0001) but `100` for JPY-quoted pairs (one pip = 0.01), and a
silently assumed scale is exactly the off-by-100× points error this module refuses to
guess about. Returns the same zero-price `Quote` as [`FX.Outright`](@ref).

# Examples

```julia
eurusd = FX.Pair(:EUR, :USD)
FX.ForwardPoints(eurusd, 25.0, 0.5; spot = 1.10, scale = 10_000)               # outright 1.1025
FX.ForwardPoints(FX.Pair(:USD, :JPY), -30.0, 1.0; spot = 150.0, scale = 100)  # outright 149.70
```
"""
ForwardPoints(pair::Pair, points, time; spot, scale) = Outright(pair, spot + points / scale, time)

# the closed-form implied base-currency discount factor for one forward quote:
# DF_f(t) = (K·DF_d(t) + price)/spot. Guarded because a corrupt quote (wrong price
# sign, points scale, or spot) implies df ≤ 0, which would otherwise surface only as
# NaNs inside a downstream curve fit, far from the offending quote.
function __implied_zcb_quote(q::Quote{<:Any, <:Forward}, spot, domestic)
    c = q.instrument
    df = (c.strike * discount(domestic, c.time) + q.price) / spot
    if !(df > 0)
        throw(ArgumentError("the FX forward quote at time $(c.time) implies a non-positive base-currency discount factor ($df); check the quote's price, points scale, and spot"))
    end
    return Quote(df, Cashflow(one(df), c.time))
end

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
meaningless curve. A quote implying a non-positive discount factor (a corrupt price,
points scale, or spot) also throws an `ArgumentError` naming the offending maturity,
rather than letting NaNs surface inside a downstream curve fit.

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
        if q.instrument.pair != pair
            throw(ArgumentError("FX quotes must share a single currency pair; got both $(pair) and $(q.instrument.pair)"))
        end
        __implied_zcb_quote(q, spot, domestic)
    end
end

function implied_zcb_quotes(m::Forwards, quotes)
    qpair = first(quotes).instrument.pair
    if qpair != m.pair
        throw(ArgumentError("quotes are for $(qpair) but the model prices $(m.pair)"))
    end
    return implied_zcb_quotes(quotes, m.spot, m.domestic)
end

"""
    FX.Converted(contract, key)

Wrap a contract whose cashflows are denominated in the *base* (foreign) currency of an
FX pair, converting each projected cashflow into the *quote* (domestic) currency at the
arbitrage-free forward exchange rate for that cashflow's time: an amount paid at time
`t` is multiplied by `forward(fx, t)`, where `fx` is the [`FX.Forwards`](@ref) model
looked up from the projection's model store under `key` — the same key/value pattern
used by [`Bond.Floating`](@ref FinanceModels.Bond.Floating) for its reference rate.

Converting at CIP forwards and then discounting on the domestic curve is identical to
discounting the unconverted cashflows on the FX model's foreign curve and converting the
result at spot:

    pv(domestic, Converted(c)) == spot * pv(foreign, c)

so collateralized (hedged) cross-currency valuation is consistent whichever way it is
sliced, and — when the FX model was fit to market forwards — automatically reflects the
cross-currency basis.

Because the wrapper works on *any* projectable contract (including transducer-modified
ones and [`Bond.Floating`](@ref FinanceModels.Bond.Floating), whose reference model is
resolved from the same store), cross-currency swaps are ordinary [`Composite`](@ref
FinanceCore.Composite)s — see the "Foreign Exchange" documentation page for a worked
fixed-for-fixed and floating-floating example.

`maturity(::Converted)` delegates to the wrapped contract, so it is undefined
(`MethodError`) when the wrapped contract is a transducer-modified `Eduction`: a
transducer may alter cashflow times, so no general delegation through it is sound.

# Examples

```julia
eurusd = FX.Pair(:EUR, :USD)
usd = Yield.Constant(Continuous(0.05))
eur = Yield.Constant(Continuous(0.03))
fx = FX.Forwards(eurusd, 1.08, usd, eur)

bond = Bond.Fixed(0.04, Periodic(1), 2.0)  # a EUR-denominated bond
p = Projection(FX.Converted(bond, "EURUSD"), Dict("EURUSD" => fx), CashflowProjection())

collect(p)  # cashflows now in USD: each amount is multiplied by forward(fx, t)
pv(usd, p)  # == 1.08 * pv(eur, bond)
```
"""
struct Converted{C, K} <: FinanceCore.AbstractContract
    # `contract` is deliberately unconstrained (like `FinanceCore.Composite`) so that
    # transducer-wrapped contracts (`Transducers.Eduction`s such as `bond |> Map(-)`)
    # can be converted too.
    contract::C
    key::K
end

FinanceCore.maturity(c::Converted) = FinanceCore.maturity(c.contract)

"""
    FX.BasisSwapLeg(pair, cashflows)

The foreign-currency leg of a constant-notional cross-currency basis swap, materialized
into deterministic cashflows: the periodic coupons — projection-curve forward plus the
quoted basis spread, times the accrual fraction — and the final unit principal. Amounts
are denominated in the *base* (foreign) currency of `pair`. Produced by
[`FX.ParBasisSwap`](@ref); see that docstring for the quoting assumptions.

Because the amounts are fixed once the projection curve is known, the leg is priceable
by any single discount curve, which is what lets par basis-swap quotes flow through the
standard curve-fitting machinery. Under an [`FX.Forwards`](@ref) model for the same
pair, `present_value(m, leg)` is the leg's value **in base-currency units**, priced on
`m.foreign` (the basis-adjusted/CSA discount curve); a mismatched pair throws an
`ArgumentError` naming both pairs.
"""
struct BasisSwapLeg{P <: Pair, C} <: FinanceCore.AbstractContract
    pair::P
    cashflows::C
end

FinanceCore.maturity(c::BasisSwapLeg) = last(c.cashflows).time

function FinanceCore.present_value(m::Forwards, c::BasisSwapLeg, cur_time = 0.0)
    if c.pair != m.pair
        throw(ArgumentError("cannot price an FX.BasisSwapLeg on $(c.pair) with an FX.Forwards model for $(m.pair)"))
    end
    return FinanceCore.present_value(m.foreign, c, cur_time)
end

"""
    FX.ParBasisSwap(pair, spread, maturity; reference, frequency=Periodic(4))

A `Quote` for a constant-notional cross-currency basis swap struck at par: receive the
foreign (base-currency) leg paying `reference`-curve forwards plus the quoted basis
`spread` (a decimal, e.g. `-0.0015` for −15bp) on a unit foreign notional, pay the
domestic (quote-currency) leg flat, with notionals exchanged at inception and
`maturity`.

Quoting assumption (standard for collateralized, OIS-discounted swaps): the domestic
leg pays the forwards of the *same* curve used for domestic discounting, so it is worth
par and drops out of the calibration. What remains is the base-currency par condition

    Σᵢ (fᵢ + spread) * δ * DF_f(tᵢ) + DF_f(T) = 1

where `fᵢ` are the (known) `reference` projection forwards, `δ` the accrual fraction,
and `DF_f` the basis-adjusted (CSA) discount curve being calibrated. The returned quote
is therefore `Quote(1.0, FX.BasisSwapLeg(pair, cashflows))`: a deterministic
base-currency cashflow strip that must price to par on the curve being fit. Each coupon
is fix-in-advance at the `reference` forward over its accrual period, exactly matching
how [`Bond.Floating`](@ref FinanceModels.Bond.Floating) projects, so the strip equals
the projected floating leg. A `maturity` that is not a whole number of periods
produces a front stub which — again exactly like `Bond.Floating` — still accrues
`1/frequency` and reads its reference forward from `t - 1/frequency`, extrapolating
the `reference` curve below time zero; prefer whole-period tenors.

Because the strip is deterministic, these quotes flow through the same `fit` methods as
[`FX.Outright`](@ref) quotes, and the two can be mixed in a single calibration —
forward points for the liquid short end, basis swaps for the long end:

```julia
eurusd = FX.Pair(:EUR, :USD)
sofr = Yield.Constant(0.05)  # domestic (USD) discount curve
estr = Yield.Constant(0.03)  # foreign (EUR) *projection* curve, already fitted

quotes = [
    FX.Outright.(eurusd, [1.1055, 1.1113], [0.25, 0.5]);
    FX.ParBasisSwap.(eurusd, [-0.0012, -0.0018], [2.0, 5.0]; reference = estr)
]
m = fit(FX.Forwards(eurusd, 1.10, sofr, Spline.Cubic()), quotes, Fit.Bootstrap())
```

This is the constant-notional representation. Interbank basis swaps are quoted on the
mark-to-market (resetting-notional) structure; under deterministic curves the two par
spreads agree to first order. The MTM contract itself is deliberately not modeled yet —
see the "Foreign Exchange" documentation page.
"""
function ParBasisSwap(pair::Pair, spread, maturity; reference, frequency = Periodic(4))
    f = frequency.frequency
    ts = Bond.coupon_times(maturity, f)
    cfs = map(ts) do t
        # fix-in-advance forward + spread, replicating `Bond.Floating`'s projection
        # (which the FX module cannot call directly: `Projection` loads after it)
        reference_rate = rate(frequency(forward(reference, t - 1 / f, t)))
        coup = (reference_rate + spread) / f
        amt = t == last(ts) ? 1.0 + coup : coup
        Cashflow(amt, t)
    end
    return Quote(1.0, BasisSwapLeg(pair, cfs))
end

# Reduce one market quote to its equivalent quote on the base-currency discount curve,
# using the spot and domestic curve carried by the model. This is the per-quote
# primitive behind the `fit` methods for `FX.Forwards`, and is what lets a single
# calibration mix instrument types: outright forwards for the short end, par basis
# swaps for the long end.
function __implied_foreign_quote(m::Forwards, q::Quote{<:Any, <:Forward})
    if q.instrument.pair != m.pair
        throw(ArgumentError("quote is for $(q.instrument.pair) but the model prices $(m.pair)"))
    end
    return __implied_zcb_quote(q, m.spot, m.domestic)
end

# a par basis-swap leg is already denominated on the base-currency curve; only check
# that it belongs to the model's pair
function __implied_foreign_quote(m::Forwards, q::Quote{<:Any, <:BasisSwapLeg})
    if q.instrument.pair != m.pair
        throw(ArgumentError("quote is for $(q.instrument.pair) but the model prices $(m.pair)"))
    end
    return q
end

end
