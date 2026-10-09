# Contracts

## **Contracts** - A composable way to represent financial instruments

Contracts are a composable way to represent financial instruments. They are, in essence, anything that is a collection of cashflows. Contracts can be combined to represent more complex instruments. For example, a bond can be represented as a collection of cashflows that correspond to the coupon payments and the principal repayment.

Examples:

- a `Cashflow`
- `Bond`s:
  - `Bond.Fixed`, `Bond.Floating`
- `Option`s:
  - `Option.EuroCall` and `Option.EuroPut`
- Compositional contracts:
  - `Forward` to represent an instrument that is relative to a forward point in time.
  - `Composite` to represent the combination of two other instruments.  

In the future, this notion may be extended to liabilities (e.g. insurance policies in LifeContingencies.jl)

## `Cashflow` - a fundamental financial type

Say you wanted to model a contract that paid quarterly payments, and those payments occurred starting 15 days from the valuation date (first payment time = 15/365 = 0.057)

Previously, you had two options:

- Choose a discrete timestep to model (e.g. monthly, quarterly, annual) and then lump the cashflows into those timesteps. E.g. with monthly timesteps  of a unit payment of our contract, it might look like: `[1,0,0,1,0,0...]`
- Keep track of two vectors: one for the payment and one for the times. In this case, that might look like: `cfs = [1,1,...];`times = `[0.057, 0.307...]`

The former has inaccuracies due to the simplified timing and logical complication related to mapping the contracts natural periodicity into an arbitrary modeling choice. The latter becomes unwieldy and fails to take advantage of Julia's type system.

The new solution: `Cashflow`s. Our example above would become: `[Cashflow(1,0.057), Cashflow(1,0.307),...]`

### Creating a new Contract

A contract is anything that creates a vector of `Cashflow`s when `collect`ed. For example, let's create a bond which only pays down principal and offers no coupons.

```julia
using FinanceModels,FinanceCore

# Transducers is used to provide a more powerful, composable way to construct collections than the basic iteration interface
import Transducers
using Transducers: __foldl__, @next, complete

"""
A bond which pays down its par (one unit) in equal payments. 
"""
struct PrincipalOnlyBond{F<:FinanceCore.Frequency} <: FinanceModels.Bond.AbstractBond
    frequency::F
    maturity::Float64
end

# We extend the interface to say what should happen as the bond is projected
# There's two parts to customize:
# 1. any initialization or state to keep track of
# 2. The loop where we decide what gets returned at each timestep
function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:PrincipalOnlyBond,M,K}
    # initialization stuff
    b = p.contract # the contract within a projection
    ts = Bond.coupon_times(b) # works since it's a FinanceModels.Bond.AbstractBond with a frequency and maturity
    pmt = 1 / length(ts)

    for t in ts
        # the loop which returns a value
        cf = Cashflow(pmt, t)
        val = @next(rf, val, cf) # the value to return is the last argument
    end
    return complete(rf, val)
end
```

That's it! We can now use this contract to fit models, create projections, quotes, etc. Here we simply collect the bond into an array of cashflows:

```julia-repl
julia> PrincipalOnlyBond(Periodic(2),5.) |> collect
10-element Vector{Cashflow{Float64, Float64}}:
 Cashflow{Float64, Float64}(0.1, 0.5)
 Cashflow{Float64, Float64}(0.1, 1.0)
 Cashflow{Float64, Float64}(0.1, 1.5)
 Cashflow{Float64, Float64}(0.1, 2.0)
 Cashflow{Float64, Float64}(0.1, 2.5)
 Cashflow{Float64, Float64}(0.1, 3.0)
 Cashflow{Float64, Float64}(0.1, 3.5)
 Cashflow{Float64, Float64}(0.1, 4.0)
 Cashflow{Float64, Float64}(0.1, 4.5)
 Cashflow{Float64, Float64}(0.1, 5.0)
```

Note that all contracts in FinanceModels.jl are currently *unit* contracts in that they assume a unit par value.

#### More complex Contracts

##### Sets of contracts

Sets of contracts can be put in an `AbstractArray` contained (e.g. a `Vector`) and then handled together. For example, we combine two bonds as a portfolio to project together:

```julia-repl
julia> c1 = Bond.Fixed(0.05, Periodic(1), 2.0);
julia> c2 = Bond.Fixed(0.04, Periodic(1), 2.0);

julia> Projection([c1, c2]) |> collect
4-element Vector{Cashflow{Float64, Float64}}:
 Cashflow{Float64, Float64}(0.05, 1.0)
 Cashflow{Float64, Float64}(1.05, 2.0)
 Cashflow{Float64, Float64}(0.04, 1.0)
 Cashflow{Float64, Float64}(1.04, 2.0)
```

##### Transformations

Contracts (`<:AbstractContract`) and [`Projection`](@ref)s can be modified to be scaled or transformed using the transformations in [Transducers.jl](https://juliafolds2.github.io/Transducers.jl/stable/#List-of-transducers) after importing that package.

Most commonly, this is likely simply chaining `Map(...)` calls. Two use-cases of this may be to (1) scale the contract by a factor or (2) change the sign of the contract to indicate an obligation/liability instead of an asset. Examples of this:

```julia-repl
julia> using Transducers, FinanceModels

julia> Bond.Fixed(0.05,Periodic(1),3) |> collect
3-element Vector{Cashflow{Float64, Float64}}:
 Cashflow{Float64, Float64}(0.05, 1.0)
 Cashflow{Float64, Float64}(0.05, 2.0)
 Cashflow{Float64, Float64}(1.05, 3.0)

julia> Bond.Fixed(0.05,Periodic(1),3) |> Map(-) |> collect
3-element Vector{Cashflow{Float64, Float64}}:
 Cashflow{Float64, Float64}(-0.05, 1.0)
 Cashflow{Float64, Float64}(-0.05, 2.0)
 Cashflow{Float64, Float64}(-1.05, 3.0)

julia> Bond.Fixed(0.05,Periodic(1),3) |> Map(-) |> Map(x->x*2) |> collect
3-element Vector{Cashflow{Float64, Float64}}:
 Cashflow{Float64, Float64}(-0.1, 1.0)
 Cashflow{Float64, Float64}(-0.1, 2.0)
 Cashflow{Float64, Float64}(-2.1, 3.0)
```

Another example of this is [`InterestRateSwap`](@ref). It's simply a `Composite` contract of a positive fixed rate bond and a negative floating rate bond. For a tenor of whole coupon periods, whose par coupon is the par yield, it is:

```julia
function my_swap(curve, tenor; frequency, model_key = "OIS")
    fixed_leg = Bond.Fixed(par(curve, tenor; frequency), frequency, tenor)
    float_leg = Bond.Floating(0.0, frequency, tenor, model_key) |> Map(-)
    return Composite(fixed_leg, float_leg)
end
```

`InterestRateSwap` also solves the par coupon of a tenor with a short first stub.

##### Cashflows are model dependent

A contract whose cashflows depend on other models, such as a floating bond's index curve, reads them
by key from the valuation context. [`Models`](@ref) holds the model that discounts and the models a
contract reads; `Models(model; index)` returns `index` for every key, and `Models(curve)` is its
single-curve case, `index = curve`. Value a contract under a
context, and list its cashflows by projecting it against the same context. For how a floating bond
projects its coupons from forward rates, see [this section in the overview](@ref Contracts-that-depend-on-the-model-(or-multiple-models)).

```julia
curve = Yield.Constant(0.04)
swap = InterestRateSwap(curve, 5.0; frequency = 1)
value(index, disc) = present_value(Models(disc; index), swap)
value(curve, curve) # approximately zero
collect(Projection(swap, Models(curve)))
```

Values add under one context: a `Composite` is worth the sum of its parts and a collection of contracts the
sum of its contracts' values. For multiple index curves or combined yield and FX models, use an
explicit store, for example `Models(ois, Dict("SOFR" => sofr, "EURUSD" => fx))`.

A contract with a closed-form value defines `present_value(ctx, contract)`, the one extension point
for a contract's value. It reads the models it needs from the context: `discount(ctx, t)` for
discounting, `ctx[key]` for an observed model, and [`valuation_model(ctx)`](@ref valuation_model)
for the model whose formula prices it. It then values inside a `Composite` or a portfolio:

```julia
struct OnePeriodFloater <: FinanceCore.AbstractContract
    key::String
end
# Pays principal plus the forward rate from t = 1 to 2 at t = 2.
FinanceCore.present_value(ctx, c::OnePeriodFloater) = discount(ctx, 2.0) / discount(ctx[c.key], 1.0, 2.0)

present_value(Models(ois, Dict("SOFR" => sofr)), [OnePeriodFloater("SOFR"), Bond.Fixed(0.04, Periodic(2), 5.0)])
```

When the formula depends on the type of the pricing model, `present_value` hands the context's
model to the model kernel [`FinanceModels.closed_form(model, contract)`](@ref FinanceModels.closed_form),
and each model that prices the contract adds a method for it. FinanceModels' options, caps, floors
and swaptions are valued this way. A custom contract opts in with one line:

```julia
# pays S_T - K at T on a unit stock
struct EquityForward <: FinanceCore.AbstractContract
    strike::Float64
    maturity::Float64
end
FinanceCore.present_value(ctx, c::EquityForward) = FinanceModels.closed_form(valuation_model(ctx), c)
FinanceModels.closed_form(m::Equity.BlackScholesMerton, c::EquityForward) =
    exp(-m.q * c.maturity) - c.strike * exp(-m.r * c.maturity)

m = Equity.BlackScholesMerton(0.03, 0.01, 0.2)
fwd = EquityForward(1.0, 2.0)
call = Option.EuroCall(CommonEquity(), 1.0, 2.0)
present_value(m, fwd)
present_value(Models(m; index = Yield.Constant(0.03)), fwd)   # the same value
present_value(m, [fwd, call])                                 # the sum of the two values
present_value(m, Composite(fwd, call))                        # the same sum
```

A model without a `closed_form` method for the contract throws a `MethodError`.

Wrappers that act on the cashflow stream (`Forward`, `FX.Converted`, `contract |> Map(f)`) need the
contract's projection, so a contract with only a closed form cannot be valued inside them.

## Quote conventions

The quote constructors encode market conventions. Pass rates in the convention of
the quoted instrument:

| Constructor | Instrument | Payment frequency | Rate convention |
|:--|:--|:--|:--|
| `ZCBYield(r, t)` | zero-coupon bond | at maturity | annual effective for a scalar `r`, or the `Rate`'s own |
| `ZCBPrice(p, t)` | zero-coupon bond | at maturity | price |
| `ParYield(r, t; frequency = 2)` | par bond | `frequency`, or a `Periodic` rate's own | nominal at that frequency |
| `CMTYield(r, t)` | US Treasury constant-maturity yield | at maturity for `t ≤ 1`; semiannual otherwise | annual effective for `t ≤ 1`; semiannual bond-equivalent otherwise |
| `OISYield(r, t)` | overnight index swap | at maturity for `t ≤ 1`; annual otherwise | annual |
| `ParSwapYield(r, t; frequency)` | par swap fixed leg | required `frequency` | nominal at that frequency |

FinanceModels measures time in year fractions, so day-count conventions are not
modeled.

## Available Contracts & Modules

See the Modules in the left navigation for details on available contracts/models/functions.
