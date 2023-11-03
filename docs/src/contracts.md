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
  - `Forward`to represent an instrument that is relative to a forward point in time.
  - `Composite` to represent the combination of two other instruments.  

In the future, this notion may be extended to liabilities (e.g. insurance policies in LifeContingencies.jl)


## `Cashflow` - a fundamental financial type

Say you wanted to model a contract that paid quarterly payments, and those payments occurred starting 15 days from the valuation date (first payment time = 15/365 = 0.057)

Previously, you had two options:
- Choose a discrete timestep to model (e.g. monthly, quarterly, annual) and then lump the cashflows into those timesteps. E.g. with monthly timesteps  of a unit payment of our contract, it might look like: `[1,0,0,1,0,0...]`
- Keep track of two vectors: one for the payment and one for the times. In this case, that might look like: `cfs = [1,1,...]; `times = `[0.057, 0.307...]`

The former has inaccuracies due to the simplified timing and logical complication related to mapping the contracts natural periodicity into an arbitrary modeling choice. The latter becomes unwieldy and fails to take advantage of Julia's type system. 

The new solution: `Cashflow`s. Our example above would become: `[Cashflow(1,0.057), Cashflow(1,0.307),...]`

### Creating a new Contract

A contract is anything that creates a vector of `Cashflow`s when `collect`ed. For example, let's create a bond which only pays down principle and offers no coupons.

```julia
using FinanceModels,FinanceCore

# Transducers is used to provide a more powerful, composible way to construct collections than the basic iteration interface
using Transducers: __foldl__, @next, complete

"""
A bond which pays down its par (one unit) in equal payments. 
"""
struct PrincipleOnlyBond{F<:FinanceCore.Frequency} <: FinanceModels.Bond.AbstractBond
    frequency::F
    maturity::Float64
end

# We extend the interface to say what should happen as the bond is projected
# There's two parts to customize:
# 1. any initialization or state to keep track of
# 2. The loop where we decide what gets returned at each timestep
function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:PrincipleOnlyBond,M,K}
    # initialization stuff
    b = p.contract # the contract within a projection
    ts = Bond.coupon_times(b) # works since it's a FinanceModels.Bond.AbstractBond with a frequency and maturity
    pmt = 1 / length(ts)

    for t in ts
        # the loop wich returns a value
        cf = Cashflow(pmt, t)
        val = @next(rf, val, cf) # the value to return is the last argument
    end
    return complete(rf, val)
end
```

That's it! then we can use this fitting models, projections, quotes, etc. Here we simply collect the bond into an array of cashflows:

```julia-repl
julia> PrincipleOnlyBond(Periodic(2),5.) |> collect
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

##### Cashflows are model dependent

An example of this is a floating bond where the coupon paid depends on a view of forward rates. See [this section in the overview](@ref Contracts-that-depend-on-the-model-(or-multiple-models)) on projections for how this is handled.

## Available Contracts & Modules

See the Modules in the left navigation for details on available contracts/models/functions.
