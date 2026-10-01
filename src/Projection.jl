abstract type AbstractProjection end

"""
    Projection(contract,model,kind)

The set of `contract`s and assumptions (`model`) to project the `kind` of output desired. Some assets require a projection in order to be valued (e.g. a floating rate bond).

If attempting to `collect` or otherwise reduce a contract (`<:AbstractContract`), by default it will get wrapped into a `Projection(contract,NullModel(),CashflowProjection())`

Use `Projection(contract; index)` to build a shared-index model store from
[`model_requirements`](@ref).
"""
struct Projection{C, M, K} <: AbstractProjection
    contract::C
    model::M
    kind::K
end

"""
    abstract type ProjectionKind

An abstract type that controls what gets produced from the model.

Subtypes of `ProjectionKind` define the level of detail in the output of the model. For example, if you just want cashflows or you want a full amortization schedule, you might define an `AmortizationSchedule` kind which shows principle, interest, etc.

After defining a new `ProjectionKind`, you need to define the how the projection works for that new output by extending either:

```julia
function Transducers.asfoldable(p::Projection{C,M,K}) where {C<:Cashflow,M,K<:CashflowProjection}
    ...
end
```
or 
```julia
function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:Cashflow,M,K<:CashflowProjection}
    ...
end
```

There are examples of this in the documentation.

# Examples
```julia
julia> struct CashflowProjection <: ProjectionKind end
CashflowProjection

julia> struct AmortizationSchedule <: ProjectionKind end
AmortizationSchedule
"""
abstract type ProjectionKind end

"""
    CashflowProjection()

A concrete subtype of `ProjectionKind` which is the projection which returns only a reducible collection of `Cashflow`s. Use in conjunction with a [`Projection`](@ref).
"""
struct CashflowProjection <: ProjectionKind end

# Collecting a Projection #######################
# Map(identity) is a Transducer, for which `collect` is defined. More on Transducers below

# collecting a Projection gives your the reducible defined below with __foldl__
Base.collect(p::P) where {P <: AbstractProjection} = p |> Map(identity) |> collect
# collecting a contract wraps the contract in with the default Projection, defined next
Base.collect(c::C) where {C <: FinanceCore.AbstractContract} = Projection(c) |> collect

# Default Projections ##########################

# the default projection is just one where we get the cashflows and assume that the contract needs
# no assumptions/model to determine the cashflows (the contract will error if a certain model is needed)
# The default and keyword index forms are defined in projection_models.jl.
# if the model is also given, assume that we want a `CashflowProjection` by default
Projection(c, m) = Projection(c, m, CashflowProjection())


# Reducibles ###################################

# a more composable, efficient way to create a collection of things that you can apply subsequent transformations to
# (and those transformations can be Transducers).
# https://juliafolds2.github.io/Transducers.jl/stable/howto/reducibles/
# https://www.youtube.com/watch?v=6mTbuzafcII


# There are two ways to define a reducible collection provided by Transducers.jl:
# `asfoldable` where you can define your reducible in terms of transducers
# `__foldl__` where you can define the collection using a `for` loop
# and `foldl__` you can also define state that is used within the loop

# this wraps a contract in a default projection and makes a contract a reducible collection of cashflows
function Transducers.asfoldable(c::C) where {C <: FinanceCore.AbstractContract}
    return Projection(c) |> Map(identity)
end
# Wrappers (a portfolio, `Composite`, `Forward`, `FX.Converted`) fold their children through
# `__foldl__` with the wrapper's own step composed into the reducing function, so they nest in any
# order: a fold that reaches a wrapper by `__foldl__` (as it does inside another wrapper) projects it
# just as a fold that starts there. A child is folded without completing, and the wrapper completes once.
# Each child goes through `asfoldable` first, as `Cat` does, so a contract defined by either extension
# point, `asfoldable` or `__foldl__`, folds inside any wrapper.
@inline function __fold_child(rf, val, child)
    rf0, coll = Transducers.retransform(rf, Transducers.asfoldable(child))
    return Transducers.foldl_nocomplete(rf0, val, coll)
end
# A portfolio is its members' projections, concatenated.
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: AbstractArray, M, K}
    members = map(c -> Projection(c, p.model, p.kind), p.contract)
    return Transducers.__foldl__(Transducers.Reduction(Cat(), rf), val, members)
end

# A cashflow is the simplest, single item reducible collection
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Cashflow, M, K}
    for i in 1:1
        val = @next(rf, val, p.contract)
    end
    return complete(rf, val)
end

# If a Transducer has been combined with a contract into an Eduction
# then unwrap the contract and apply the transducer to the projection
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Transducers.Eduction, M, K}
    rf = __rewrap(p.contract.rf, rf)             # compose the xform with any othe existing transducers
    p_alt = @set p.contract = p.contract.coll    # reset the contract to the underlying contract without transducers
    return Transducers.__foldl__(rf, val, p_alt)        # project with a newly combined reduction
end

#
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Bond.Fixed, M, K}
    b = p.contract
    ts = Bond.coupon_times(b)
    coup = b.coupon_rate / b.frequency.frequency
    # a maturity that is not a whole number of periods leaves a short first stub
    # (the schedule anchors at maturity and counts backward), which accrues its
    # actual length [0, t₁] rather than paying a full period's coupon
    first_coup = if Bond.__regular_schedule(b.maturity, b.frequency.frequency)
        coup
    else
        b.coupon_rate * first(ts)
    end
    for t in ts
        c = t == first(ts) ? first_coup : coup
        amt = if t == last(ts)
            1.0 + c
        else
            c
        end
        cf = Cashflow(amt, t)
        val = @next(rf, val, cf)
    end
    return complete(rf, val)
end


# here a floating bond references the projections's model to determine
# what the refernece rate is at that point in time
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Bond.Floating, M, K}
    b = p.contract
    ts = Bond.coupon_times(b)
    freq = b.frequency # e.g. `Periodic(2)`
    freq_scalar = freq.frequency  # the 2 from `Periodic(2)`
    model = p.model[b.key]
    regular = Bond.__regular_schedule(b.maturity, freq_scalar)
    for t in ts
        # Fix-in-advance: the coupon paid at t reflects the forward rate
        # observed at the START of the accrual period — [t - 1/freq, t] for a
        # whole period, [0, t₁] for the short first stub of a non-whole maturity.
        # This matches the standard market convention for FRNs and Ibor coupons,
        # and avoids referencing a rate beyond the bond's maturity for the final coupon.
        coup = if !regular && t == first(ts)
            # the stub accrues its actual length: the reference part is the
            # discount-factor ratio over the true window (never a lookback to
            # before issue), and the spread accrues simply over the stub
            (1 / discount(model, zero(t), t) - 1) + b.coupon_rate * t
        else
            reference_rate = rate(freq(forward(model, t - 1 / freq_scalar, t)))
            (reference_rate + b.coupon_rate) / freq_scalar
        end
        amt = if t == last(ts)
            1.0 + coup
        else
            coup
        end
        cf = Cashflow(amt, t)
        val = @next(rf, val, cf)
    end
    return complete(rf, val)
end

# A composite's cashflows are its first part's, then its second's
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: FinanceCore.Composite, M, K}
    val = __fold_child(rf, val, @set p.contract = p.contract.a)
    Transducers.@return_if_reduced val
    val = __fold_child(rf, val, @set p.contract = p.contract.b)
    Transducers.@return_if_reduced val
    return complete(rf, val)
end

# A forward contract's instrument runs on a clock that starts at the forward time: it observes its
# models and pays at times measured from there. Project it against the models seen from that start,
# then move its payments onto the projection's clock, so projecting `Forward(s, c)` is projecting `c`
# on the models seen from `s`, with every cashflow shifted by `s`.
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Forward, M, K <: CashflowProjection}
    s = p.contract.time
    inner = Projection(p.contract.instrument, __Rebased(p.model, s), p.kind)
    val = __fold_child(Transducers.Reduction(Map(cf -> @set cf.time += s), rf), val, inner)
    Transducers.@return_if_reduced val
    return complete(rf, val)
end

# The models a forward-starting instrument reads, seen from its start: each is rebased when the
# projection reads it by key, so what it never reads (a discount rate standing in for the model,
# an unused entry) passes through untouched, as does a store of any type that indexes by key.
struct __Rebased{M, T}
    models::M
    time::T
end
Base.getindex(r::__Rebased, key) = __rebase(r.models[key], r.time)
# A model seen from time `t`: a yield curve from `t` on; an FX model with its forward rate at `t` as
# spot and its curves from `t` on; a flat rate (a `Rate`, or a number read as one) is the same from any
# start.
__rebase(m::Yield.AbstractYieldModel, t) = Yield.ForwardStarting(m, t)
__rebase(m::Union{Real, FinanceCore.Rate}, t) = m
__rebase(m::FX.Forwards, t) = FX.Forwards(m.pair, forward(m, t), __rebase(m.domestic, t), __rebase(m.foreign, t))

# an `FX.Converted` contract's cashflows are denominated in the base (foreign) currency
# of an FX pair; convert each amount into the quote (domestic) currency at the
# arbitrage-free forward exchange rate for its payment time. The FX model is looked up
# from the projection's model store by key, mirroring `Bond.Floating`'s reference-rate
# lookup, so the wrapped contract still sees the same store (a converted floating leg
# resolves its own reference curve from it).
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: FX.Converted, M, K <: CashflowProjection}
    fx = p.model[p.contract.key]
    # a mis-keyed store (e.g. a yield curve where the FX model belongs) would otherwise
    # fail deep inside the transducer pipeline — `forward(curve, t)` returns a `Rate`,
    # which `Rate`-scales the amount and only errors at `Cashflow` reconstruction; name
    # the actual problem at its source instead
    if !(fx isa FX.AbstractFXModel)
        throw(ArgumentError("`FX.Converted` expects an FX model (e.g. `FX.Forwards`) under the key $(repr(p.contract.key)), but the projection's model store holds a $(typeof(fx))"))
    end
    # the declared pair guards the conversion's direction: an inverted (or crossed)
    # model under an otherwise-valid key would silently multiply by the reciprocal rate
    if fx.pair != p.contract.pair
        throw(ArgumentError("`FX.Converted` declared $(p.contract.pair) but the model under key $(repr(p.contract.key)) prices $(fx.pair)"))
    end
    inner = @set p.contract = p.contract.contract
    val = __fold_child(Transducers.Reduction(Map(cf -> @set cf.amount *= forward(fx, cf.time)), rf), val, inner)
    Transducers.@return_if_reduced val
    return complete(rf, val)
end

# an `FX.BasisSwapLeg` is a materialized strip of base-currency cashflows; emit them
# directly (they need no model to resolve). A contract needs `__foldl__`, not `asfoldable`:
# wrappers fold their children through `__foldl__`, which does not re-apply `asfoldable`.
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: FX.BasisSwapLeg, M, K}
    for cf in p.contract.cashflows
        val = @next(rf, val, cf)
    end
    return complete(rf, val)
end

@inline function Transducers.asfoldable(p::Projection{C, M, K}) where {C <: Cashflow, M, K <: CashflowProjection}
    return Ref(p.contract) |> Map(identity)
end

"""
    __rewrap(from::Transducers.Reduction, to)
    __rewrap(from, to)

Used to unwrap a Reduction which is a composition of contracts and a transducer and apply the transducers to the associated projection instead of the transducer.

For example, on its own a contract is not project-able, but wrapped in a (default) [`Projection`](@ref) it can be. But it may also be a lot more convienent 
to construct contracts which have scaling or negated modifications and let that flow into a projection.

# Examples

```julia-repl
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
"""
function __rewrap(from::Transducers.Reduction, to)
    rfx = from.xform                    # get the transducer's "xform" from the projection's contract
    return __rewrap(from.inner, Transducers.Reduction(rfx, to))          # compose the xform with any othe existing transducers
end
function __rewrap(from, to)
    # we've hit bottom, so return `to`
    return to
end

# The projection kind is dispatched internally, so a contract with a closed form can define
# `present_value(model, p::Projection{MyContract})` without an ambiguity with this method.
FinanceCore.present_value(model, p::FinanceModels.Projection, cur_time = 0.0) = __present_value(model, p, p.kind, cur_time)

function __present_value(model, p, ::CashflowProjection, cur_time)
    xf = p |> Filter(cf -> cf.time >= cur_time) |> Map(cf -> FinanceCore.discount(model, cur_time, cf.time) * cf.amount)
    # init: an empty fold (e.g. valuing past maturity) is worth 0, not an error
    return foldxl(+, xf; init = 0.0)
end
