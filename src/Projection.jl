abstract type AbstractProjection end

"""
    Projection(contract,model,kind)

The set of `contract`s and assumptions (`model`) to project the `kind` of output desired. Some assets require a projection in order to be valued (e.g. a floating rate bond).

If attempting to `collect` or otherwise reduce a contract (`<:AbstractContract`), by default it will get wrapped into a `Projection(contract,NullModel(),CashflowProjection())`

`Projection(contract, models)` projects against a valuation context such as
[`Models`](@ref), whose models the contract reads by key: `collect(Projection(swap, Models(curve; index = curve)))`.
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

After defining a new `ProjectionKind`, you need to define how the projection works for that new output by extending either:

```julia
function Transducers.asfoldable(p::Projection{C,M,K}) where {C<:MyContract,M,K<:MyKind}
    ...
end
```
or
```julia
function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:MyContract,M,K<:MyKind}
    ...
end
```
Either extension point folds inside every wrapper (a portfolio, `Composite`, `Forward`,
`FX.Converted`, a transducer over the contract).

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

# Collecting folds the projection into a vector that widens to the types it receives. A projection
# that emits nothing gives `Union{}[]`: it is not iterable, so no element type can be inferred.
Base.collect(p::P) where {P <: AbstractProjection} = foldxl(push!!, p; init = Union{}[])
# collecting a contract wraps the contract in with the default Projection, defined next
Base.collect(c::C) where {C <: FinanceCore.AbstractContract} = Projection(c) |> collect
# Transducers over a projection or a contract (`Projection(c) |> Take(2)`) collect the same way.
# Transducers' own `collect` infers the element type of an empty result from an iterable, which
# neither is.
Base.collect(e::Transducers.Eduction{<:Any, <:Union{AbstractProjection, FinanceCore.AbstractContract}}) =
    foldxl(push!!, e; init = Union{}[])

# Default Projections ##########################

# A projection defaults to cashflows. Without a model it lists fixed cashflows; a contract that needs
# a model errors.
Projection(c, m) = Projection(c, m, CashflowProjection())
Projection(c) = Projection(c, NullModel())


# Reducibles ###################################

# a more composable, efficient way to create a collection of things that you can apply subsequent transformations to
# (and those transformations can be Transducers).
# https://juliafolds2.github.io/Transducers.jl/stable/howto/reducibles/
# https://www.youtube.com/watch?v=6mTbuzafcII


# Transducers.jl defines a reducible in one of two ways: `asfoldable`, in terms of transducers, or
# `__foldl__`, with a loop that can carry state.

# A bare contract folds as its default projection, so `collect(contract)` and a transducer over a
# contract reach the same `Projection` methods, whichever extension point the contract defines.
Transducers.asfoldable(c::FinanceCore.AbstractContract) = Transducers.asfoldable(Projection(c))

# A wrapper (a portfolio, `Composite`, `Forward`, `FX.Converted`) folds its children's projections
# through `Cat`, with its own step composed into `rf`. `Cat` applies `asfoldable` to each child, so
# either extension point folds inside any wrapper, and wrappers nest in any order.
@inline __concat(rf, val, children) = Transducers.__foldl__(Transducers.Reduction(Cat(), rf), val, children)

# A portfolio is its members' projections, concatenated.
@inline Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: AbstractArray, M, K} =
    __concat(rf, val, map(c -> Projection(c, p.model, p.kind), p.contract))

# A cashflow is the simplest, single item reducible collection
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Cashflow, M, K}
    val = @next(rf, val, p.contract)
    return complete(rf, val)
end

# A transducer over a contract (`contract |> xf`, an `Eduction`): the contract's projection, through
# `xf` in the order written, then into `rf`. `xf` keeps its own state (`Take`, `Scan`), started and
# completed here without restarting `rf`, whose input it reaches through `__Emit`. An early stop of
# `xf` (an inner `Take`) ends this stream only; an early stop of `rf` ends the whole fold, so `__Emit`
# marks it (`__OuterStop`) on its way back out through `xf`.
struct __OuterStop{A}
    acc::A
end
struct __Emit{F}
    rf::F
end
@inline function (e::__Emit)(acc, x)
    r = Transducers.next(e.rf, acc, x)
    return r isa Transducers.Reduced ? Transducers.reduced(__OuterStop(Transducers.unreduced(r))) : r
end
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Transducers.Eduction, M, K}
    xf = Transducers.Reduction(Transducers.Transducer(p.contract), Transducers.BottomRF(Transducers.Completing(__Emit(rf))))
    child = Projection(p.contract.coll, p.model, p.kind)
    acc = Transducers.foldl_nocomplete(Transducers.Reduction(Cat(), xf), Transducers.start(xf, val), (child,))
    # A stop arrives with `xf` already completed. Otherwise `xf` completes here, and may still emit
    # (a flushed `Partition`, the last group of `PartitionBy`) and so stop `rf`.
    acc isa Transducers.Reduced || (acc = complete(xf, acc))
    # Either way the result is `rf`'s: stopped and completed by `rf` itself, or to complete here.
    u = Transducers.unreduced(acc)
    return u isa __OuterStop ? Transducers.reduced(u.acc) : complete(rf, u)
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
# what the reference rate is at that point in time
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
@inline Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: FinanceCore.Composite, M, K} =
    __concat(rf, val, (Projection(p.contract.a, p.model, p.kind), Projection(p.contract.b, p.model, p.kind)))

# `Forward(s, c)` projects `c` on the models seen from `s`, then shifts every cashflow by `s`.
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Forward, M, K <: CashflowProjection}
    s = p.contract.time
    inner = Projection(p.contract.instrument, __Rebased(p.model, s), p.kind)
    return __concat(Transducers.Reduction(Map(cf -> @set cf.time += s), rf), val, (inner,))
end

# The models seen from `time`. A model is rebased only when the instrument reads it (by key, or as the
# valuation model), so unread entries, a discounting rate and any keyed store pass through.
struct __Rebased{M, T}
    models::M
    time::T
end
Base.getindex(r::__Rebased, key) = __rebase(r.models[key], r.time)
valuation_model(r::__Rebased) = __rebase(valuation_model(r.models), r.time)
# A model seen from time `t`: a yield curve from `t` on; an FX model with its forward rate at `t` as
# spot and its curves from `t` on; a flat rate (a `Rate`, or a number read as one) is the same from any
# start.
__rebase(m::Yield.AbstractYieldModel, t) = Yield.ForwardStarting(m, t)
__rebase(m::Union{Real, FinanceCore.Rate}, t) = m
__rebase(m::FX.Forwards, t) = FX.Forwards(m.pair, forward(m, t), __rebase(m.domestic, t), __rebase(m.foreign, t))

# `FX.Converted` multiplies each base-currency amount by the forward for its payment time. It reads
# its FX model from the store by key; its contract reads the same store, declared in the base currency
# (`__InBase`), so a converted floating leg finds its own reference curve.
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: FX.Converted, M, K <: CashflowProjection}
    fx = p.model[p.contract.key]
    # the declared pair guards the conversion's direction: an inverted (or crossed)
    # model under an otherwise-valid key would silently multiply by the reciprocal rate
    if fx.pair != p.contract.pair
        throw(ArgumentError("`FX.Converted` declared $(p.contract.pair) but the model under key $(repr(p.contract.key)) prices $(fx.pair)"))
    end
    # a context that reports in its own currency admits it; an enclosing conversion must convert from it
    q = p.contract.pair.quote_currency
    d = __denomination(p.model)
    d === nothing || __admits(d, q) || FX.__throw_misconverted(p.contract, q, d)
    inner = Projection(p.contract.contract, __InBase(p.model, p.contract.pair.base), p.kind)
    return __concat(Transducers.Reduction(Map(cf -> @set cf.amount *= forward(fx, cf.time)), rf), val, (inner,))
end

# The currency a projection's cashflows are declared in: any, under a curve or a listing (`NullModel`);
# an enclosing `FX.Converted`'s base currency (`__InBase`), kept through `Forward`; `nothing` under any
# other context, which reports in its own currency.
struct __AnyCurrency end
struct __InBase{S, B}
    store::S
    base::B
end
Base.getindex(m::__InBase, key) = m.store[key]
__denomination(::Union{Yield.AbstractYieldModel, NullModel}) = __AnyCurrency()
__denomination(m::__InBase) = m.base
__denomination(m::__Rebased) = __denomination(m.models)
__denomination(_) = nothing
__admits(d, currency) = d isa __AnyCurrency || isequal(d, currency)

# An `FX.Forward` is one quote-currency cashflow, priced on the context's own FX model (there is none
# inside an `FX.Converted`).
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: FX.Forward, M, K <: CashflowProjection}
    val = @next(rf, val, FX.__forward_cashflow(valuation_model(p.model), p.contract))
    return complete(rf, val)
end

# An `FX.BasisSwapLeg` is a strip of base-currency cashflows, projected only where that currency is
# admitted: under a curve (its native route) or inside an `FX.Converted` from it.
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: FX.BasisSwapLeg, M, K}
    __admits(__denomination(p.model), p.contract.pair.base) || FX.__throw_unconverted(p.contract)
    return Transducers.__foldl__(rf, val, p.contract.cashflows)
end


"""
    present_value(ctx, contract)

The value of `contract` at time 0 in the valuation context `ctx`: a yield curve or rate that
discounts its cashflows, a model with a closed form for it (`Equity.BlackScholesMerton` for a
European option), or a [`Models`](@ref) context that also holds the models it reads by key.

- A contract with a cashflow projection is worth its projected cashflows, each discounted by
  `ctx` (`discount(ctx, t)`). A cashflow paid before time 0 accumulates.
- A `Composite` is worth the sum of its parts, and a collection of contracts the sum of its
  contracts' values, so each part keeps its own valuation (a closed form or its projection).
- `Forward`, `FX.Converted` and transducers (`contract |> Map(f)`) act on the cashflow stream: a
  contract under one of them is valued by its projection.

There is no valuation-time argument. For a deterministic curve, the value as of `t` of the
cashflows at or after `t` discounts each from `t`:

```julia
foldxl(+, Projection(c, curve) |> Filter(cf -> cf.time >= t) |> Map(cf -> cf.amount * discount(curve, t, cf.time)); init = zero(discount(curve, t, t)))
```

# Examples

```julia
m = Equity.BlackScholesMerton(0.01, 0.02, 0.15)
a = Option.EuroCall(CommonEquity(), 1.0, 1.0)
pv(m, a) # ≈ 0.05410094201902403
```
"""
FinanceCore.present_value(ctx, c::FinanceCore.AbstractContract) = __stream_value(ctx, c)
# a transducer over a contract (a swap leg is one) is an ordered cashflow stream, as a wrapper is
FinanceCore.present_value(ctx, c::Transducers.Eduction{<:Any, <:FinanceCore.AbstractContract}) = __stream_value(ctx, c)

# The sum of the discounted cashflows, from a positive zero in the type of a present value under
# `ctx` (FinanceCore's empty-collection value), so the accumulator is concretely typed. The zero time
# is an `Int`: a DecFP curve rejects a `Bool` time.
__stream_value(ctx, c) = foldxl(
    +, Projection(c, ctx, CashflowProjection()) |> Map(cf -> FinanceCore.present_value(ctx, cf));
    init = zero(FinanceCore.present_value(ctx, false, 0))
)
