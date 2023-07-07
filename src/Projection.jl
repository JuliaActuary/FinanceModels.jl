
abstract type AbstractProjection end

struct Projection{C,M,K} <: AbstractProjection
    contract::C
    model::M
    kind::K
end

## ProjectionKind ###############################
# controls what gets produced from the model,
# e.g. if you just want cashflows or you want full amortization schedule you might define an AmortizationSchedule kind which shows principle, interest, etc.

abstract type ProjectionKind end

struct CashflowProjection <: ProjectionKind end

# Collecting a Projection #######################
# Map(identity) is a Transducer, for which `collect` is defined. More on Transducers below

# collecting a Projection gives your the reducable defined below with __foldl__
Base.collect(p::P) where {P<:AbstractProjection} = p |> Map(identity) |> collect
# collecting a contract wraps the contract in with the default Proejction, defined next
Base.collect(c::C) where {C<:FinanceCore.AbstractContract} = Projection(c) |> Map(identity) |> collect

# Default Projections ##########################

# the default projection is just one where we get the cashflows and assume that the contract needs 
# no assumptions/model to determine the cashflows (the contract will error if a certain model is needed)
Projection(c) = Projection(c, NullModel(), CashflowProjection())
# if the model is also given, assume that we want a `CashflowProjection` by default
Projection(c, m) = Projection(c, m, CashflowProjection())


# Reducibles ###################################

# a more composible, efficient way to create a collection of things that you can apply subsequent transformations to 
# (and those transformations can be Transducers).
# https://juliafolds2.github.io/Transducers.jl/stable/howto/reducibles/
# https://www.youtube.com/watch?v=6mTbuzafcII


# There are two ways to define a reducible collection provided by Transducers.jl:
# `asfoldable` where you can define your reducible in terms of transducers
# `__foldl__` where you can define the collection using a `for` loop
# and `foldl__` you can also define state that is used within the loop

# this wraps a contract in a default proejction and makes a contract a reducible collection of cashflows
function Transducers.asfoldable(c::C) where {C<:FinanceCore.AbstractContract}
    Projection(c) |> Map(identity)
end

# A cashflow is the simplest, single item reducible collection
@inline function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:Cashflow,M,K}
    for i in 1:1
        val = @next(rf, val, p.contract)
    end
    return complete(rf, val)
end

# 
@inline function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:Bond.Fixed,M,K}
    b = p.contract
    ts = Bond.coupon_times(b)
    coup = b.coupon_rate / b.frequency.frequency
    for t in ts
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


# here a floating bond references the projections's model to determine 
# what the refernece rate is at that point in time
@inline function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:Bond.Floating,M,K}
    b = p.contract
    ts = Bond.coupon_times(b)
    for t in ts
        freq = b.frequency # e.g. `Periodic(2)`
        freq_scalar = freq.frequency  # the 2 from `Periodic(2)`

        # get the rate from the current time to next payment 
        # out of the model and convert it to the contract's periodicity
        model = p.model[b.key]
        reference_rate = rate(freq(forward(model, t, t + 1 / freq_scalar)))
        coup = (reference_rate + b.coupon_rate) / freq_scalar
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

# we simply concatenate two reducible collections to create a composite contract
@inline function Transducers.asfoldable(p::Projection{C,M,K}) where {C<:FinanceCore.Composite,M,K}
    # creates two sub-projections where the contract projected is decomposed to a non-composite contract
    # and then conctanate the two projections together
    ap = @set p.contract = p.contract.a
    bp = @set p.contract = p.contract.b
    (ap, bp) |> Cat()
end

# forward contract defines a set of cashflows that are relative to a future point in time,
# so we adjust the resulting cashflows `time`s by the forward start date
@inline function Transducers.asfoldable(p::Projection{C,M,K}) where {C<:Forward,M,K<:CashflowProjection}
    fwd_start = p.contract.time
    p_alt = @set p.contract = p.contract.instrument
    p_alt |> Map(cf -> @set cf.time += fwd_start)
end

@inline function Transducers.asfoldable(p::Projection{C,M,K}) where {C<:Cashflow,M,K<:CashflowProjection}
    Ref(p.contract) |> Map(identity)
end