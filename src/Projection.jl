
abstract type AbstractProjection end

struct Projection{C,M,K} <: AbstractProjection
    contract::C
    model::M
    kind::K
end

Projection(c) = Projection(c, NullModel(), CashflowProjection())
Projection(c, m) = Projection(c, m, CashflowProjection())

function Transducers.asfoldable(c::C) where {C<:AbstractContract}
    Projection(c) |> Map(identity)
end
Base.collect(p::P) where {P<:AbstractProjection} = p |> Map(identity) |> collect
Base.collect(c::C) where {C<:AbstractContract} = Projection(c) |> Map(identity) |> collect


# controls what gets produced from the model,
# e.g. if you just want cashflows or you want full amortization schedule, etc
abstract type ProjectionKind end

struct CashflowProjection <: ProjectionKind end

function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:Cashflow,M,K}
    for i in 1:1
        val = @next(rf, val, p.contract)
    end
    return complete(rf, val)
end


function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:Bond.Fixed,M,K}
    b = p.contract
    ts = Bond.coupon_times(b)
    for t in ts
        coup = b.coupon_rate / b.frequency.frequency
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

function Transducers.__foldl__(rf, val, p::Projection{C,M,K}) where {C<:Bond.Floating,M,K}
    b = p.contract
    ts = Bond.coupon_times(b)
    for t in ts
        reference_rate = Periodic(b.frequency.frequency)(rate(p.model[b.key], t))
        coup = (forward(reference_rate, t, t + b.frequency.frequency) + b.coupon_rate) / b.frequency.frequency
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


function Transducers.asfoldable(p::Projection{C,M,K}) where {C<:Composite,M,K}
    ap = @set p.contract = p.contract.a
    bp = @set p.contract = p.contract.b
    (ap, bp) |> Cat()
end

function Transducers.asfoldable(p::Projection{C,M,K}) where {C<:Forward,M,K<:CashflowProjection}
    fwd_start = p.contract.time
    p_alt = @set p.contract = p.contract.instrument
    p_alt |> Map(cf -> @set cf.time += fwd_start)
end