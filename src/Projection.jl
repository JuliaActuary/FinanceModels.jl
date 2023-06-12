
abstract type AbstractProjection end

struct Projection{M,C,K} <: AbstractProjection
    model::M
    contract::C
    kind::K
    function Projection(model::M, contract::C, kind::K) where {M,C,K}
        new{M,C,K}(model, contract, kind)
    end
end




# controls what gets produced from the model,
# e.g. if you just want cashflows or you want full amortization schedule, etc
abstract type ProjectionKind end

Base.collect(p::P) where {P<:AbstractProjection} = p |> Map(identity) |> collect

struct CashflowProjection <: ProjectionKind end

function Transducers.__foldl__(rf, val, p::Projection{M,C,K}) where {M,C<:Cashflow,K}
    for i in 1:1
        val = @next(rf, val, p.contract)
    end
    return complete(rf, val)
end


function Transducers.__foldl__(rf, val, p::Projection{M,C,K}) where {M,C<:Bond.Fixed,K}
    b = p.contract
    ts = Bond.coupon_times(b)
    for t in ts
        amt = if t == last(ts)
            1.0 + b.coupon_rate / b.frequency.frequency
        else
            b.coupon_rate / b.frequency.frequency
        end
        cf = Cashflow(amt, t)
        val = @next(rf, val, cf)
    end
    return complete(rf, val)
end

function Transducers.__foldl__(rf, val, p::Projection{M,C,K}) where {M,C<:Bond.Floating,K}
    b = p.contract
    ts = Bond.coupon_times(b)
    for t in ts
        reference_rate = Periodic(b.frequency.frequency)(rate(p.model[b.key], t))
        amt = if t == last(ts)
            1.0 + (rate(reference_rate) + b.coupon_rate) / b.frequency.frequency
        else
            (rate(reference_rate) + b.coupon_rate) / b.frequency.frequency
        end
        cf = Cashflow(amt, t)
        val = @next(rf, val, cf)
    end
    return complete(rf, val)
end


function Transducers.asfoldable(p::Projection{M,C,K}) where {M,C<:Composite,K}
    ap = @set p.contract = p.contract.a
    bp = @set p.contract = p.contract.b
    (ap, bp) |> Cat()
end

function Transducers.asfoldable(p::Projection{M,C,K}) where {M,C<:Forward,K<:CashflowProjection}
    fwd_start = p.contract.time
    p_alt = @set p.contract = p.contract.instrument
    p_alt |> Map(cf -> @set cf.time += fwd_start)
end