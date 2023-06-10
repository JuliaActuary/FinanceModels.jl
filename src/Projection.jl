
abstract type AbstractProjection end

struct Projection{M,K,C} <: AbstractProjection 
    model::M
    kind::K
    contract::C
end




# controls what gets produced from the model,
# e.g. if you just want cashflows or you want full amortization schedule, etc
abstract type ProjectionKind end

Base.collect(p::P) where {P<:AbstractProjection} = p |> Map(identity) |> collect

struct CashflowProjection <: ProjectionKind end

function Transducers.__foldl__(rf, val, p::Projection{M,K,C}) where {M,K,C<:Cashflow}
    for i in 1:1
        val = @next(rf, val, p.contract)
    end
    return complete(rf, val)
end


function Transducers.__foldl__(rf, val, p::Projection{M,K,C}) where {M,K,C<:Bond.Fixed}
    b = p.contract
    ts = Bond.coupon_times(b)
    for t in ts
        amt = if t == last(ts)
             1. + b.coupon_rate/b.frequency.frequency
             else
                b.coupon_rate/b.frequency.frequency
             end
        cf = Cashflow(amt,t)
        val = @next(rf, val, cf)
    end
    return complete(rf, val)
end

function Transducers.__foldl__(rf, val, p::Projection{M,K,C}) where {M,K,C<:Bond.Floating}
    b = p.contract
    ts = Bond.coupon_times(b)
    for t in ts
        reference_rate = Periodic(b.frequency.frequency)(rate(p.model[b.key],t))
        amt = if t == last(ts)
             1. + (rate(reference_rate)+ b.coupon_rate)/b.frequency.frequency
             else
                 (rate(reference_rate)+ b.coupon_rate)/b.frequency.frequency
             end
        cf = Cashflow(amt,t)
        val = @next(rf, val, cf)
    end
    return complete(rf, val)
end


function Transducers.asfoldable(p::Projection{M,K,C}) where {M,K,C<:Composite}
    ap = @set p.contract = p.contract.a
    bp = @set p.contract = p.contract.b
    a = ap |> Map(identity) |> collect 
    b = bp |> Map(identity) |> collect
    [a,b] |> Cat()
end