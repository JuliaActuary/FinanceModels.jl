
abstract type AbstractProjection end

struct Projection{M,K,C} <: AbstractProjection 
    model::M
    kind::K
    contract::C
end




# controls what gets produced from the model,
# e.g. if you just want cashflows or you want full amortization schedule, etc
abstract type ProjectionKind end

struct CashflowProjection <: ProjectionKind end
Base.eltype(p::Projection) = eltype(p.kind,p.contract)


# TODO this needs to be smarter about what the output type is for the flexibility in 
# the amount and Timepoint type
Base.eltype(pk::CashflowProjection,c::AbstractContract) = Cashflow{Float64,Float64}

function Base.iterate(p::Projection{M,K,C},state=1) where {M,K,C<:Cashflow}
    if state > 1
        return nothing
    else
        return (p.contract,state+1)
    end
end

function Base.getindex(p::Projection{M,K,C},i::Int) where {M,K,C<:Cashflow}
    if i > 1
        throw(BoundsError(p,i))
    else
        return p.contract
    end
end

Base.length(p::Projection{M,K,C}) where {M,K,C<:Cashflow}= 1

function Base.iterate(p::Projection{M,K,C},state=(Bond.coupon_times(p.contract),1)) where {M,K,C<:Bond.Fixed}
    b = p.contract
    if state[2] > lastindex(state[1])
        return nothing
    elseif state[2] == lastindex(state[1])
        return (Cashflow(1. + b.coupon_rate/b.frequency.frequency,state[1][end]), (state[1],state[2]+1))
    else
        return (Cashflow(b.coupon_rate/b.frequency.frequency,state[1][state[2]]), (state[1],state[2] + 1))
    end
end

function Base.getindex(p::Projection{M,K,C},i::Int) where {M,K,C<:Bond.Fixed}
    ct = Bond.coupon_times(p.contract)
    l  = length(ct)
    b = p.contract
    if i > l
        throw(BoundsError(p,i))
    elseif i == l
        return Cashflow(1. + b.coupon_rate/b.frequency.frequency,ct[i])
    else
        return Cashflow(b.coupon_rate/b.frequency.frequency,ct[i])
    end
end

Base.length(p::Projection{M,K,C}) where {M,K,C<:Bond.AbstractBond}  = length(Bond.coupon_times(p.contract))

function Base.iterate(p::Projection{M,K,C},state=(Bond.coupon_times(p.contract),1)) where {M,K,C<:Bond.Floating}
    b = p.contract
    if state[2] > lastindex(state[1])
        return nothing
    elseif state[2] == lastindex(state[1])
        # last period with principle payment
        reference_rate = Periodic(b.frequency.frequency)(rate(p.model[b.key],state[1][state[2]]))
        return (
            Cashflow(1. + (rate(reference_rate)+ b.coupon_rate)/b.frequency.frequency, state[1][end]),
            (state[1],state[2]+1)
            )
    else
        reference_rate = Periodic(b.frequency.frequency)(rate(p.model[b.key],state[1][state[2]]))
        return (
            Cashflow((rate(reference_rate)+ b.coupon_rate)/b.frequency.frequency, state[2]),
            (state[1],state[2] + 1)
            )
    end
end

function Base.getindex(p::Projection{M,K,C},i::Int) where {M,K,C<:Bond.Floating}
    ct = Bond.coupon_times(p.contract)
    l  = length(ct)
    b = p.contract
    if i > l
        throw(BoundsError(p,i))
    elseif i == l
        reference_rate = Periodic(b.frequency.frequency)(rate(p.model[b.key]))
        return Cashflow(1. + (rate(reference_rate)+ b.coupon_rate)/b.frequency.frequency, ct[i])
    else
        reference_rate = Periodic(b.frequency.frequency)(rate(p.model[b.key]))
        return Cashflow((rate(reference_rate)+ b.coupon_rate)/b.frequency.frequency, ct[i])
    end
end

function Base.iterate(p::Projection{M,K,C}) where {M,K,C<:Composite}
    itera = iterb =1
    iter = 2
    ap=@set p.contract = p.contract.a
    bp=@set p.contract = p.contract.b
    l=length(p)
    cfa = ap[itera]
    cfb = bp[iterb]
    if cfa.time == cfb.time
        res = Cashflow(cfa.amount + cfb.amount, cfa.time)
        itera += 1
        iterb += 1
    elseif cfa.time < cfb.time
        res = cfa
        itera += 1
    else
        res = cfb
        iterb += 1
    end


    state=(;iter,itera,iterb,ap,bp,l)


    return (res,state)
end
function Base.iterate(p::Projection{M,K,C},state) where {M,K,C<:Composite}
    @show state.iter,state.itera,state.iterb
    if state.iter > state.l
        return nothing
    else
        cfa = state.ap[state.itera]
        cfb = state.bp[state.iterb]
        @show cfa, cfb
        res = if cfa.time == cfb.time
            @reset state.itera += 1
            @reset state.iterb += 1
            Cashflow(cfa.amount + cfb.amount, cfa.time)
        elseif cfa.time < cfb.time
            @reset state.itera += 1
            cfa
        else
            @reset state.iterb += 1
            cfb
        end
        @reset state.iter += 1
        return (res,state)
    end
end

function Base.length(p::Projection{M,K,C}) where {M,K,C<:Composite}  
    a = Bond.coupon_times(p.contract.a)
    b = Bond.coupon_times(p.contract.b)
    length(merge_range(a,b))
end
