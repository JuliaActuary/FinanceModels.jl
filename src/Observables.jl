### TODO
# - if explicit price, should have par to calculate unit price
# - can we just dispatch on the type of price/yield to simplify constructor?
# - collect cashflows into vector

abstract type Instrument end

struct Quote{N<:Real,T<:Instrument} 
    price::N
    instrument::T
end
Base.isapprox(a::Quote,b::Quote) = isapprox(a.price,b.price) && isapprox(a.instrument,b.instrument)

struct Cashflow{N<:Real,T<:Real} <: Instrument
    amount::N
    time::T
end
Base.isapprox(a::Cashflow,b::Cashflow) = isapprox(a.amount,b.amount) && isapprox(a.time,b.time)


abstract type AbstractBond <: Instrument end

struct Composite{T,U} <: Instrument 
    first::T
    second::U
end

Base.:+(a::Instrument,b::Instrument) = Composite(a,b)

function Composite(a::Cashflow,b::Cashflow)
    if a.time == b.time
        return Cashflow(a.amount + b.amount, a.time)
    else
        return Composite{Cashflow,Cashflow}(a,b)
    end
end

### Bonds 
ZCBPrice(price,time) = Quote(price,Cashflow(1.,time))
ZCBPrice(prices) =  ZCBPrice.(prices,eachindex(prices))
ZCBYield(yield,time) = Quote(discount(yield,time),Cashflow(1.,time)) 
ZCBYield(yields) = ZCBYield.(yields,eachindex(yields))

struct Bond{F<:FinanceCore.CompoundingFrequency,N<:Real,M<:Real} <: AbstractBond
    coupon_rate::N # coupon_rate / frequency is the actual payment amount
    frequency::F
    maturity::M
end

function Base.isapprox(a::Bond,b::Bond)
    isapprox(a.coupon_rate,b.coupon_rate) && ==(a.frequency,b.frequency) && isapprox(a.maturity,b.maturity)
end

coupon_times(b::Bond) = coupon_times(b.maturity,b.frequency.frequency)

function Base.iterate(b::Bond,state=(coupon_times(b),1)) 
    if state[2] > lastindex(state[1])
        return nothing
    elseif state[2] == lastindex(state[1])
        return (Cashflow(1. + b.coupon_rate/b.frequency.frequency,state[1][end]), (state[1],state[2]+1))
    else
        return (Cashflow(b.coupon_rate/b.frequency.frequency,state[1][state[2]]), (state[1],state[2] + 1))
    end
end

Base.length(b::Bond) = length(coupon_times(b))
Base.eltype(::Type{Bond{F,N,M}}) where {F,N,M} = Cashflow{N,M}

function _pv(y,b::Bond)
    return sum(cf.amount * discount(y,cf.time) for cf in b)
end

__coerce_periodic(y::Periodic) = y
__coerce_periodic(y::T) where {T<:Int} = Periodic(y)

# assume the frequency is two or infer it from the yield
function ParYield(yield,maturity;frequency=Periodic(2))
    frequency = __coerce_periodic(frequency)
    price = 1. # by definition for a par bond
    coupon_rate = rate(frequency(yield))
    return Quote(price,Bond(coupon_rate,frequency,maturity)) 
end

function ParYield(yields;frequency=Periodic(2))
    frequency = __coerce_periodic(frequency)
    return ParYield.(yields,eachindex(yields);frequency=frequency)
end

# the fixed leg of the swap
function ParSwapYield(yield,maturity;frequency=Periodic(4))
    frequency = __coerce_periodic(frequency)
    ParYield(yield,maturity;frequency=frequency)
end

# assume maturity < 1 don't pay coupons and are therefore discount bonds
# assume maturity > 1 pay coupons and are therefore par bonds
# floating point yield assumed to be annual effective (`Periodic(1)`) unless otherwise specified
function CMTYield(yield,maturity)
    frequency =  Periodic(2)
    r, v = if maturity ≤ 1
        Periodic(0.,1), discount(yield,maturity) 
    else
        # coupon paying par bond 
        frequency(yield), 1.0
    end
    return Quote(v,Bond(rate(r),r.compounding,maturity))
end

function OISYield(yield,maturity)
    frequency = if m <= 1
        Periodic(1 / m)
    else
        Periodic(4)
    end
    r = frequency(yield)
    return Quote(discount(r,maturity),Bond(rate(r),r.compounding,maturity))
end


# the instrument is relative to the Forward time.
# e.g. if you have a Forward(1.0, Cashflow(1.0, 3.0)) then the instrument is a cashflow that pays 1.0 at time 4.0
struct Forward{N<:Real,I<:Instrument} <: Instrument
    time::N
    instrument::I
end


ForwardYield(yield,start=0.0,duration=1.0) = Quote(discount(yield,duration),Forward(start,Cashflow(1.,duration)))
