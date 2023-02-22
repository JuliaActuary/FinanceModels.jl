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
function Base.iterate(cf::Cashflow,state=cf) 
    if isnothing(state)
        return nothing
    else
        return (state,nothing)
    end
end
Base.length(cf::Cashflow) = 1
Base.eltype(x::Type{Cashflow{N,T}}) where{N,T} = x

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
"""
    ZCBYield(yield,maturity)
    ZCBYield(yield::Vector)

Takes zero (sometimes called "spot") rates. Assumes annual effective compounding (`Periodic(1)``) unless given a `Rate` with a different compounding frequency.

Use broadcasting to create a set of quotes given a collection of yields and maturities, e.g. `ZCBYield.(yields,maturities)`.
"""
ZCBPrice(price,time) = Quote(price,Cashflow(1.,time))

"""
    ZCBPrice(discount,maturity)
    ZCBPrice(yield::Vector)

Takes discount factors. 

Use broadcasting to create a set of quotes given a collection of prices and maturities, e.g. `ZCBPrice.(yields,maturities)`.
"""
ZCBYield(yield,time) = Quote(discount(yield,time),Cashflow(1.,time)) 

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

"""
    ParYield(yield,maturity)
    ParYield(yield::Vector)

Takes bond equivalent yields, and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual. Alternative, you may pass a `Rate` as the yield and the coupon frequency will be inferred from the `Rate`'s frequency. 

Use broadcasting to create a set of quotes given a collection of yields and maturities, e.g. `ParYield.(yields,maturities)`.
"""
function ParYield(yield,maturity;frequency=Periodic(2))
    # assume the frequency is two or infer it from the yield
    frequency = __coerce_periodic(frequency)
    price = 1. # by definition for a par bond
    coupon_rate = rate(frequency(yield))
    return Quote(price,Bond(coupon_rate,frequency,maturity)) 
end
function ParYield(yield::Rate{N,T},maturity;frequency=Periodic(2)) where {T<:Periodic,N}
    frequency = yield.compounding
    price = 1. # by definition for a par bond
    coupon_rate = rate(frequency(yield))
    return Quote(price,Bond(coupon_rate,frequency,maturity)) 
end

# the fixed leg of the swap
function ParSwapYield(yield,maturity;frequency=Periodic(4))
    frequency = __coerce_periodic(frequency)
    ParYield(yield,maturity;frequency=frequency)
end

"""
    CMTYield(yield,maturity)
    CMTYield(yield::Vector)
Takes constant maturity (treasury) yields (bond equivalent), and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual.

Use broadcasting to create a set of quotes given a collection of yields and maturities, e.g. `CMTYield.(yields,maturities)`.
"""
function CMTYield(yield,maturity)
    # Assume maturity < 1 don't pay coupons and are therefore discount bonds
    # Assume maturity > 1 pay coupons and are therefore par bonds
    frequency =  Periodic(2)
    r, v = if maturity ≤ 1
        Periodic(0.,1), discount(yield,maturity) 
    else
        # coupon paying par bond 
        frequency(yield), 1.0
    end
    return Quote(v,Bond(rate(r),r.compounding,maturity))
end

"""
    OISYield(yield,maturity)
    OISYield(yield::Vector)

Assumes that maturities less than or equal to 12 months are settled once (per Hull textbook, 4.7), otherwise quarterly and that the yields given are bond equivalent.

Use broadcasting to create a set of quotes given a collection of yields and maturities, e.g. `OISYield.(yields,maturities)`.

"""
function OISYield(yield,maturity)
        
    if maturity <= 1
        return Quote(discount(yield,maturity),Bond(0.,Periodic(1),maturity))
    else
        frequency = Periodic(4)
        r = frequency(yield)
        return Quote(1.0,Bond(rate(r),frequency,maturity))
    end
end


"""
    Forward(time,instrument)

The instrument is relative to the Forward time.
e.g. if you have a `Forward(1.0, Cashflow(1.0, 3.0))` then the instrument is a cashflow that pays 1.0 at time 4.0
"""
struct Forward{N<:Real,I<:Instrument} <: Instrument
    time::N
    instrument::I
end


"""
    ForwardYield(yield,to=1.0,from=to-1.0) 
Returns a `Quote`d price for a future cashflow. 

# Examples
```julia
fy = ForwardYield.([0.01,0.02],[1.,2.])
first(fy) == Quote(1/1.01,Forward(0.0,Cashflow(1.,1.)))
last(fy) == Quote(1/1.02,Forward(1.0,Cashflow(1.,1.)))
```
"""
ForwardYield(yield,to=1.0,from=to-1.0) = Quote(discount(yield,to-from),Forward(from,Cashflow(1.,to-from)))

# convert ZCB to non-forward versions
function __process_forwards(qs::Vector{Quote{U,Forward{N,T}}}) where {N,T<:Cashflow,U}
    v = 1.0
    t = 0.0
    map(qs) do q
        v *= (q.price / q.instrument.instrument.amount)
        t = q.instrument.time + q.instrument.instrument.time
        Quote(v,Cashflow(1.,t))
    end
end


for op = (:ZCBPrice, :ZCBYield, :ParYield, :ParSwapYield, :CMTYield, :ForwardYield)
    eval(quote
        $op(x::Vector{T};y...) where {T} = $op.(x,eachindex(x);y...)
    end)
end


# cashflows should be a vector of a vector of cashflows
function cashflow_matrix(instruments::Vector{Q};resolution=1000) where {Q<:Instrument} 
    vcf = collect.(instruments) # a vector of vector of cashflows
    ts = timesteps(instruments;resolution=resolution)
    m = zeros(round(Int,last(ts) ÷ step(ts)), length(vcf))
    for (i, cf) in enumerate(vcf)
        for c in cf
            m[ round(Int,c.time ÷ step(ts)),i] = c.amount
        end
    end
    return m
    # for each obs determine the closest integer multiple of the gcd
    # fill in the matrix
end


function timesteps(instruments::Vector{Q};resolution=1000) where {Q<:Instrument} 
    # calculate the gcd of the timesteps 
    vcf = collect.(instruments) # a vector of vector of cashflows
    scaled_step = gcd([Int(cf.time * resolution) for cf in Iterators.flatten(vcf) if cf.time > 0])
    last = Int(maximum(cf.time * resolution for cf in Iterators.flatten(vcf)) ÷ resolution)
    (scaled_step / resolution): (scaled_step / resolution) : (last)
end