### TODO
# - if explicit price, should have par to calculate unit price
# - can we just dispatch on the type of price/yield to simplify constructor?
# - collect cashflows into vector

abstract type Instrument end

struct Quote{T<:Instrument} 
    price::Float64
    instrument::T
end

struct Cashflow <: Instrument
    amount::Float64
    time::Float64
end


abstract type AbstractBond <: Instrument end
struct Bond
end

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
ZCBYield(yield,time) = Quote(discount(yield,time),Cashflow(1.,time)) 
