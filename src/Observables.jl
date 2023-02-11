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

struct Bond{F<:FinanceCore.CompoundingFrequency} <: AbstractBond
    coupon_rate::Float64 # coupon_rate / frequency is the actual payment amount
    frequency::F
    maturity::Float64
end

# force conversion to float
Bond(c,f,m::I) where {I<:Integer} = Bond(c,f,Float64(m))

get_frequency(a::FinanceCore.Rate; default) = a.compounding
get_frequency(a; default) = default

# assume the frequency is two or infer it from the yield
function ParYield(yield,maturity;frequency=nothing)
    if isnothing(frequency)
        frequency = get_frequency(yield;default=Periodic(2))
    end
    price = 1. # by definition for a par bond
    coupon_rate = rate(frequency(yield))
    return Quote(price,Bond(coupon_rate,frequency,maturity)) 

end

# the fixed leg of the swap
function ParSwapYield(yield,maturity;frequency=nothing)
    ParYield(yield,maturity;frequency=frequency)
end

function CMTYield(yield,maturity)
    frequency = maturity <= 1 ? Periodic(1) : Periodic(2)
    r = frequency(yield)
    return Quote(discount(r,maturity),Bond(rate(r),r.compounding,maturity))
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