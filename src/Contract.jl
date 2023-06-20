### TODO
# - if explicit price, should have par to calculate unit price
# - can we just dispatch on the type of price/yield to simplify constructor?
# - collect cashflows into vector

# Allow Dates or real timesteps
const Timepoint{T} = Union{T,Dates.Date} where {T<:Real}

abstract type AbstractContract end

struct Quote{N<:Real,T}
    price::N
    instrument::T
end

maturity(q::Quote) = maturity(q.instrument)
Base.isapprox(a::Quote, b::Quote) = isapprox(a.price, b.price) && isapprox(a.instrument, b.instrument)

struct Cashflow{N<:Real,T<:Timepoint} <: AbstractContract
    amount::N
    time::T
end

maturity(c::Cashflow) = c.time
Base.isapprox(a::Cashflow, b::Cashflow) = isapprox(a.amount, b.amount) && isapprox(a.time, b.time)
Base.convert(::Type{Cashflow{A,B}}, y::Cashflow{C,D}) where {A,B,C,D} = Cashflow(A(y.amount), B(y.time))

struct Composite{A,B} <: AbstractContract
    a::A
    b::B
end

maturity(c::Composite) = max(maturity(c.a), maturity(c.b))

### Bonds 
module Bond
import ..AbstractContract
import ..Timepoint
import ..Cashflow, ..Quote
import ..FinanceModels: maturity
using ..FinanceCore

using FinanceCore: Periodic, Continuous, Rate

export ZCBYield, ZCBPrice, ParSwapYield, ParYield, CMTYield

abstract type AbstractBond <: AbstractContract end
maturity(b::AbstractBond) = b.maturity


"""
ZCBPrice(discount,maturity)
ZCBPrice(yield::Vector)

Takes discount factors. 

Use broadcasting to create a set of quotes given a collection of prices and maturities, e.g. `ZCBPrice.(FinanceModels,maturities)`.
"""
ZCBPrice(price, time) = Quote(price, Cashflow(1.0, time))


"""
ZCBYield(yield,maturity)
ZCBYield(yield::Vector)

Takes zero (sometimes called "spot") rates. Assumes annual effective compounding (`Periodic(1)``) unless given a `Rate` with a different compounding frequency.

Use broadcasting to create a set of quotes given a collection of FinanceModels and maturities, e.g. `ZCBYield.(FinanceModels,maturities)`.
"""
ZCBYield(yield, time) = Quote(discount(yield, time), Cashflow(1.0, time))

struct Fixed{F<:FinanceCore.CompoundingFrequency,N<:Real,M<:Timepoint} <: AbstractBond
    coupon_rate::N # coupon_rate / frequency is the actual payment amount
    frequency::F
    maturity::M
end

function Base.isapprox(a::Fixed, b::Fixed)
    isapprox(a.coupon_rate, b.coupon_rate) && ==(a.frequency, b.frequency) && isapprox(a.maturity, b.maturity)
end

# function timesteps(b::AbstractBond)
#     f = 1 / b.frequency.frequency
#     f:f:b.maturity
# end


struct Floating{F<:FinanceCore.CompoundingFrequency,N<:Real,M<:Timepoint,K} <: AbstractBond
    coupon_rate::N # coupon_rate / frequency is the actual payment amount
    frequency::F
    maturity::M
    key::K
end

__coerce_periodic(y::Periodic) = y
__coerce_periodic(y::T) where {T<:Int} = Periodic(y)

"""
ParYield(yield,maturity)
ParYield(yield::Vector)

Takes bond equivalent FinanceModels, and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual. Alternative, you may pass a `Rate` as the yield and the coupon frequency will be inferred from the `Rate`'s frequency. 

Use broadcasting to create a set of quotes given a collection of FinanceModels and maturities, e.g. `ParYield.(FinanceModels,maturities)`.
"""
function ParYield(yield, maturity; frequency=Periodic(2))
    # assume the frequency is two or infer it from the yield
    frequency = __coerce_periodic(frequency)
    price = 1.0 # by definition for a par bond
    coupon_rate = rate(frequency(yield))
    return Quote(price, Fixed(coupon_rate, frequency, maturity))
end
function ParYield(yield::Rate{N,T}, maturity; frequency=Periodic(2)) where {T<:Periodic,N}
    frequency = yield.compounding
    price = 1.0 # by definition for a par bond
    coupon_rate = rate(frequency(yield))
    return Quote(price, Fixed(coupon_rate, frequency, maturity))
end

# the fixed leg of the swap
function ParSwapYield(yield, maturity; frequency=Periodic(4))
    frequency = __coerce_periodic(frequency)
    ParYield(yield, maturity; frequency=frequency)
end

"""
CMTYield(yield,maturity)
CMTYield(yield::Vector)
Takes constant maturity (treasury) FinanceModels (bond equivalent), and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual.

Use broadcasting to create a set of quotes given a collection of FinanceModels and maturities, e.g. `CMTYield.(FinanceModels,maturities)`.
"""
function CMTYield(yield, maturity)
    # Assume maturity < 1 don't pay coupons and are therefore discount bonds
    # Assume maturity > 1 pay coupons and are therefore par bonds
    frequency = Periodic(2)
    r, v = if maturity ≤ 1
        Periodic(0.0, 1), discount(yield, maturity)
    else
        # coupon paying par bond 
        frequency(yield), 1.0
    end
    return Quote(v, Fixed(rate(r), r.compounding, maturity))
end

"""
OISYield(yield [, maturity=eachindex(yield)]))

Assumes that maturities less than or equal to 12 months are settled once (per Hull textbook, 4.7), otherwise quarterly and that the FinanceModels given are bond equivalent.

Use broadcasting to create a set of quotes given a collection of FinanceModels and maturities, e.g. `OISYield.(FinanceModels,maturities)`.

"""
function OISYield(yield, maturity=eachindex(yield))

    if maturity <= 1
        return Quote(discount(yield, maturity), Fixed(0.0, Periodic(1), maturity))
    else
        frequency = Periodic(4)
        r = frequency(yield)
        return Quote(1.0, Fixed(rate(r), frequency, maturity))
    end
end

"""
ForwardYields(yields,times) 
Returns a vector of `Quote` corresponding to the . 
    
# Examples
```julia-repl
julia> FinanceModels.Bond.ForwardYields([0.01,0.02],[1.,3.])
2-element Vector{Quote{Float64, Cashflow{Float64, Float64}}}:
 Quote{Float64, Cashflow{Float64, Float64}}(0.9900990099009901, Cashflow{Float64, Float64}(1.0, 1.0))
 Quote{Float64, Cashflow{Float64, Float64}}(0.9423223345470445, Cashflow{Float64, Float64}(1.0, 3.0))
```
"""
function ForwardYields(yields, times=eachindex(yields))
    df = 1.0
    t_prior = 0.0
    map(zip(yields, times)) do (y, t)
        df *= discount(y, t - t_prior)
        t_prior = t
        Quote(
            df,
            Cashflow(1.0, t)
        )
    end
end


# Bond utility funcs

function coupon_times(maturity, frequency)
    Δt = min(1 / frequency, maturity)
    times = maturity:-Δt:0
    f = last(times)
    f += iszero(f) ? Δt : zero(f)
    l = first(times)
    return f:Δt:l
end
coupon_times(b::AbstractBond) = coupon_times(b.maturity, b.frequency.frequency)


for op = (:ZCBPrice, :ZCBYield, :ParYield, :ParSwapYield, :CMTYield, :ForwardYield)
    eval(quote
        $op(x::Vector; kwargs...) = $op.(x, eachindex(x); kwargs...)
    end)
end


end

struct CommonEquity <: AbstractContract end

module Option
import ..AbstractContract
import ..Timepoint

struct EuroCall{S<:AbstractContract,K<:Real,M<:Timepoint} <: AbstractContract
    underlying::S
    strike::K
    maturity::M
end
end

"""
Forward(time,instrument)

The instrument is relative to the Forward time.
e.g. if you have a `Forward(1.0, Cashflow(1.0, 3.0))` then the instrument is a cashflow that pays 1.0 at time 4.0
"""
struct Forward{T<:Timepoint,I<:AbstractContract} <: AbstractContract
    time::T
    instrument::I
end



# convert ZCB to non-forward versions
# function __process_forwards(qs::Vector{Quote{U,Forward{N,T}}}) where {N,T<:Cashflow,U}
#     v = 1.0
#     t = 0.0
#     map(qs) do q
#         v *= (q.price / q.instrument.instrument.amount)
#         t = q.instrument.time + q.instrument.instrument.time
#         Quote(v, Cashflow(1.0, t))
#     end
# end



# # cashflows should be a vector of a vector of cashflows
# function cashflow_matrix(instruments::Vector{Q}; resolution=1000) where {Q<:AbstractContract}
#     vcf = collect.(instruments) # a vector of vector of cashflows
#     ts = timepoints(instruments; resolution=resolution)
#     m = zeros(round(Int, last(ts) ÷ step(ts)), length(vcf))
#     for (i, cf) in enumerate(vcf)
#         for c in cf
#             m[round(Int, c.time ÷ step(ts)), i] = c.amount
#         end
#     end
#     return m
#     # for each obs determine the closest integer multiple of the gcd
#     # fill in the matrix
# end

# function cashflow_matrix(quotes::Vector{Quote{T,U}}; resolution=1000) where {T,U}
#     cashflow_matrix([q.instrument for q in quotes]; resolution=resolution)
# end

# function timepoints(instruments::Vector{Q}; resolution=1000) where {Q<:AbstractContract}
#     # calculate the gcd of the timepoints 
#     mapreduce(timepoints, merge_range, instruments)
# end

# function timepoints(quotes::Vector{Q}; resolution=1000) where {Q<:Quote}
#     timepoints([q.instrument for q in quotes]; resolution=resolution)
# end

# timepoints(c::Cashflow) = c.time:c.time
# timepoints(c::Bond.AbstractBond) = Bond.coupon_times(c)



# function merge_range(a, b; resolution=1000)
#     start = min(minimum(a), minimum(b))
#     last = max(maximum(a), maximum(b))
#     delta = gcd(Int(step(a) * resolution), Int(step(b) * resolution)) / resolution
#     return start:delta:last
# end

# create a matrix of cashflows and a vector of timepoints
# timepoints need not be spaced evenly
function cashflows_timepoints(qs)
    cfs = map(q -> collect(q.instrument), qs)
    times = map(cfs) do cf
                map(c -> c.time, cf)
            end |> Iterators.flatten |> unique |> sort!

    m = zeros(length(qs), length(times))

    for t in 1:length(times)
        for q in 1:length(qs)
            if times[t] == qs[q].instrument.time
                m[q, t] += qs[q].instrument.amount
            end
        end
    end
    m
    return m, times
end