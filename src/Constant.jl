
"""
    Constant(rate::Real, cf::CompoundingFrequency=Periodic(1))
    Constant(r::Rate)

Construct a yield object where the spot rate is constant for all maturities. If `rate` is not a `Rate` type, will assume `Periodic(1)` for the compounding frequency

# Examples

```julia-repl
julia> y = Yields.Constant(0.05)
julia> FinanceCore.discount(y,2)
0.9070294784580498     # 1 / (1.05) ^ 2
```
"""
struct Constant{T} <: AbstractYieldCurve
    rate::T
    Constant(rate::T) where {T<:Rate} = new{T}(rate)
end

__ratetype(::Type{Constant{T}}) where {T} = T
__default_rate_interpretation(::Type{Constant},r) = Periodic(r,1)
FinanceCore.CompoundingFrequency(c::Constant{T}) where {T} = c.rate.compounding

function Constant(rate::T) where {T<:Real}
    r = __default_rate_interpretation(Constant,rate)
    return Constant(r)
end

Base.zero(c::Constant, time) = c.rate
Base.zero(c::Constant, time, cf::FinanceCore.CompoundingFrequency) = convert(cf, c.rate)
FinanceCore.rate(c::Constant) = c.rate
FinanceCore.rate(c::Constant, time) = c.rate
FinanceCore.discount(r::Constant, time) = FinanceCore.discount(r.rate, time)
FinanceCore.discount(r::Constant, from, to) = FinanceCore.discount(r.rate, to - from)
FinanceCore.accumulation(r::Constant, time) = accumulation(r.rate, time)
FinanceCore.accumulation(r::Constant, from, to) = accumulation(r.rate, to - from)
