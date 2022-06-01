"""
An AbstractYield is an object which can be used as an argument to:

- zero-coupon spot rates via [`zero`](@ref)
- discount factor via [`discount`](@ref)
- accumulation factor via [`accumulation`](@ref)

It can be be constructed via:

- zero rate curve with [`Zero`](@ref)
- forward rate curve with [`Forward`](@ref)
- par rate curve with [`Par`](@ref)
- typical OIS curve with [`OIS`](@ref)
- typical constant maturity treasury (CMT) curve with [`CMT`](@ref)
"""
abstract type AbstractYield end

# make interest curve broadcastable so that you can broadcast over multiple`time`s in `interest_rate`
Base.Broadcast.broadcastable(ic::T) where {T<:AbstractYield} = Ref(ic)

struct YieldCurve{T,U,V} <: AbstractYield
    rates::T
    maturities::U
    zero::V # function time -> continuous zero rate
end

# internal function (will be used in EconomicScenarioGenerators)
# defines the rate output given just the type of curve
__ratetype(::Type{YieldCurve{T,U,V}}) where {T,U,V}= Yields.Rate{Float64, typeof(DEFAULT_COMPOUNDING)}