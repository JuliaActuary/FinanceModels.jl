"""
An AbstractYield is an object which can be used as an argument to:

- zero-coupon spot rates via [`zero`](@ref)
- forward zero rates via [`forward`](@ref)
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

abstract type AbstractYieldCurve <: AbstractYield end

abstract type YieldCurveFitParameters end
Base.Broadcast.broadcastable(x::T) where {T<:YieldCurveFitParameters} = Ref(x)