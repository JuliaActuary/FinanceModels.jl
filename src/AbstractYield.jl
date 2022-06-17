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

"""
A `YieldCurveFitParameters` is a structure which contains associated parameters for a yield curve fitting procedure. The type of the object determines the method, and the values of the object determine the parameters.

If the fitting data and the rates are passed as `<:Real` numbers instead of a type of `Rate`s, the default interpretation may vary depending on the fitting type/parameter. See the individual docstrings of the types for more information.

Available types are:

- [`Bootstrap`](@ref)
- [`NelsonSiegel`](@ref)
- [`NelsonSiegelSvensson`](@ref)
"""
abstract type YieldCurveFitParameters end

Base.Broadcast.broadcastable(x::T) where {T<:YieldCurveFitParameters} = Ref(x)