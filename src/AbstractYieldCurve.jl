abstract type AbstractYieldCurve <: FinanceCore.AbstractYield end
# make interest curve broadcastable so that you can broadcast over multiple`time`s in `interest_rate`
Base.Broadcast.broadcastable(ic::T) where {T<:AbstractYieldCurve} = Ref(ic)

"""
A `CurveMethod` is a structure which contains associated parameters for a yield curve fitting procedure. The type of the object determines the method, and the values of the object determine the parameters.

If the fitting data and the rates are passed as `<:Real` numbers instead of a type of `Rate`s, the default interpretation may vary depending on the fitting type/parameter. See the individual docstrings of the types for more information.

Available types are:

- [`Bootstrap`](@ref)
- [`NelsonSiegel`](@ref)
- [`NelsonSiegelSvensson`](@ref)
"""
abstract type CurveMethod end

Base.Broadcast.broadcastable(x::T) where {T<:CurveMethod} = Ref(x)