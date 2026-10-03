"""
    Constant(rate)

A yield curve representing a flat term structure. `rate` can be a [`Rate`](@ref) object or a `Real` object.


If [`fit`](@ref FinanceModels.fit)ing with the default FinanceModels.jl settings, the solver will attempt to fit a discount rate with the range of: `-1.0 .. 1.0`
"""
struct Constant{R} <: AbstractYieldModel
    rate::R
end

function Constant(rate::R) where {R <: Real}
    return Constant(FinanceCore.Rate(rate))
end

Constant() = Constant(0.0)

FinanceCore.discount(c::Constant, t) = FinanceCore.discount(c.rate, t)

# The continuous zero rate of a flat curve is its (continuous) rate at every tenor, including t = 0.
# Defining `zero` directly skips the generic L(t)/t and keeps curves composed from a `Constant` in
# zero-rate space (see `CompositeYield`/`ScaledYield`).
Base.zero(c::Constant, t) = convert(Continuous(), c.rate)
__log_discount(c::Constant, t) = __zero_log_discount(c, t)
__log_native(::Constant) = true
__instantaneous_forward(c::Constant, t) = FinanceCore.rate(Base.zero(c, t))
function __log_tail(c::Constant)
    z = FinanceCore.rate(Base.zero(c, Inf))
    return __LogTail(zero(z), z, zero(z))
end
