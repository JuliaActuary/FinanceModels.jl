
## Generic and Fallbacks
"""
    discount(yc, to)
    discount(yc, from,to)

The discount factor for the yield curve `yc` for times `from` through `to`.
"""
discount(yc::T, from, to) where {T<:AbstractYieldCurve}= discount(yc, to) / discount(yc, from)

"""
    forward(yc, from, to, CompoundingFrequency=Periodic(1))

The forward `Rate` implied by the yield curve `yc` between times `from` and `to`.
"""
function forward(yc::T, from, to) where {T<:AbstractYieldCurve}
    return forward(yc, from, to, DEFAULT_COMPOUNDING)
end

function forward(yc::T, from, to, cf::CompoundingFrequency) where {T<:AbstractYieldCurve}
    r = Periodic((accumulation(yc, to) / accumulation(yc, from))^(1 / (to - from)) - 1, 1)
    return convert(cf, r)
end

function forward(yc::T, from) where {T<:AbstractYieldCurve}
    to = from + 1
    return forward(yc, from, to)
end

function CompoundingFrequency(curve::T) where {T<:AbstractYieldCurve}
    return DEFAULT_COMPOUNDING
end

"""
    zero(curve,time)
    zero(curve,time,CompoundingFrequency)

Return the zero rate for the curve at the given time.
"""
function Base.zero(c::YC, time) where {YC<:AbstractYieldCurve} 
     zero(c, time, CompoundingFrequency(c))
end

function Base.zero(c::YC, time, cf::C) where {YC<:AbstractYieldCurve,C<:CompoundingFrequency}
    df = discount(c, time)
    r = -log(df)/time
    return convert(cf, Continuous(r)) # c.zero is a curve of continuous rates represented as floats. explicitly wrap in continuous before converting
end

"""
    accumulation(yc, from, to)

The accumulation factor for the yield curve `yc` for times `from` through `to`.
"""
function accumulation(yc::AbstractYieldCurve, time)
    return 1 ./ discount(yc, time)
end

function accumulation(yc::AbstractYieldCurve, from, to)
    return 1 ./ discount(yc, from, to)
end
