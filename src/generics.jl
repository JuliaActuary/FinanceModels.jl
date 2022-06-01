
## Generic and Fallbacks
"""
    discount(yc, to)
    discount(yc, from,to)

The discount factor for the yield curve `yc` for times `from` through `to`.
"""
discount(yc::YieldCurve, time) = exp(-yc.zero(time) * time)
discount(yc::AbstractYield, from, to) = discount(yc, to) / discount(yc, from)

"""
    forward(yc, from, to, CompoundingFrequency=Periodic(1))

The forward `Rate` implied by the yield curve `yc` between times `from` and `to`.
"""
function forward(yc::T, from, to) where {T<:AbstractYield}
    return forward(yc, from, to, DEFAULT_COMPOUNDING)
end

function forward(yc::T, from, to, cf::CompoundingFrequency) where {T<:AbstractYield}
    r = Periodic((accumulation(yc, to) / accumulation(yc, from))^(1 / (to - from)) - 1, 1)
    return convert(cf, r)
end

function forward(yc::T, from) where {T<:AbstractYield}
    to = from + 1
    return forward(yc, from, to)
end

"""
    accumulation(yc, from, to)

The accumulation factor for the yield curve `yc` for times `from` through `to`.
"""
function accumulation(yc::AbstractYield, time)
    return 1 ./ discount(yc, time)
end

function accumulation(yc::AbstractYield, from, to)
    return 1 ./ discount(yc, from, to)
end
