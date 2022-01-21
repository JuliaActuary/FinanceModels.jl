
## Generic and Fallbacks
"""
    discount(rate,to)
    discount(rate,from,to)

The discount factor for the `rate` for times `from` through `to`. If rate is a `Real` number, will assume a `Constant` interest rate.
"""
discount(yc, time) = exp(-yc.zero(time) * time)
discount(rate::Rate{<:Real,<:CompoundingFrequency}, from, to) = discount(Constant(rate), from, to)
discount(rate::Rate{<:Real,<:CompoundingFrequency}, to) = discount(Constant(rate), to)



discount(yc, from, to) = discount(yc, to) / discount(yc, from)

"""
    forward(curve,from,to,CompoundingFrequency=Periodic(1))

The forward `Rate` implied by the curve between times `from` and `to`.
"""
function forward(yc, from, to)
    return forward(yc, from, to, Periodic(1))
end

function forward(yc, from, to, cf::T) where {T<:CompoundingFrequency}

    r = Periodic((accumulation(yc, to) / accumulation(yc, from))^(1 / (to - from)) - 1, 1)
    return convert(cf, r)
end

function forward(yc, from)
    to = from + 1
    return forward(yc, from, to)
end

"""
    accumulation(rate,from,to)

The accumulation factor for the `rate` for times `from` through `to`. If rate is a `Real` number, will assume a `Constant` interest rate.
"""
function accumulation(y::T, time) where {T<:AbstractYield}
    return 1 ./ discount(y, time)
end
accumulation(rate::Rate{<:Real,<:CompoundingFrequency}, to) = accumulation(Constant(rate), to)

function accumulation(y::T, from, to) where {T<:AbstractYield}
    return 1 ./ discount(y, from, to)
end
accumulation(rate::Rate{<:Real,<:CompoundingFrequency}, from, to) = accumulation(Constant(rate), from, to)