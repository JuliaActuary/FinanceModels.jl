## Generic and Fallbacks
"""
    discount(yc, to)
    discount(yc, from,to)

The discount factor for the yield curve `yc` for times `from` through `to`.
"""
FinanceCore.discount(yc::T, from, to) where {T<:AbstractYieldCurve}= discount(yc, to) / discount(yc, from)

"""
    forward(yc, from, to, CompoundingFrequency=Periodic(1))

The forward `Rate` implied by the yield curve `yc` between times `from` and `to`.
"""
function FinanceCore.forward(yc::T, from, to) where {T<:AbstractYieldCurve}
    return forward(yc, from, to, DEFAULT_COMPOUNDING)
end

function FinanceCore.forward(yc::T, from, to, cf::FinanceCore.CompoundingFrequency) where {T<:AbstractYieldCurve}
    r = Periodic((accumulation(yc, to) / accumulation(yc, from))^(1 / (to - from)) - 1, 1)
    return convert(cf, r)
end

function FinanceCore.forward(yc::T, from) where {T<:AbstractYieldCurve}
    to = from + 1
    return forward(yc, from, to)
end

function FinanceCore.CompoundingFrequency(curve::T) where {T<:AbstractYieldCurve}
    return DEFAULT_COMPOUNDING
end


"""
    par(curve,time;frequency=2)

Calculate the par yield for maturity `time` for the given `curve` and `frequency`. Returns a `Rate` object with periodicity corresponding to the `frequency`. The exception to this is if `time` is less than what the payments allowed by frequency (e.g. a time `0.5` but with frequency `1`) will effectively assume frequency equal to 1 over `time`.

# Examples

```julia-repl
julia> c = Yields.Constant(0.04);

julia> Yields.par(c,4)
Yields.Rate{Float64, Yields.Periodic}(0.03960780543711406, Yields.Periodic(2))

julia> Yields.par(c,4;frequency=1)
Yields.Rate{Float64, Yields.Periodic}(0.040000000000000036, Yields.Periodic(1))

julia> Yields.par(c,0.6;frequency=4)
Yields.Rate{Float64, Yields.Periodic}(0.039413626195875295, Yields.Periodic(4))

julia> Yields.par(c,0.2;frequency=4)
Yields.Rate{Float64, Yields.Periodic}(0.039374942589460726, Yields.Periodic(5))

julia> Yields.par(c,2.5)
Yields.Rate{Float64, Yields.Periodic}(0.03960780543711406, Yields.Periodic(2))
```
"""
function par(curve, time; frequency=2)
    mat_disc = discount(curve, time)
    coup_times = coupon_times(time,frequency)
    coupon_pv = sum(discount(curve,t) for t in coup_times)
    Δt = step(coup_times)
    r = (1-mat_disc) / coupon_pv
    cfs = [t == last(coup_times) ? 1+r : r for t in coup_times]
    # `sign(r)`` is used instead of `1` because there are times when the coupons are negative so we want to flip the sign
    cfs = [-1;cfs]
    r = internal_rate_of_return(cfs,[0;coup_times])
    frequency_inner = min(1/Δt,max(1 / Δt, frequency))
    r = convert(Periodic(frequency_inner),r)
    return r
end

"""
    zero(curve,time)

Return the zero rate for the curve at the given time.
"""
function Base.zero(c::YC, time) where {YC<:AbstractYieldCurve}
    df = discount(c, time)
    r = -log(df)/time
    return Continuous(r) # c.zero is a curve of continuous rates represented as floats. explicitly wrap in continuous before converting
end

"""
    accumulation(yc, from, to)

The accumulation factor for the yield curve `yc` for times `from` through `to`.
"""
function FinanceCore.accumulation(yc::AbstractYieldCurve, time)
    return 1 ./ discount(yc, time)
end

function FinanceCore.accumulation(yc::AbstractYieldCurve, from, to)
    return 1 ./ discount(yc, from, to)
end


