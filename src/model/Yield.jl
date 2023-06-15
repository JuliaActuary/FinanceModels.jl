
module Yield
import ..AbstractModel
import ..FinanceCore
using FinanceCore: Continuous, discount

export discount, zero, forward, par, pv

abstract type AbstractYieldModel <: AbstractModel end


struct Constant{R} <: AbstractYieldModel
    rate::R
end

function Constant(rate::R) where {R<:Real}
    Constant(FinanceCore.Rate(rate))
end

Constant() = Constant(0.0)

FinanceCore.discount(c::Constant, t) = FinanceCore.discount(c.rate, t)

function pv(model::M, c::AbstractContract; cur_time=0.0) where {M<:Yield.AbstractYieldModel}
    p = Projection(c, model, CashflowProjection())
    xf = p |> Filter(cf -> cf.time >= cur_time) |> Map(cf -> FinanceCore.discount(model, cf.time - cur_time) * cf.amount)
    foldxl(+, xf)
end

## Generic and Fallbacks
"""
    discount(yc, to)
    discount(yc, from,to)

The discount factor for the yield curve `yc` for times `from` through `to`.
"""
FinanceCore.discount(yc::T, from, to) where {T<:AbstractYieldModel} = discount(yc, to) / discount(yc, from)

"""
    forward(yc, from, to)

The forward `Rate` implied by the yield curve `yc` between times `from` and `to`.
"""
function FinanceCore.forward(yc::T, from, to=from + 1) where {T<:AbstractYieldModel}
    return forward(yc, from, to,)
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
    coup_times = coupon_times(time, frequency)
    coupon_pv = sum(discount(curve, t) for t in coup_times)
    Δt = step(coup_times)
    r = (1 - mat_disc) / coupon_pv
    cfs = [t == last(coup_times) ? 1 + r : r for t in coup_times]
    # `sign(r)`` is used instead of `1` because there are times when the coupons are negative so we want to flip the sign
    cfs = [-1; cfs]
    r = internal_rate_of_return(cfs, [0; coup_times])
    frequency_inner = min(1 / Δt, max(1 / Δt, frequency))
    r = convert(Periodic(frequency_inner), r)
    return r
end

"""
    zero(curve,time)
    zero(curve,time,CompoundingFrequency)

Return the zero rate for the curve at the given time.
"""
function Base.zero(c::YC, time) where {YC<:AbstractYieldModel}
    zero(c, time, FinanceCore.CompoundingFrequency(c))
end

function Base.zero(c::YC, time, cf::C) where {YC<:AbstractYieldModel,C<:FinanceCore.CompoundingFrequency}
    df = discount(c, time)
    r = -log(df) / time
    return Continuous(r)
end

"""
    accumulation(yc, from, to)

The accumulation factor for the yield curve `yc` for times `from` through `to`.
"""
function FinanceCore.accumulation(yc::AbstractYieldModel, time)
    return 1 ./ discount(yc, time)
end

function FinanceCore.accumulation(yc::AbstractYieldModel, from, to)
    return 1 ./ discount(yc, from, to)
end

end
