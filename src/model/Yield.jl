
module Yield
import ..AbstractModel
import ..FinanceCore
import ..AbstractContract

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

Return the zero rate for the curve at the given time.
"""
function Base.zero(c::YC, time) where {YC<:AbstractYieldModel}
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

## Curve Manipulations
"""
    CompositeYield(curve1,curve2,operation)

Creates a datastructure that will perform the given `operation` after independently calculating the effects of the two curves. 
Can only be created via the public API by using the `+`, `-`, `*`, and `/` operatations on `AbstractYield` objects.

As this is double the normal operations when performing calculations, if you are using the curve in performance critical locations, you should consider transforming the inputs and 
constructing a single curve object ahead of time.
"""
struct CompositeYield{T,U,V} <: AbstractYieldModel
    r1::T
    r2::U
    op::V
end

FinanceCore.rate(rc::CompositeYield, time) = rc.op(rate(rc.r1, time), rate(rc.r2, time))

function FinanceCore.discount(rc::CompositeYield, time)
    a1 = discount(rc.r1, time)^(-1 / time) - 1
    a2 = discount(rc.r2, time)^(-1 / time) - 1
    return 1 / (1 + rc.op(a1, a2))^time
end

Base.zero(rc::CompositeYield, time) = zero(rc, time)
function Base.zero(rc::CompositeYield, time, cf::C) where {C<:FinanceCore.CompoundingFrequency}
    d = discount(rc, time)
    return Periodic(1 / d^(1 / time) - 1, 1)
end

"""
    Yields.AbstractYieldModel + Yields.AbstractYieldModel

The addition of two yields will create a `CompositeYield`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will be added together.
"""
function Base.:+(a::AbstractYieldModel, b::AbstractYieldModel)
    return CompositeYield(a, b, +)
end

function Base.:+(a::Constant, b::Constant)
    return Constant(a.rate + b.rate)
end

function Base.:+(a::T, b) where {T<:AbstractYieldModel}
    return a + Constant(b)
end

function Base.:+(a, b::T) where {T<:AbstractYieldModel}
    return Constant(a) + b
end

"""
    Yields.AbstractYieldModel * Yields.AbstractYieldModel

The multiplication of two yields will create a `CompositeYield`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will be added together. This can be useful, for example, if you wanted to after-tax a yield.

# Examples

```julia-repl
julia> m = Yields.Constant(0.01) * 0.79;

julia> accumulation(m,1)
1.0079

julia> accumulation(.01*.79,1)
1.0079
```
"""
function Base.:*(a::AbstractYieldModel, b::AbstractYieldModel)
    return CompositeYield(a, b, *)
end

function Base.:*(a::Constant, b::Constant)
    a_kind = rate(a).compounding
    rate_new_basis = rate(convert(a_kind, rate(b)))
    return Constant(
        Rate(
            rate(a.rate) * rate_new_basis,
            a_kind
        )
    )
end

function Base.:*(a::T, b) where {T<:AbstractYieldModel}
    return a * Constant(b)
end

function Base.:*(a, b::T) where {T<:AbstractYieldModel}
    return Constant(a) * b
end

"""
    Yields.AbstractYieldModel - Yields.AbstractYieldModel

The subtraction of two yields will create a `CompositeYield`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the second curves will be subtracted from the first.
"""
function Base.:-(a::AbstractYieldModel, b::AbstractYieldModel)
    return CompositeYield(a, b, -)
end

function Base.:-(a::Constant, b::Constant)
    Constant(a.rate - b.rate)
end

function Base.:-(a::T, b) where {T<:AbstractYieldModel}
    return a - Constant(b)
end

function Base.:-(a, b::T) where {T<:AbstractYieldModel}
    return Constant(a) - b
end

"""
    Yields.AbstractYieldModel / Yields.AbstractYieldModel

The division of two yields will create a `CompositeYield`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will have the first divided by the second. This can be useful, for example, if you wanted to gross-up a yield to be pre-tax.

# Examples

```julia-repl
julia> m = Yields.Constant(0.01) / 0.79;

julia> accumulation(d,1)
1.0126582278481013

julia> accumulation(.01/.79,1)
1.0126582278481013
```
"""
function Base.:/(a::AbstractYieldModel, b::AbstractYieldModel)
    return CompositeYield(a, b, /)
end

function Base.:/(a::Constant, b::Constant)
    a_kind = rate(a).compounding
    rate_new_basis = rate(convert(a_kind, rate(b)))
    return Constant(
        Rate(
            rate(a.rate) / rate_new_basis,
            a_kind
        )
    )
end

function Base.:/(a::T, b) where {T<:AbstractYieldModel}
    return a / Constant(b)
end

function Base.:/(a, b::T) where {T<:AbstractYieldModel}
    return Constant(a) / b
end

end
