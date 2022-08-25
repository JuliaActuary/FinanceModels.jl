## Curve Manipulations
"""
    RateCombination(curve1,curve2,operation)

Creates a datastructure that will perform the given `operation` after independently calculating the effects of the two curves. 
Can only be created via the public API by using the `+`, `-`, `*`, and `/` operatations on `AbstractYield` objects.

As this is double the normal operations when performing calculations, if you are using the curve in performance critical locations, you should consider transforming the inputs and 
constructing a single curve object ahead of time.
"""
struct RateCombination{T,U,V} <: AbstractYieldCurve
    r1::T
    r2::U
    op::V
end
__ratetype(::Type{RateCombination{T,U,V}}) where {T,U,V}= __ratetype(T)

FinanceCore.rate(rc::RateCombination, time) = rc.op(rate(rc.r1, time), rate(rc.r2, time))
function FinanceCore.discount(rc::RateCombination, time)
    a1 = discount(rc.r1, time)^(-1 / time) - 1
    a2 = discount(rc.r2, time)^(-1 / time) - 1
    return 1 / (1 + rc.op(a1, a2))^time
end

Base.zero(rc::RateCombination, time) = zero(rc,time,Periodic(1))
function Base.zero(rc::RateCombination, time, cf::C) where {C<:FinanceCore.CompoundingFrequency}
    d = discount(rc,time)
    i = Periodic(1/d^(1/time)-1,1)
    return convert(cf, i) # c.zero is a curve of continuous rates represented as floats. explicitly wrap in continuous before converting
end

"""
    Yields.AbstractYieldCurve + Yields.AbstractYieldCurve

The addition of two yields will create a `RateCombination`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will be added together.
"""
function Base.:+(a::AbstractYieldCurve, b::AbstractYieldCurve)
    return RateCombination(a, b, +)
end

function Base.:+(a::Constant, b::Constant)
    a_kind = rate(a).compounding
    rate_new_basis = rate(convert(a_kind, rate(b)))
    return Constant(
        Rate(
            rate(a.rate) + rate_new_basis,
            a_kind
        )
    )
end

function Base.:+(a::T, b) where {T<:AbstractYieldCurve}
    return a + Constant(b)
end

function Base.:+(a, b::T) where {T<:AbstractYieldCurve}
    return Constant(a) + b
end

"""
    Yields.AbstractYieldCurve * Yields.AbstractYieldCurve

The multiplication of two yields will create a `RateCombination`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will be added together. This can be useful, for example, if you wanted to after-tax a yield.

# Examples

```julia-repl
julia> m = Yields.Constant(0.01) * 0.79;

julia> accumulation(m,1)
1.0079

julia> accumulation(.01*.79,1)
1.0079
```
"""
function Base.:*(a::AbstractYieldCurve, b::AbstractYieldCurve)
    return RateCombination(a, b, *)
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

function Base.:*(a::T, b) where {T<:AbstractYieldCurve}
    return a * Constant(b)
end

function Base.:*(a, b::T) where {T<:AbstractYieldCurve}
    return Constant(a) * b
end

"""
    Yields.AbstractYieldCurve - Yields.AbstractYieldCurve

The subtraction of two yields will create a `RateCombination`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the second curves will be subtracted from the first.
"""
function Base.:-(a::AbstractYieldCurve, b::AbstractYieldCurve)
    return RateCombination(a, b, -)
end

function Base.:-(a::Constant, b::Constant)
    a_kind = rate(a).compounding
    rate_new_basis = rate(convert(a_kind, rate(b)))
    return Constant(
        Rate(
            rate(a.rate) - rate_new_basis,
            a_kind
        )
    )
end

function Base.:-(a::T, b) where {T<:AbstractYieldCurve}
    return a - Constant(b)
end

function Base.:-(a, b::T) where {T<:AbstractYieldCurve}
    return Constant(a) - b
end

"""
    Yields.AbstractYieldCurve / Yields.AbstractYieldCurve

The division of two yields will create a `RateCombination`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will have the first divided by the second. This can be useful, for example, if you wanted to gross-up a yield to be pre-tax.

# Examples

```julia-repl
julia> m = Yields.Constant(0.01) / 0.79;

julia> accumulation(d,1)
1.0126582278481013

julia> accumulation(.01/.79,1)
1.0126582278481013
```
"""
function Base.:/(a::AbstractYieldCurve, b::AbstractYieldCurve)
    return RateCombination(a, b, /)
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

function Base.:/(a::T, b) where {T<:AbstractYieldCurve}
    return a / Constant(b)
end

function Base.:/(a, b::T) where {T<:AbstractYieldCurve}
    return Constant(a) / b
end