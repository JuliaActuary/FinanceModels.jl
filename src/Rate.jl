
abstract type CompoundingFrequency end
Base.Broadcast.broadcastable(x::T) where {T<:CompoundingFrequency} = Ref(x)

""" 
    Continuous()

A type representing continuous interest compounding frequency.

# Examples

```julia-repl
julia> Rate(0.01,Continuous())
Rate(0.01, Continuous())
```

See also: [`Periodic`](@ref)
"""
struct Continuous <: CompoundingFrequency end


""" 


```julia-repl
julia> Continuous(0.01)
Rate(0.01, Continuous())
```

See also: [`Periodic`](@ref)
"""
Continuous(rate) = Rate(rate, Continuous())

""" 
    Periodic(frequency)

A type representing periodic interest compounding with the given frequency

# Examples

Creating a semi-annual bond equivalent yield:

```julia-repl
julia> Rate(0.01,Periodic(2))
Rate(0.01, Periodic(2))
```

See also: [`Continuous`](@ref)
"""
struct Periodic <: CompoundingFrequency
    frequency::Int
end

""" 
    Periodic(rate,frequency)

A convenience constructor for Rate(rate,Periodic(frequency)).

# Examples

Creating a semi-annual bond equivalent yield:

```julia-repl
julia> Periodic(0.01,2)
Rate(0.01, Periodic(2))
```

See also: [`Continuous`](@ref)
"""
Periodic(x, frequency) = Rate(x, Periodic(frequency))

struct Rate{N<:Real,T<:CompoundingFrequency}
    value::N
    compounding::T
end

# Base.:==(r1::Rate,r2::Rate) = (r1.value == r2.value) && (r1.compounding == r2.compounding)

"""
    Rate(rate[,frequency=1])
    Rate(rate,frequency::CompoundingFrequency)

Rate is a type that encapsulates an interest `rate` along with its compounding `frequency`.

Periodic rates can be constructed via `Rate(rate,frequency)` or `Rate(rate,Periodic(frequency))`.

Continuous rates can be constructed via `Rate(rate, Inf)` or `Rate(rate,Continuous())`.

# Examples

```julia-repl
julia> Rate(0.01,Continuous())
Rate(0.01, Continuous())

julia> Rate(0.01,Periodic(2))
Rate(0.01, Periodic(2))

julia> Rate(0.01)
Rate(0.01, Periodic(1))

julia> Rate(0.01,2)
Rate(0.01, Periodic(2))

julia> Rate(0.01,Periodic(4))
Rate(0.01, Periodic(4))

julia> Rate(0.01,Inf)
Rate(0.01, Continuous())

julia> Rate(0.01,Continuous())
Rate(0.01, Continuous())
```
"""
Rate(rate) = Rate(rate, Periodic(1))
Rate(x, frequency::T) where {T<:Real} = isinf(frequency) ? Rate(x, Continuous()) : Rate(x, Periodic(frequency))

"""
    convert(T::CompoundingFrequency,r::Rate)

Returns a `Rate` with an equivalent discount but represented with a different compounding frequency.

# Examples

```julia-repl
julia> r = Rate(Periodic(12),0.01)
Rate(0.01, Periodic(12))

julia> convert(Periodic(1),r)
Rate(0.010045960887181016, Periodic(1))

julia> convert(Continuous(),r)
Rate(0.009995835646701251, Continuous())
```
"""
function Base.convert(T::CompoundingFrequency, r::Rate{<:Real,<:CompoundingFrequency})
    convert(T, r, r.compounding)
end
function Base.convert(to::Continuous, r, from::Continuous)
    return r
end

function Base.convert(to::Periodic, r, from::Continuous)
    return Rate(to.frequency * (exp(r.value / to.frequency) - 1), to)
end

function Base.convert(to::Continuous, r, from::Periodic)
    return Rate(from.frequency * log(1 + r.value / from.frequency), to)
end

function Base.convert(to::Periodic, r, from::Periodic)
    c = convert(Continuous(), r, from)
    return convert(to, c, Continuous())
end

function rate(r::Rate{<:Real,<:CompoundingFrequency})
    r.value
end

function Base.isapprox(a::Rate{N,T}, b::Rate{N,T}; atol::Real = 0, rtol::Real = atol > 0 ? 0 : √eps()) where {T<:Periodic,N<:Real}
    return (a.compounding.frequency == b.compounding.frequency) && isapprox(rate(a), rate(b); atol, rtol)
end

function Base.isapprox(a::Rate{N,T}, b::Rate{N,T}; atol::Real = 0, rtol::Real = atol > 0 ? 0 : √eps()) where {T<:Continuous,N<:Real}
    return isapprox(rate(a), rate(b); atol, rtol)
end

# the fallback for rates not of the same type
function Base.isapprox(a::T, b::N; atol::Real = 0, rtol::Real = atol > 0 ? 0 : √eps()) where {T<:Rate,N<:Rate}
    return isapprox(convert(b.compounding, a), b; atol, rtol)
end