
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

# make rate a broadcastable type
Base.Broadcast.broadcastable(ic::T) where {T<:Rate} = Ref(ic)

"""
    Rate(rate[,frequency=1])
    Rate(rate,frequency::CompoundingFrequency)

Rate is a type that encapsulates an interest `rate` along with its compounding `frequency`.

Periodic rates can be constructed via `Rate(rate,frequency)` or `Rate(rate,Periodic(frequency))`. If not given a second argument, `Rate(rate)` is equivalent to `Rate(rate,Periodic(1))`.

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
    c = convert(a.compounding,b)
    return (a.compounding.frequency == c.compounding.frequency) && isapprox(rate(a), rate(c); atol, rtol)
end

function Base.isapprox(a::Rate{N,T}, b::Rate{N,T}; atol::Real = 0, rtol::Real = atol > 0 ? 0 : √eps()) where {T<:Continuous,N<:Real}
    return isapprox(rate(a), rate(b); atol, rtol)
end

# the fallback for rates not of the same type
function Base.isapprox(a::T, b::N; atol::Real = 0, rtol::Real = atol > 0 ? 0 : √eps()) where {T<:Rate,N<:Rate}
    return isapprox(convert(b.compounding, a), b; atol, rtol)
end

"""
    discount(rate, t)
    discount(rate, from, to)

Discount `rate` for a time `t` or for an interval `(from, to)`. If `rate` is not a `Rate`, it will be assumed to be a `Periodic` rate compounded once per period, i.e. `Periodic(rate,1)`. 

# Examples

```julia-repl
julia> discount(0.03, 10)
0.7440939148967249

julia> discount(Periodic(0.03, 2), 10)
0.7424704182237725

julia> discount(Continuous(0.03), 10)
0.7408182206817179

julia> discount(0.03, 5, 10)
0.8626087843841639
```
"""

discount(rate::T, t) where {T<:Real} = discount(Rate(rate), t)
discount(rate::Rate{<:Real, <:Continuous}, t) = exp(-rate.value * t)
discount(rate::Rate{<:Real, <:Periodic}, t) = (1 + rate.value / rate.compounding.frequency)^(-rate.compounding.frequency * t)
discount(rate, from, to) = discount(rate, to - from)

"""
    accumulation(rate, t)
    accumulation(rate, from, to)

Accumulate `rate` for a time `t` or for an interval `(from, to)`. If `rate` is not a `Rate`, it will be assumed to be a `Periodic` rate compounded once per period, i.e. `Periodic(rate,1)`. 

    # Examples

```julia-repl
julia> accumulation(0.03, 10)
1.3439163793441222

julia> accumulation(Periodic(0.03, 2), 10)
1.3468550065500535

julia> accumulation(Continuous(0.03), 10)
1.3498588075760032

julia> accumulation(0.03, 5, 10)
1.1592740743
```
"""

accumulation(rate, t) = accumulation(Rate(rate), t)
accumulation(rate::Rate{<:Real, <:Continuous}, t) = exp(rate.value * t)
accumulation(rate::Rate{<:Real, <:Periodic}, t) = (1 + rate.value / rate.compounding.frequency)^(rate.compounding.frequency * t)
accumulation(rate, from, to) = accumulation(rate, to - from)


"""
    +(Yields.Rate, T<:Real)
    +(T<:Real, Yields.Rate)
    +(Yields.Rate,Yields.Rate)

The addition of a rate with a number will inherit the type of the `Rate`, or the first argument's type if both are `Rate`s.

# Examples

```julia-repl
julia> Yields.Periodic(0.01,2) + Yields.Periodic(0.04,2)
Yields.Rate{Float64, Yields.Periodic}(0.05000000000000004, Yields.Periodic(2))

julia> Yields.Periodic(0.04,2) + 0.01
Yields.Rate{Float64, Yields.Periodic}(0.05, Yields.Periodic(2))
```
"""
function Base.:+(a::Rate{N,T}, b::Real) where {N<:Real,T<:Continuous}
    return Continuous(a.value + b)
end
function Base.:+(a::Real, b::Rate{N,T}) where {N<:Real,T<:Continuous}
    return Continuous(b.value + a)
end

function Base.:+(a::Rate{N,T}, b::Real) where {N<:Real, T<:Periodic}
    return Periodic(a.value + b,a.compounding.frequency)
end
function Base.:+(a::Real, b::Rate{N,T}) where {N<:Real, T<:Periodic}
    return Periodic(b.value + a,b.compounding.frequency)
end

function Base.:+(a::Rate{N,T},b::Rate{N,T}) where {N<:Real,T<:CompoundingFrequency}
    a_rate = rate(a)
    b_rate = rate(convert(a.compounding,b))
    r = Rate(a_rate + b_rate,a.compounding)
    return r
end

"""
    -(Yields.Rate, T<:Real)
    -(T<:Real, Yields.Rate)
    -(Yields.Rate, Yields.Rate)


The addition of a rate with a number will inherit the type of the `Rate`, or the first argument's type if both are `Rate`s.

# Examples

```julia-repl
julia> Yields.Periodic(0.04,2) - Yields.Periodic(0.01,2)
Yields.Rate{Float64, Yields.Periodic}(0.030000000000000214, Yields.Periodic(2))

julia> Yields.Periodic(0.04,2) - 0.01
Yields.Rate{Float64, Yields.Periodic}(0.03, Yields.Periodic(2))

```
"""
function Base.:-(a::Rate{N,T}, b::Real) where {N<:Real,T<:Continuous}
    return Continuous(a.value - b)
end
function Base.:-(a::Real, b::Rate{N,T}) where {N<:Real,T<:Continuous}
    return Continuous( a - b.value)
end

function Base.:-(a::Rate{N,T}, b::Real) where {N<:Real, T<:Periodic}
    return Periodic(a.value - b, a.compounding.frequency)
end
function Base.:-(a::Real, b::Rate{N,T}) where {N<:Real, T<:Periodic}
    return Periodic(a - b.value, b.compounding.frequency)
end
function Base.:-(a::Rate{N,T},b::Rate{N,T}) where {N<:Real,T<:CompoundingFrequency}
    a_rate = rate(a)
    b_rate = rate(convert(a.compounding,b))
    r = Rate(a_rate - b_rate,a.compounding)
    return r
end

"""
    *(Yields.Rate, T<:Real)
    *(T<:Real, Yields.Rate)

The multiplication of a Rate with a scalar will inherit the type of the `Rate`, or the first argument's type if both are `Rate`s.
"""
function Base.:*(a::Rate{N,T}, b::Real) where {N<:Real,T<:Continuous}
    return Continuous(a.value * b)
end
function Base.:*(a::Real, b::Rate{N,T}) where {N<:Real,T<:Continuous}
    return Continuous( a * b.value)
end

function Base.:*(a::Rate{N,T}, b::Real) where {N<:Real, T<:Periodic}
    return Periodic(a.value * b, a.compounding.frequency)
end
function Base.:*(a::Real, b::Rate{N,T}) where {N<:Real, T<:Periodic}
    return Periodic(a * b.value, b.compounding.frequency)
end


"""
    /(Yields.Rate, T<:Real)
    /(T<:Real, Yields.Rate)

The division of a Rate with a scalar will inherit the type of the `Rate`, or the first argument's type if both are `Rate`s.
"""
function Base.:/(a::Rate{N,T}, b::Real) where {N<:Real,T<:Continuous}
    return Continuous(a.value / b)
end
function Base.:/(a::Real, b::Rate{N,T}) where {N<:Real,T<:Continuous}
    return Continuous( a / b.value)
end

function Base.:/(a::Rate{N,T}, b::Real) where {N<:Real, T<:Periodic}
    return Periodic(a.value / b, a.compounding.frequency)
end
function Base.:/(a::Real, b::Rate{N,T}) where {N<:Real, T<:Periodic}
    return Periodic(a / b.value, b.compounding.frequency)
end


"""
    <(Rate,Rate)

Convert the second argument to the periodicity of the first and compare the scalar rate values to determine if the first argument has a lower force of interest than the second.

# Examples

```julia-repl
julia> Yields.Periodic(0.03,100) < Yields.Continuous(0.03)
true
```
"""
function Base.:<(a::Rate{N,T},b::Rate{N,U}) where {N<:Real,T<:CompoundingFrequency,U<:CompoundingFrequency}
    bc = convert(a.compounding,b)
    return rate(a) < rate(bc)
end

"""
    >(Rate,Rate)

Convert the second argument to the periodicity of the first and compare the scalar rate values to determine if the first argument has a greater force of interest than the second.

# Examples

```julia-repl
julia> Yields.Periodic(0.03,100) > Yields.Continuous(0.03)
false
```
"""
function Base.:>(a::Rate{N,T},b::Rate{N,U}) where {N<:Real,T<:CompoundingFrequency,U<:CompoundingFrequency}
    bc = convert(a.compounding,b)
    return rate(a) > rate(bc)
end