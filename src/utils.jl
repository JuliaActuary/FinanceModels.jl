# all we need is cumulative normal, so avoid Distributions.jl dependency
# https://www.johndcook.com/blog/cpp_phi/

function ϕ(x)
    return 0.5 * SpecialFunctions.erfc(-x * √(0.5))
end

N(x) = ϕ(x)

function d1(S, K, τ, r, σ, q)
    return (log(S / K) + (r - q + σ^2 / 2) * τ) / (σ * √(τ))
end

function d2(S, K, τ, r, σ, q)
    return d1(S, K, τ, r, σ, q) - σ * √(τ)
end

"""
    eurocall(;S=1.,K=1.,τ=1,r,σ,q=0.)

Calculate the Black-Scholes implied option price for a european call, where:

- `S` is the current asset price
- `K` is the strike or exercise price
- `τ` is the time remaining to maturity (can be typed with \\tau[tab])
- `r` is the continuously compounded risk free rate
- `σ` is the (implied) volatility (can be typed with \\sigma[tab])
- `q` is the continuously paid dividend rate

Rates should be input as rates (not percentages), e.g.: `0.05` instead of `5` for a rate of five percent.

!!! warning "Experimental"
    This function is well-tested, but the derivatives functionality (API) may change in a future version of FinanceModels.

# Extended Help

This is the same as the formulation presented in the [dividend extension of the BS model in Wikipedia](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model#Black%E2%80%93Scholes_equation).

## Other general comments:

- Swap/OIS curves are generally better sources for `r` than government debt (e.g. US Treasury) due to the collateralized nature of swap instruments.
- (Implied) volatility is characterized by a curve that is a function of the strike price (among other things), so take care when using
- FinanceModels.jl can assist with converting rates to continuously compounded if you need to perform conversions (e.g. `convert(Continuous(), r)`).

"""
function eurocall(; S = 1.0, K = 1.0, τ = 1, r, σ, q = 0.0)
    iszero(τ) && return max(zero(S), S - K)
    d₁ = d1(S, K, τ, r, σ, q)
    d₂ = d2(S, K, τ, r, σ, q)
    return (N(d₁) * S * exp(τ * (r - q)) - N(d₂) * K) * exp(-r * τ)
end

"""
    europut(;S=1.,K=1.,τ=1,r,σ,q=0.)

Calculate the Black-Scholes implied option price for a european put, where:

- `S` is the current asset price
- `K` is the strike or exercise price
- `τ` is the time remaining to maturity (can be typed with \\tau[tab])
- `r` is the continuously compounded risk free rate
- `σ` is the (implied) volatility (can be typed with \\sigma[tab])
- `q` is the continuously paid dividend rate

Rates should be input as rates (not percentages), e.g.: `0.05` instead of `5` for a rate of five percent.


!!! warning "Experimental"
    This function is well-tested, but the derivatives functionality (API) may change in a future version of FinanceModels.

# Extended Help

This is the same as the formulation presented in the [dividend extension of the BS model in Wikipedia](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model#Black%E2%80%93Scholes_equation).

## Other general comments:

- Swap/OIS curves are generally better sources for `r` than government debt (e.g. US Treasury) due to the collateralized nature of swap instruments.
- (Implied) volatility is characterized by a curve that is a function of the strike price (among other things), so take care when using
- FinanceModels.jl can assist with converting rates to continuously compounded if you need to perform conversions (e.g. `convert(Continuous(), r)`).

"""
function europut(; S = 1.0, K = 1.0, τ = 1, r, σ, q = 0.0)
    iszero(τ) && return max(zero(S), K - S)
    d₁ = d1(S, K, τ, r, σ, q)
    d₂ = d2(S, K, τ, r, σ, q)
    return (N(-d₂) * K - N(-d₁) * S * exp(τ * (r - q))) * exp(-r * τ)
end

"""
    ReadOnlyVector(v::Vector)

Internal read-only view over an owned `Vector`. It defines no `setindex!`, so indexed assignment
— and therefore `.=`, `fill!`, `sort!`, `reverse!`, and writes through `view` — throws Base's
`CanonicalIndexError`; reads behave as a normal `AbstractVector` (indexing, iteration,
`searchsortedlast`, broadcasting, `==`/`isequal`/`hash` identical to the equivalent `Vector`).
`copy`/`collect` return a mutable `Vector`. To change a curve's knots, use
`Accessors.@set curve.rates[i] = x` or `reconstruct`, which rebuild it.

Used by the knot curves (`Yield.Spline`, `Yield.MonotoneConvex`), which cache state derived
from their knot vectors, so that ordinary array operations on the public fields cannot
desynchronise the cache. This is Julia's conventional privacy, not literal immutability: the
backing `Vector` is the internal field `_data`, and code that mutates it is unsupported.
"""
struct ReadOnlyVector{T} <: AbstractVector{T}
    _data::Vector{T}
end

Base.size(v::ReadOnlyVector) = size(getfield(v, :_data))
Base.IndexStyle(::Type{<:ReadOnlyVector}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(v::ReadOnlyVector, i::Int) = getfield(v, :_data)[i]

# The polynomial with exact (Rational) coefficients `cs` at `x`, each coefficient rounded once to
# the precision of `x` (for a dual number, of its primal value): Float32 stays Float32 and BigFloat
# keeps its precision. Inlined, so that for a concrete float type the rounding of a constant table
# folds at compile time.
@inline __evalpoly_exact(x, cs) = evalpoly(x, map(c -> convert(typeof(float(__primal(x))), c), cs))

# The float element type of a collection of reals, for a copy or an accumulator that keeps the values'
# numeric type (BigFloat, dual numbers): the declared element type when concrete, otherwise the
# promotion of the values' types; `Bool`, which every real type absorbs, for an untyped empty
# collection, so that it is `Float64`.
__float_eltype(v) = float(isconcretetype(eltype(v)) ? eltype(v) : mapreduce(typeof, promote_type, v; init = Bool))
