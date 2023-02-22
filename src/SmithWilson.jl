@Base.kwdef struct SmithWilson{U,A} <: CurveMethod
    ufr::U
    α::A
end

abstract type ObservableQuote end

"""
    ZeroCouponQuote(price, maturity)

Quote for a set of zero coupon bonds with given `price` and `maturity`. 

# Examples

```julia-repl
julia> prices = [1.3, 0.1, 4.5]
julia> maturities = [1.2, 2.5, 3.6]
julia> swq = Yields.ZeroCouponQuote.(prices, maturities)
```
"""
struct ZeroCouponQuote <: ObservableQuote
    price
    maturity
end

"""
    SwapQuote(yield, maturity, frequency)

Quote for a set of interest rate swaps with the given `yield` and `maturity` and a given payment `frequency`.

# Examples

```julia-repl
julia> maturities = [1.2, 2.5, 3.6]
julia> interests = [-0.02, 0.3, 0.04]
julia> prices = [1.3, 0.1, 4.5]
julia> frequencies = [2,1,2]
julia> swq = Yields.SwapQuote.(interests, maturities, frequencies)
```
"""
struct SwapQuote <: ObservableQuote
    yield
    maturity
    frequency
    function SwapQuote(yield, maturity, frequency)
        frequency <= 0 && throw(DomainError("Payment frequency must be positive"))
        return new(yield, maturity, frequency)
    end
end


"""
    BulletBondQuote(yield, price, maturity, frequency)

Quote for a set of fixed interest bullet bonds with given `yield`, `price`, `maturity` and a given payment frequency `frequency`.

Construct a vector of quotes for use with SmithWilson methods, e.g. by broadcasting over an array of inputs.

# Examples

```julia-repl
julia> maturities = [1.2, 2.5, 3.6]
julia> interests = [-0.02, 0.3, 0.04]
julia> prices = [1.3, 0.1, 4.5]
julia> frequencies = [2,1,2]
julia> bbq = Yields.BulletBondQuote.(interests, maturities, prices, frequencies)
```
"""
struct BulletBondQuote <: ObservableQuote
    yield
    price
    maturity
    frequency

    function BulletBondQuote(yield, maturity, price, frequency)
        frequency <= 0 && throw(DomainError("Payment frequency must be positive"))
        return new(yield, maturity, price, frequency)
    end
end


"""
    SmithWilson(zcq::Vector{ZeroCouponQuote}; ufr, α)
    SmithWilson(swq::Vector{SwapQuote}; ufr, α)
    SmithWilson(bbq::Vector{BulletBondQuote}; ufr, α)
    SmithWilson(times<:AbstractVector, cashflows<:AbstractMatrix, prices<:AbstractVector; ufr, α)
    SmithWilson(u, qb; ufr, α)
    
Create a yield curve object that implements the Smith-Wilson interpolation/extrapolation scheme.

Positional arguments to construct a curve:

- Quoted instrument as the first argument: either a `Vector` of `ZeroCouponQuote`s, `SwapQuote`s, or `BulletBondQuote`s, or 
- A set of `times`, `cashflows`, and `prices`, or
- A curve can be with `u` is the timepoints coming from the calibration, and `qb` is the internal parameterization of the curve that ensures that the calibration is correct. Users may prefer the other constructors but this mathematical constructor is also available.

Required keyword arguments:

- `ufr` is the Ultimate Forward Rate, the forward interest rate to which the yield curve tends, in continuous compounding convention. 
- `α` is the parameter that governs the speed of convergence towards the Ultimate Forward Rate. It can be typed with `\\alpha[TAB]`
"""
struct SmithWilsonCurve{TU<:AbstractVector,TQb<:AbstractVector} <: AbstractYieldCurve
    u::TU
    qb::TQb
    ufr
    α

    # Inner constructor ensures that vector lengths match
    function SmithWilsonCurve{TU,TQb}(u, qb; ufr, α) where {TU<:AbstractVector,TQb<:AbstractVector}
        if length(u) != length(qb)
            throw(DomainError("Vectors u and qb in SmithWilson must have equal length"))
        end
        return new(u, qb, ufr, α)
    end
end

# SmithWilson(u::TU, qb::TQb; ufr, α) where {TU<:AbstractVector,TQb<:AbstractVector} = SmithWilson{TU,TQb}(u, qb; ufr = ufr, α = α)

__ratetype(::Type{SmithWilsonCurve{TU,TQb}}) where {TU,TQb}= Yields.Rate{Float64, Yields.Continuous}

"""
    H_ordered(α, t_min, t_max)

The Smith-Wilson H function with ordered arguments (for better performance than using min and max).
"""
function H_ordered(α, t_min, t_max)
    return α * t_min + exp(-α * t_max) * sinh(-α * t_min)
end

"""
    H(α, t1, t2)

The Smith-Wilson H function implemented in a faster way.
"""
function H(α, t1::T, t2::T) where {T}
    return t1 < t2 ? H_ordered(α, t1, t2) : H_ordered(α, t2, t1)
end

H(α, t1, t2) = H(α, promote(t1, t2)...)

H(α, t1vec::AbstractVector, t2) = [H(α, t1, t2) for t1 in t1vec]
H(α, t1vec::AbstractVector, t2vec::AbstractVector) = [H(α, t1, t2) for t1 in t1vec, t2 in t2vec]
# This can be optimized by going to H_ordered directly, but it might be a bit cumbersome 
H(α, tvec::AbstractVector) = H(α, tvec, tvec)


FinanceCore.discount(sw::SmithWilsonCurve, t) = exp(-sw.ufr * t) * (1.0 + H(sw.α, sw.u, t) ⋅ sw.qb)
Base.zero(sw::SmithWilsonCurve, t) = Continuous(sw.ufr - log(1.0 + H(sw.α, sw.u, t) ⋅ sw.qb) / t)

function SmithWilson(times::AbstractVector, cashflows::AbstractMatrix, prices::AbstractVector; ufr, α)
    Q = Diagonal(exp.(-ufr * times)) * cashflows
    q = vec(sum(Q, dims = 1))  # We want q to be a column vector
    QHQ = Q' * H(α, times) * Q
    b = QHQ \ (prices - q)
    Qb = Q * b
    return SmithWilson(times, Qb; ufr = ufr, α = α)
end

function cashflows(qs::Vector{Q}) where {Q<:ObservableQuote}
    yield = [q.yield for q in qs]
    maturity = [q.maturity for q in qs]
    frequency = [q.frequency for q in qs]
    return cashflows(yield, maturity, frequency)
end



function curve() 
end

