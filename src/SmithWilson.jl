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
struct SmithWilson{TU<:AbstractVector,TQb<:AbstractVector} <: AbstractYield
    u::TU
    qb::TQb
    ufr
    α

    # Inner constructor ensures that vector lengths match
    function SmithWilson{TU,TQb}(u, qb; ufr, α) where {TU<:AbstractVector,TQb<:AbstractVector}
        if length(u) != length(qb)
            throw(DomainError("Vectors u and qb in SmithWilson must have equal length"))
        end
        return new(u, qb, ufr, α)
    end
end

SmithWilson(u::TU, qb::TQb; ufr, α) where {TU<:AbstractVector,TQb<:AbstractVector} = SmithWilson{TU,TQb}(u, qb; ufr = ufr, α = α)

__ratetype(::SmithWilson{TU,TQb}) where {TU,TQb}= Yields.Rate{Float64, Yields.Continuous}

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


discount(sw::SmithWilson, t) = exp(-sw.ufr * t) * (1.0 + H(sw.α, sw.u, t) ⋅ sw.qb)
Base.zero(sw::SmithWilson, t) = Continuous(sw.ufr - log(1.0 + H(sw.α, sw.u, t) ⋅ sw.qb) / t)
Base.zero(sw::SmithWilson, t, cf::CompoundingFrequency) = convert(cf, zero(sw, t))

function SmithWilson(times::AbstractVector, cashflows::AbstractMatrix, prices::AbstractVector; ufr, α)
    Q = Diagonal(exp.(-ufr * times)) * cashflows
    q = vec(sum(Q, dims = 1))  # We want q to be a column vector
    QHQ = Q' * H(α, times) * Q
    b = QHQ \ (prices - q)
    Qb = Q * b
    return SmithWilson(times, Qb; ufr = ufr, α = α)
end

""" 
    timepoints(zcq::Vector{ZeroCouponQuote})
    timepoints(bbq::Vector{BulletBondQuote})

Return the times associated with the `cashflows` of the instruments.
"""
function timepoints(qs::Vector{Q}) where {Q<:ObservableQuote}
    frequency = maximum(q.frequency for q in qs)
    timestep = 1 / frequency
    maturity = maximum(q.maturity for q in qs)
    return [timestep:timestep:maturity...]
end


"""
    cashflows(interests, maturities, frequency)
    timepoints(zcq::Vector{ZeroCouponQuote})
    timepoints(bbq::Vector{BulletBondQuote})

Produce a cash flow matrix for a set of instruments with given `interests` and `maturities`
and a given payment frequency `frequency`. All instruments are assumed to have their first payment at time 1/`frequency`
and have their last payment at the largest multiple of 1/`frequency` less than or equal to the input maturity.
"""
function cashflows(interests, maturities, frequencies)
    frequency = lcm(frequencies)
    fq = inv.(frequencies)
    timestep = 1 / frequency
    floored_mats = floor.(maturities ./ timestep) .* timestep
    times = timestep:timestep:maximum(floored_mats)
    # we need to determine the coupons in relation to the payment date, not time zero
    time_adj = floored_mats .% fq

    cashflows = [
        # if on a coupon date and less than maturity, pay coupon
        ((((t + time_adj[instrument]) % fq[instrument] ≈ 0) && t <= floored_mats[instrument]) ? interests[instrument] / frequencies[instrument] : 0.0) +
        (t ≈ floored_mats[instrument] ? 1.0 : 0.0) # add maturity payment
        for t in times, instrument = 1:length(interests)
    ]

    return cashflows
end

function cashflows(qs::Vector{Q}) where {Q<:ObservableQuote}
    yield = [q.yield for q in qs]
    maturity = [q.maturity for q in qs]
    frequency = [q.frequency for q in qs]
    return cashflows(yield, maturity, frequency)
end

# Utility methods for calibrating Smith-Wilson directly from quotes
function SmithWilson(zcq::Vector{ZeroCouponQuote}; ufr, α)
    n = length(zcq)
    maturities = [q.maturity for q in zcq]
    prices = [q.price for q in zcq]
    return SmithWilson(maturities, Matrix{Float64}(I, n, n), prices; ufr = ufr, α = α)
end

function SmithWilson(swq::Vector{SwapQuote}; ufr, α)
    times = timepoints(swq)
    cfs = cashflows(swq)
    ones(length(swq))
    return SmithWilson(times, cfs, ones(length(swq)), ufr = ufr, α = α)
end

function SmithWilson(bbq::Vector{BulletBondQuote}; ufr, α)
    times = timepoints(bbq)
    cfs = cashflows(bbq)
    prices = [q.price for q in bbq]
    return SmithWilson(times, cfs, prices, ufr = ufr, α = α)
end