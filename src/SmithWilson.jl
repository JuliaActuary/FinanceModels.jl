@Base.kwdef struct SmithWilson{U,A} <: CurveMethod
    ufr::U
    α::A
end

"""
    SmithWilsonCurve(zcq::Vector{ZeroCouponQuote}; ufr, α)
    SmithWilsonCurve(swq::Vector{SwapQuote}; ufr, α)
    SmithWilsonCurve(bbq::Vector{BulletBondQuote}; ufr, α)
    SmithWilsonCurve(times<:AbstractVector, cashflows<:AbstractMatrix, prices<:AbstractVector; ufr, α)
    SmithWilsonCurve(u, qb; ufr, α)
    
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
    function SmithWilsonCurve(sw::S, u::TU, qb::TQb) where {S<:SmithWilson,TU<:AbstractVector,TQb<:AbstractVector}
        if length(u) != length(qb)
            throw(DomainError("Vectors u and qb in SmithWilson must have equal length"))
        end
        return new{TU,TQb}(u, qb, sw.ufr, sw.α)
    end
end

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

function __SW_inner(sw::SmithWilson,times, cashflows, prices)
    Q = Diagonal(exp.(-sw.ufr * times)) * cashflows
    q = vec(sum(Q, dims = 1))  # We want q to be a column vector
    QHQ = Q' * H(sw.α, times) * Q
    b = QHQ \ (prices - q)
    Qb = Q * b
    return SmithWilsonCurve(sw, times, Qb)
end

function (sw::SmithWilson)(quotes::Vector{Quote{T,C}}) where {T,C<:Cashflow}
    n = length(quotes)
    maturities = [q.instrument.time for q in quotes]
    prices = [q.price for q in quotes]
    return __SW_inner(sw,maturities, Matrix{Float64}(I, n, n), prices)
end

function (sw::SmithWilson)(quotes::Vector{Quote{T,C}}) where {T,C<:Cashflow}
    n = length(quotes)
    time
    prices = [q.price for q in quotes]
    return __SW_inner(sw,maturities, Matrix{Float64}(I, n, n), prices)
end

function (sw::SmithWilson)(quotes::Vector{Quote{T,B}}) where {T,B<:Bond}
    @show ts = timesteps(quotes)
    @show prices = [q.price for q in quotes]
    @show cfs = cashflow_matrix(quotes)
    return __SW_inner(sw,ts, cfs, prices)
end