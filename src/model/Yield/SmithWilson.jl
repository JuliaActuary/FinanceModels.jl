# SmithWilson
## Originally developed by kasperrisager


using ..LinearAlgebra
using ..FinanceCore

"""
    Yield.SmithWilson(u, qb; ufr=ufr, α=α)
    Yield.SmithWilson(;ufr=ufr, α=α)

    
Create a yield curve object that implements the Smith-Wilson interpolation/extrapolation scheme.

To calibrate a curve, you generally want to construct the object without the `u` and `qb` arguments and call [`fit`](@ref FinanceModels.fit) in conjunction with Quotes (`fit` requires no third parameter for SmithWilson curves). See **Examples** for what this looks like.
Positional arguments to construct a curve:
- A curve can be with `u` is the timepoints coming from the calibration, and `qb` is the internal parameterization of the curve that ensures that the calibration is correct. Users may prefer the other constructors but this mathematical constructor is also available.

Required keyword arguments:

- `ufr` is the Ultimate Forward Rate, the forward interest rate to which the yield curve tends: a number, read as continuously compounded, or a `Rate`, converted to its continuously compounded value. The `ufr` field holds that number.
- `α` is the parameter that governs the speed of convergence towards the Ultimate Forward Rate. It can be typed with `\\alpha[TAB]`

# Examples

```julia
times = [1.0, 2.5, 5.6]
prices = [0.9, 0.7, 0.5]
qs = ZCBPrice.(prices, times)

ufr = 0.03
α = 0.1

model = fit(Yield.SmithWilson(ufr=ufr, α=α), qs)
```

# Extended Help

## References

- [Smith-Wilson Yields Curves](http://gli.lu/2017/12/smith-wilson-yield-curves/)
- [A Technical Note on the Smith-Wilson Method](http://www.ressources-actuarielles.net/EXT/ISFA/fp-isfa.nsf/2b0481298458b3d1c1256f8a0024c478/bd689cce9bb2aeb5c1257998001ede2b/\$FILE/A_Technical_Note_on_the_Smith-Wilson_Method_100701.pdf)

"""
struct SmithWilson{TU <: AbstractVector, TQb <: AbstractVector, U, A} <: AbstractYieldModel
    u::TU
    qb::TQb
    ufr::U
    α::A

    # Inner constructor ensures that vector lengths match, and stores the ufr as a continuous number
    function SmithWilson(u::TU, qb::TQb, ufr, α::A) where {TU <: AbstractVector, TQb <: AbstractVector, A}
        if length(u) != length(qb)
            throw(DomainError("Vectors u and qb in SmithWilson must have equal length"))
        end
        z = __continuous(ufr)
        return new{TU, TQb, typeof(z), A}(u, qb, z, α)
    end
end


function SmithWilson(u, qb; ufr, α)
    return SmithWilson(u, qb, ufr, α)
end

# uninitialized rates used for `fit`
function SmithWilson(; ufr, α)
    return SmithWilson(Float64[], Float64[]; ufr, α)
end


function SmithWilson(times::AbstractVector, cashflows::AbstractMatrix, prices::AbstractVector; ufr, α)
    Q = Diagonal(exp.(-__continuous(ufr) * times)) * cashflows
    q = vec(sum(Q, dims = 1))  # We want q to be a column vector
    QHQ = Q' * H(α, times) * Q
    b = QHQ \ (prices - q)
    Qb = Q * b
    return SmithWilson(times, Qb; ufr = ufr, α = α)
end

# D(t) = exp(-ufr·t)·(1 + s(t)) with s(t) = Σ H(α, uᵢ, t)·qbᵢ (0 at t = 0), accumulated in a fused
# loop: the vector form `H(α, u, t) ⋅ qb` allocated a fresh length(u) array on every call. `H` is
# inlined into the loop (about a third faster). Fitted to arbitrary prices, 1 + s can be negative,
# and so can the discount factor.
function __smith_wilson_s(sw::SmithWilson, t)
    s = H(sw.α, sw.u[1], t) * sw.qb[1]
    @inbounds for i in 2:length(sw.u)
        s += H(sw.α, sw.u[i], t) * sw.qb[i]
    end
    return s
end

function FinanceCore.discount(sw::SmithWilson, t)
    isempty(sw.u) && return exp(-sw.ufr * t)
    return exp(-sw.ufr * t) * (1 + __smith_wilson_s(sw, t))
end

# The interval factor exp(-ufr·(to − from))·(1 + s(to))/(1 + s(from)) is exact for either sign of
# 1 + s and stays finite where both discount factors underflow. s(0) = 0, so from an exact 0 (every
# present value) it is `discount(sw, to)` and needs one sum; a time dual at 0 still takes the ratio.
function FinanceCore.discount(sw::SmithWilson, from, to)
    d = exp(-sw.ufr * (to - from))
    if !isempty(sw.u)
        g = 1 + __smith_wilson_s(sw, to)
        d *= iszero(from) ? g : g / (1 + __smith_wilson_s(sw, from))
    end
    return from == to ? one(d) : d
end

# L(t) = ufr·t − log1p(s(t)), for `zero` and `forward`; it is undefined where 1 + s ≤ 0 (a DomainError).
__log_discount(sw::SmithWilson, t) = isempty(sw.u) ? sw.ufr * t : sw.ufr * t - log1p(__smith_wilson_s(sw, t))

# L(to) − L(from) as the log of the interval factor above: real wherever that factor is positive,
# including where both discount factors are negative, and finite in the far tail. From an exact 0 it
# needs one sum, as the factor does.
function __log_interval(sw::SmithWilson, from, to)
    L = sw.ufr * (to - from)
    isempty(sw.u) && return L
    g = 1 + __smith_wilson_s(sw, to)
    return L - log(iszero(from) ? g : g / (1 + __smith_wilson_s(sw, from)))
end


"""
    H_ordered(α, t_min, t_max)

The Smith-Wilson H function with ordered arguments (for better performance than using min and max).
"""
@inline function H_ordered(α, t_min, t_max)
    return α * t_min + exp(-α * t_max) * sinh(-α * t_min)
end

"""
    H(α, t1, t2)

The Smith-Wilson H function implemented in a faster way.
"""
@inline function H(α, t1::T, t2::T) where {T}
    return t1 < t2 ? H_ordered(α, t1, t2) : H_ordered(α, t2, t1)
end

H(α, t1, t2) = H(α, promote(t1, t2)...)

H(α, t1vec::AbstractVector, t2) = [H(α, t1, t2) for t1 in t1vec]
H(α, t1vec::AbstractVector, t2vec::AbstractVector) = [H(α, t1, t2) for t1 in t1vec, t2 in t2vec]
# This can be optimized by going to H_ordered directly, but it might be a bit cumbersome
H(α, tvec::AbstractVector) = H(α, tvec, tvec)
