# The long-run behavior of a cumulative log-discount: L(t) = a2·t² + a1·t + a0 + o(1) as t → ∞.
# A composite curve combines its components' tails before taking the limit, since the limits
# alone lose it: flat forwards of 4% and −2% have L = Inf and −Inf, but their sum tends to +Inf.
# A NaN coefficient is one the curve doesn't determine (see the fallback in Interface.jl). A known
# forward a1 is finite, so ±a1 = Inf marks growth faster than t of unknown order and that sign: it
# survives addition, scaling and rebasing (and Inf − Inf is NaN, unknown), but a known quadratic
# term of the other sign could outgrow it.
struct __LogTail{T}
    a2::T
    a1::T
    a0::T
end
__LogTail(a2, a1, a0) = __LogTail(promote(a2, a1, a0)...)

# Tails add, subtract and scale coefficientwise, as their log-discounts do.
Base.:+(x::__LogTail, y::__LogTail) = __LogTail(x.a2 + y.a2, x.a1 + y.a1, x.a0 + y.a0)
Base.:-(x::__LogTail, y::__LogTail) = __LogTail(x.a2 - y.a2, x.a1 - y.a1, x.a0 - y.a0)
Base.:*(k::Real, x::__LogTail) = __LogTail(k * x.a2, k * x.a1, k * x.a0)
# L(τ + t) − L(τ), given L(τ): the tail of a curve rebased to start at τ.
__shift_tail(x::__LogTail, τ, L_τ) = __LogTail(x.a2, x.a1 + 2 * x.a2 * τ, x.a0 + (x.a1 + x.a2 * τ) * τ - L_τ)

# `discount` at t = Inf: the limit of exp(-L(t)). The sign of a2, and then of a1, decides it: 0 when
# positive, Inf when negative. When both are zero it is finite, exp(-a0), which z(Inf)·Inf = 0·Inf
# would lose as NaN. The signs are those of the primal values. A deciding coefficient that is zero
# there but carries partials has no derivative: a bump that moves it can move the limit to 0 or Inf.
function __discount_at_infinity(tail::__LogTail)
    d = exp(-tail.a0)
    limit = __tail_limit_sign(tail)
    return limit > 0 ? zero(d) : limit < 0 ? oftype(d, Inf) : limit == 0 ? d : oftype(d, NaN)
end

# The same limit as a cumulative log-discount: +Inf where the discount factor tends to 0, -Inf
# where it tends to Inf, and a0 where it is finite. It is not formed as
# `-log(__discount_at_infinity(tail))`, so each result keeps its own numeric type and partials.
function __log_discount_at_infinity(tail::__LogTail)
    L = tail.a0
    limit = __tail_limit_sign(tail)
    return limit > 0 ? oftype(L, Inf) : limit < 0 ? oftype(L, -Inf) : limit == 0 ? L : oftype(L, NaN)
end

# +1 if the discount factor at infinity is 0, -1 if it is Inf, 0 if it is finite, NaN if the tail
# doesn't decide it. An unknown forward (NaN) or growth of unknown order (±Inf) is read first, since
# either could outgrow the quadratic term: the unknown growth decides the limit unless a2 has the
# other sign.
function __tail_limit_sign(tail::__LogTail)
    a1 = __primal(tail.a1)
    isnan(a1) && return NaN
    if isinf(a1)
        s = a1 > 0 ? 1.0 : -1.0
        a2 = __primal(tail.a2)
        a2 * s > 0 && return s
        a2 == 0 || return NaN
        # a2 = 0 decides too: a bump to the other sign could outgrow the unknown growth
        __check_limit_derivative(tail.a2, "slope")
        return s
    end
    for (c, name) in ((tail.a2, "slope"), (tail.a1, "forward"))
        v = __primal(c)
        v > 0 && return 1.0
        v < 0 && return -1.0
        isnan(v) && return NaN
        __check_limit_derivative(c, name)
    end
    return 0.0
end

# A deciding coefficient whose primal value is zero has no derivative of the limit unless its partials
# are zero too (ForwardDiff 1.x `iszero` checks both): a bump that moves it can move the limit.
__check_limit_derivative(c, name) = iszero(c) || throw(
    ArgumentError(
        "discount(curve, Inf) has no derivative here: the tail's $name is zero, so a bump " *
            "that moves it can send the discount factor at infinity to 0 or Inf. " *
            "Differentiate at a finite time."
    )
)

# A tail scaled by a blend weight, as a strong zero per coefficient: a zero weight gives a zero tail,
# also against a coefficient that is infinite or NaN.
__strong_zero_mul(k::Real, x::__LogTail) =
    __LogTail(__strong_zero_mul(k, x.a2), __strong_zero_mul(k, x.a1), __strong_zero_mul(k, x.a0))

# L(t) - L(0) of a curve with tail x: its tail from its interval factors.
__rebase_tail(x::__LogTail, L0) = __LogTail(x.a2, x.a1, x.a0 - L0)

__tail_type(tails) = mapreduce(x -> typeof(x.a0), promote_type, tails)
__nan_tail(tails) = (T = __tail_type(tails); __LogTail(T(NaN), T(NaN), T(NaN)))

# The tail of a discount-factor mixture Σ w̃ᵢ·exp(-Lᵢ), with shares w̃ᵢ summing to 1: the tail of the
# curve whose L grows slowest, among those with a positive share, ordered by (a2, a1). A curve whose
# (a2, a1) is larger is negligible, its intercept included; curves tied in (a2, a1) combine their
# intercepts as a0 = -log Σ w̃ᵢ·e^{-a0ᵢ}. The order is unknown, and the tail NaN, where a curve's a2
# or a1 is NaN, or where an unknown growth faster than t (a1 = ±Inf, see `__LogTail`) meets a
# different quadratic term that it could outgrow or be outgrown by.
#
# A zero share with partials on a curve that would dominate or tie, or tied coefficients with
# different partials, would change which curves decide the limit under a bump, so the limit has no
# derivative there.
function __mixture_tail(w̃, tails)
    T = promote_type(mapreduce(typeof, promote_type, w̃), __tail_type(tails))
    nan = __LogTail(T(NaN), T(NaN), T(NaN))
    pos(i) = __primal(w̃[i]) > 0
    key(i) = (__primal(tails[i].a2), __primal(tails[i].a1))
    any(i -> pos(i) && any(isnan, key(i)), eachindex(tails)) && return nan
    j = argmin(i -> pos(i) ? key(i) : (oftype(key(i)[1], Inf), oftype(key(i)[2], Inf)), eachindex(tails))
    a2j, a1j = key(j)
    for i in eachindex(tails)
        a2i, a1i = key(i)
        if pos(i)
            a2i > a2j && (a1j == Inf || a1i == -Inf) && return nan
        elseif !iszero(w̃[i]) && (any(isnan, key(i)) || !(key(j) < key(i)))
            __throw_mixture_limit_derivative()
        end
    end
    tied(i) = pos(i) && key(i) == key(j)
    for i in eachindex(tails)
        tied(i) && !(tails[i].a2 == tails[j].a2 && tails[i].a1 == tails[j].a1) && __throw_mixture_limit_derivative()
    end
    r = argmin(i -> tied(i) ? __primal(tails[i].a0) : oftype(__primal(tails[i].a0), Inf), eachindex(tails))
    a0r = T(tails[r].a0)
    a0 = a0r - log(sum(i -> tied(i) ? T(w̃[i]) * exp(a0r - T(tails[i].a0)) : zero(T), eachindex(tails)))
    return __LogTail(T(tails[j].a2), T(tails[j].a1), a0)
end
@noinline __throw_mixture_limit_derivative() = throw(
    ArgumentError(
        "discount(blend, Inf) has no derivative here: a bump of a weight or of a curve's tail changes " *
            "which curves decide the discount factor at infinity. Differentiate at a finite time."
    )
)
