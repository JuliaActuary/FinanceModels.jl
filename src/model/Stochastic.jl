"""
Stochastic short-rate models (Vasicek, Cox-Ingersoll-Ross, Hull-White) that
implement the `AbstractYieldModel` interface via closed-form zero-coupon bond
prices.  They also support Monte Carlo simulation via `simulate` and `pv_mc`.
"""

"""
    AbstractStochasticModel <: Yield.AbstractYieldModel

Abstract supertype for stochastic short-rate models.
"""
abstract type AbstractStochasticModel <: Yield.AbstractYieldModel end

module ShortRate

    import ..Yield
    import ..AbstractStochasticModel
    import ..FinanceCore
    using ..FinanceCore: Continuous, Rate, rate

    """
        Vasicek(a, b, σ, initial)

    Vasicek (1977) mean-reverting short-rate model:

        dr = a(b - r) dt + σ dW

    # Arguments
    - `a`: speed of mean reversion
    - `b`: long-term mean rate. A `Rate` is converted to its continuously compounded equivalent.
    - `σ`: volatility
    - `initial`: initial short rate `r₀` (a `Rate` object or `Real`)

    !!! note
        The Vasicek model allows negative rates. For very negative rates or long horizons,
        discount factors may exceed 1.
    """
    struct Vasicek{A, B, S, T} <: AbstractStochasticModel
        a::A
        b::B
        σ::S
        initial::T
        function Vasicek(a::A, b::B, σ::S, initial::T) where {A, B, S, T}
            σ >= 0 || throw(ArgumentError("volatility σ must be non-negative, got $σ"))
            return new{A, B, S, T}(a, b, σ, initial)
        end
    end

    function Vasicek(a::Real, b::Real, σ::Real, initial::Real)
        return Vasicek(a, b, σ, Continuous(initial))
    end

    # Store the long-term mean as a continuous scalar, regardless of the input
    # compounding convention.
    Vasicek(a, b::Rate, σ, initial) = Vasicek(a, rate(Continuous(b)), σ, initial)

    """
        CoxIngersollRoss(a, b, σ, initial)

    Cox-Ingersoll-Ross (1985) mean-reverting short-rate model:

        dr = a(b - r) dt + σ √r dW

    # Arguments
    - `a`: speed of mean reversion
    - `b`: long-term mean rate. A `Rate` is converted to its continuously compounded equivalent.
    - `σ`: volatility
    - `initial`: initial short rate `r₀` (a `Rate` object or `Real`)

    !!! note "Feller condition"
        The condition `2ab > σ²` is required for the variance process to stay strictly
        positive. When violated, the short rate can reach zero; simulation uses the
        full truncation scheme (Lord, Koekkoek & Van Dijk, 2010), which lets an
        auxiliary process go negative while the observed rate is floored at zero.
        Discretisation bias grows with `σ²/(2ab)` — use a finer `timestep` when the
        Feller condition is strongly violated.
    """
    struct CoxIngersollRoss{A, B, S, T} <: AbstractStochasticModel
        a::A
        b::B
        σ::S
        initial::T
        function CoxIngersollRoss(a::A, b::B, σ::S, initial::T) where {A, B, S, T}
            σ >= 0 || throw(ArgumentError("volatility σ must be non-negative, got $σ"))
            initial_value = rate(Continuous(initial))
            initial_value >= 0 || throw(ArgumentError("initial rate must be non-negative for CIR, got $initial"))
            return new{A, B, S, T}(a, b, σ, initial)
        end
    end

    function CoxIngersollRoss(a::Real, b::Real, σ::Real, initial::Real)
        return CoxIngersollRoss(a, b, σ, Continuous(initial))
    end

    # Store the long-term mean as a continuous scalar, regardless of the input
    # compounding convention.
    CoxIngersollRoss(a, b::Rate, σ, initial) =
        CoxIngersollRoss(a, rate(Continuous(b)), σ, initial)

    """
        HullWhite(a, σ, curve)

    Hull-White (1990) one-factor model:

        dr = (θ(t) - a r) dt + σ dW

    where `θ(t)` is calibrated to fit the initial term structure `curve`.

    # Arguments
    - `a`: speed of mean reversion
    - `σ`: volatility
    - `curve`: an existing yield model providing the initial term structure
    """
    struct HullWhite{A, S, C} <: AbstractStochasticModel
        a::A
        σ::S
        curve::C
        function HullWhite(a::A, σ::S, curve::C) where {A, S, C}
            σ >= 0 || throw(ArgumentError("volatility σ must be non-negative, got $σ"))
            return new{A, S, C}(a, σ, curve)
        end
    end

end # module ShortRate

# ─── Closed-form discount (zero-coupon bond prices) ──────────────────────────

function _initial_rate(m::ShortRate.Vasicek)
    return rate(Continuous(m.initial))
end

function _initial_rate(m::ShortRate.CoxIngersollRoss)
    return rate(Continuous(m.initial))
end

# The affine short-rate models' factor ∫₀^τ e^{-as} ds = (1 - e^{-aτ})/a: the bond sensitivity
# B(τ) = -∂log P/∂r of Vasicek and Hull-White, and (at 2a) their variance factors. It is τ at a = 0.
# The closed form is accurate for every other a, but its ForwardDiff derivative in `a` cancels as
# aτ → 0, so |aτ| < 0.2 uses the Taylor polynomial of φ₂(x) = (e^{-x} - 1 + x)/x² (the factor is
# τ(1 - xφ₂(x))), whose omitted terms are below Float64 rounding there. The series' coefficients are
# exact in every precision, so their values at x = 0 are too; their truncation is set for Float64,
# which bounds a wider type inside the band to about 4e-18 relative.
const _AFFINE_SERIES_X = 0.2
const _φ2_COEFFS = ntuple(k -> (-1)^(k - 1) // factorial(k + 1), 11)
_φ2(x) = __evalpoly_exact(x, _φ2_COEFFS)
@inline function _decay_integral(a, τ)
    x = a * τ
    return abs(x) < _AFFINE_SERIES_X ? τ * (1 - x * _φ2(x)) : -expm1(-x) / a
end

# Vasicek ZCB price P = A(τ) exp(-B(τ) r):
#     -log P(τ) = B·r + (τ - B)·(b - σ²/(2a²)) + σ²B²/(4a) = B·r + b·(τ - B) - σ²τ³/2·h(aτ),
# with h(x) = (x - 3/2 + 2e^{-x} - e^{-2x}/2)/x³ = Σₖ (-x)ᵏ (2ᵏ⁺² - 2)/(k + 3)!. The first form, with
# B = -m/a and τ - B = (x + m)/a for m = expm1(-x), keeps the terms linear in τ in one product, so
# -log P follows the long rate b - σ²/(2a²) to ±∞; its σ² terms cancel as x = aτ → 0, so |x| < 0.2
# uses the second form with the Taylor polynomials of φ₂ and h, exact to Float64 rounding there.
const _VASICEK_H_COEFFS = ntuple(k -> (-1)^(k - 1) * (2^(k + 1) - 2) // factorial(k + 2), 13)
function _vasicek_log_zcb(a, b, σ, r, τ)
    x = a * τ
    if abs(x) < _AFFINE_SERIES_X
        p = x * _φ2(x)
        return τ * (1 - p) * r + b * τ * p - σ^2 * τ^3 / 2 * __evalpoly_exact(x, _VASICEK_H_COEFFS)
    end
    m = expm1(-x)
    # The drift term vanishes with its coefficient. Skipping it then keeps τ = ∞ finite for a > 0
    # (it is 0·∞ otherwise); `iszero` of a dual number also requires zero partials.
    c = b - σ^2 / (2a^2)
    drift = iszero(c) ? zero(c) : (x + m) / a * c
    return -m / a * r + drift + σ^2 * m^2 / (4a^3)
end
_vasicek_zcb(a, b, σ, r, τ) = exp(-_vasicek_log_zcb(a, b, σ, r, τ))

function FinanceCore.discount(m::ShortRate.Vasicek, T)
    return _vasicek_zcb(m.a, m.b, m.σ, _initial_rate(m), T)
end
Yield.__log_discount(m::ShortRate.Vasicek, T) = _vasicek_log_zcb(m.a, m.b, m.σ, _initial_rate(m), T)
Yield.__log_native(::ShortRate.Vasicek) = true
Yield.instantaneous_forward(m::ShortRate.Vasicek, T) = _vasicek_forward(m.a, m.b, m.σ, _initial_rate(m), T)

# The instantaneous forward f(τ) = d(-log P)/dτ from the Riccati equations B′ = 1 - aB and
# (-log A)′ = abB - σ²B²/2: f = r·e^{-aτ} + b·(1 - e^{-aτ}) - σ²B²/2, with 1 - aB written as e^{-aτ}.
function _vasicek_forward(a, b, σ, r, τ)
    return r * exp(-a * τ) - b * expm1(-a * τ) - σ^2 / 2 * _decay_integral(a, τ)^2
end

# CIR ZCB price P = A(τ) exp(-B(τ) r), with γ = √(a² + 2σ²) (Cox, Ingersoll & Ross 1985):
#     B = 2(e^{γτ} - 1)/D,  log A = (2ab/σ²)·log(2γ e^{(a+γ)τ/2}/D),  D = (γ + a)(e^{γτ} - 1) + 2γ.
# That form overflows e^{γτ} at long maturities, and as σ → 0 it raises a number within σ² of 1 to the
# power 2ab/σ² (σ = 1e-10 gave 5.6e57). With s = γ + a and d = γ - a, so that s·d = 2σ² and s + d = 2γ,
# dividing D by e^{γτ} gives D′ = s + d·e^{-γτ}, and for either sign of a
#     B = 2(1 - e^{-γτ})/D′,   log A = 2ab·G.
# The larger of s and d is γ + |a|; the other is formed as 2σ²/(γ + |a|) = (2σ/(γ + |a|))·σ, which
# neither cancels nor underflows σ². G is written so that nothing cancels as σ → 0:
#     G = w·ℓ(σ²w) - τ/s,             w = 2(1 - e^{-γτ})/(s·D′)            for a ≥ 0,
#     G = τ/d - log1p(σ²v)/σ²,        v = (e^{γτ} - 1)/(d·γ)               for a < 0,
# where ℓ(u) = log1p(u)/u. Where σ²v overflows (a < 0), log1p(σ²v) comes from log(σ²v) instead, with
# G's two terms combined into one slope in τ, so that τ = ∞ gives ±∞ rather than ∞ - ∞.
# Both reach the deterministic price continuously as σ → 0, and σ = 0 needs no case of its own: the
# smaller factor is then 0 and they give the deterministic (volatility-free Vasicek) price. For small γτ,
# where the divisions by γ, s and σ² are near 0/0 (exactly so at a = σ = 0, leaving derivatives NaN),
# B and the log A term come from series in X = (γτ/2)² (`_cir_series`); outside that band only a = σ = 0 remains,
# at τ = ∞, where γτ is 0·∞ and the rate stays constant.
function _cir_pieces(a, σ, τ)
    γ = hypot(a, σ, σ)
    big = γ + abs(a)
    small = 2σ / big * σ
    s, d = a >= 0 ? (big, small) : (small, big)
    q = -expm1(-γ * τ)
    D = s + d * exp(-γ * τ)
    return (; γ, s, d, q, D, B = 2q / D)
end
function _cir_log_zcb(a, b, σ, r, τ)
    if _cir_short(a, σ, τ)
        B, g = _cir_series(a, σ, τ)
        return B * r + (a * τ) * (b * τ) * g
    end
    # a = σ = 0 at τ = ∞ (γτ = 0·∞): the rate stays constant
    iszero(a) && iszero(σ) && return r * τ
    (; γ, s, d, q, D, B) = _cir_pieces(a, σ, τ)
    # Without a drift towards b, A = 1, also at τ = ∞ where τ/s is infinite
    iszero(a * b) && return B * r
    return B * r - 2a * b * _cir_log_a_ratio(a, σ, τ, γ, s, d, q, D)
end
# G = log A/(2ab), as described above
function _cir_log_a_ratio(a, σ, τ, γ, s, d, q, D)
    if a >= 0
        w = 2q / (s * D)
        return w * _log1p_ratio(σ^2 * w) - τ / s
    end
    γτ = γ * τ
    v = expm1(γτ) / (d * γ)
    isfinite(σ^2 * v) && return τ / d - v * _log1p_ratio(σ^2 * v)
    # log1p(σ²v) = ℓc + γτ + log1p(e^{-(ℓc + γτ)}) with ℓc = log(σ²q/(dγ)), and 1/d - γ/σ² =
    # -(a² + σ² - aγ)/(dσ²), a sum of positive terms for a < 0; one numerator over σ², so that σ² = 0
    # (underflow) or τ = ∞ gives an infinity of the right sign
    ℓc = 2log(σ) - log(d * γ) + log(q)
    return -(τ * (a^2 + σ^2 - a * γ) / d + ℓc + log1p(exp(-(ℓc + γτ)))) / σ^2
end

# For small γτ: in the dimensionless X = (γτ/2)² = (aτ/2)² + (στ)²/2, C = cosh(γτ/2) and
# Ŝ = sinh(γτ/2)/(γτ/2) are entire, and B = τŜ/(C + (aτ/2)Ŝ). C + (aτ/2)Ŝ is e^{aτ/2} at σ = 0
# (X = X₀ = (aτ/2)²), so with Δ̂ its divided difference between X₀ and X, the log A term is
#     log A = 2ab·G = -(aτ)(bτ)·g,   g = Δ̂e^{-aτ/2}·ℓ((στ)²/2·Δ̂e^{-aτ/2}),
# which never divides by σ² or γ. Only dimensionless powers are formed, so a long time with small
# rates (or the reverse) can't overflow one factor while underflowing another. With (γτ)² < 0.16 the
# series' omitted terms are below Float64 rounding.
const _CIR_SERIES_X2 = 0.16
const _CIR_C = ntuple(k -> 1 // factorial(2k - 2), 9)   # 1/(2j)!, j = 0, …, 8
const _CIR_S = ntuple(k -> 1 // factorial(2k - 1), 9)   # 1/(2j + 1)!
_cir_short(a, σ, τ) = (a * τ)^2 + 2(σ * τ)^2 < _CIR_SERIES_X2
function _cir_series(a, σ, τ)
    α = a * τ / 2
    y = (σ * τ)^2 / 2
    X0 = α^2
    X = X0 + y
    C = __evalpoly_exact(X, _CIR_C)
    Ŝ = __evalpoly_exact(X, _CIR_S)
    T = typeof(float(__primal(X)))
    Δ, dd, X0k = zero(X), zero(X), one(X0)   # Δ̂, (Xᵏ - X₀ᵏ)/(X - X₀), X₀ᵏ⁻¹
    for k in 2:length(_CIR_C)
        dd = X * dd + X0k
        X0k *= X0
        Δ += (convert(T, _CIR_C[k]) + α * convert(T, _CIR_S[k])) * dd
    end
    e = exp(-α)
    return τ * Ŝ / (C + α * Ŝ), Δ * e * _log1p_ratio(y * Δ * e)
end
_cir_zcb(a, b, σ, r, τ) = exp(-_cir_log_zcb(a, b, σ, r, τ))

# The instantaneous forward f(τ) = d(-log P)/dτ = r·B′ + ab·B, from the Riccati equation
# (log A)′ = -abB, with B from the same pieces as `_cir_log_zcb`. B′ = (2γ e^{-γτ/2}/D′)², with the
# exponential inside the square so that a huge 2γ/D′ (a < 0, tiny σ) meets e^{-γτ} before it is squared;
# for small γτ, B′ = 1 - aB - σ²B²/2 from the Riccati equation itself, smooth through a = σ = 0. As for
# the price, σ = 0 is an ordinary case, and only a = σ = 0 at τ = ∞ is handled apart.
function _cir_forward(a, b, σ, r, τ)
    if _cir_short(a, σ, τ)
        B, _ = _cir_series(a, σ, τ)
        return r * (1 - a * B - σ^2 * B^2 / 2) + a * b * B
    end
    iszero(a) && iszero(σ) && return r + zero(τ)
    (; γ, D, B) = _cir_pieces(a, σ, τ)
    return r * (2γ * exp(-γ * τ / 2) / D)^2 + a * b * B
end

# log1p(u)/u. Where |u| < 1e-3 its Taylor polynomial, whose omitted terms are below Float64 rounding,
# so it is 1 at u = 0 (σ² underflows, or σ is a dual zero) and differentiable there.
const _LOG1P_RATIO_COEFFS = ntuple(k -> (-1)^(k - 1) // k, 7)
_log1p_ratio(u) = abs(u) < 1.0e-3 ? __evalpoly_exact(u, _LOG1P_RATIO_COEFFS) : log1p(u) / u

function FinanceCore.discount(m::ShortRate.CoxIngersollRoss, T)
    return _cir_zcb(m.a, m.b, m.σ, _initial_rate(m), T)
end
Yield.__log_discount(m::ShortRate.CoxIngersollRoss, T) = _cir_log_zcb(m.a, m.b, m.σ, _initial_rate(m), T)
Yield.__log_native(::ShortRate.CoxIngersollRoss) = true
Yield.instantaneous_forward(m::ShortRate.CoxIngersollRoss, T) = _cir_forward(m.a, m.b, m.σ, _initial_rate(m), T)

# Hull-White is calibrated to match the initial term structure exactly.
# The model parameters (a, σ) affect derivative pricing and simulation,
# not the initial curve discount factors.
function FinanceCore.discount(m::ShortRate.HullWhite, T)
    return FinanceCore.discount(m.curve, T)
end
Yield.__log_discount(m::ShortRate.HullWhite, T) = Yield.__log_discount(m.curve, T)
Yield.__log_interval(m::ShortRate.HullWhite, from, to) = Yield.__log_interval(m.curve, from, to)
Yield.__log_tail(m::ShortRate.HullWhite) = Yield.__log_tail(m.curve)
Yield.instantaneous_forward(m::ShortRate.HullWhite, T) = Yield.instantaneous_forward(m.curve, T)
Base.zero(m::ShortRate.HullWhite, T) = Base.zero(m.curve, T)
FinanceCore.discount(m::ShortRate.HullWhite, from, to) = FinanceCore.discount(m.curve, from, to)

# ─── Conditional discount P(t,T|r(t)) ────────────────────────────────────────

"""
    discount(m::ShortRate.Vasicek, t, T, r_t)

Conditional zero-coupon bond price ``P(t,T \\mid r(t) = r_t)`` under the Vasicek model.
Since the model is time-homogeneous, ``P(t,T|r) = P(0, T-t | r)``.
"""
FinanceCore.discount(m::ShortRate.Vasicek, t, T, r_t) = _vasicek_zcb(m.a, m.b, m.σ, r_t, T - t)

"""
    discount(m::ShortRate.CoxIngersollRoss, t, T, r_t)

Conditional zero-coupon bond price ``P(t,T \\mid r(t) = r_t)`` under the CIR model.
Since the model is time-homogeneous, ``P(t,T|r) = P(0, T-t | r)``.
"""
FinanceCore.discount(m::ShortRate.CoxIngersollRoss, t, T, r_t) = _cir_zcb(m.a, m.b, m.σ, r_t, T - t)

"""
    discount(m::ShortRate.HullWhite, t, T, r_t)

Conditional zero-coupon bond price ``P(t,T \\mid r(t) = r_t)`` under the Hull-White model.
Unlike Vasicek/CIR, this depends on `t` and `T` separately (not just `T-t`)
because the model is calibrated to an initial term structure.

Formula (Brigo & Mercurio 2006, Proposition 3.2.2):
```math
\\ln P(t,T) = \\ln\\frac{P(0,T)}{P(0,t)} + B(t,T) f(0,t) - \\frac{\\sigma^2}{4a} B(t,T)^2 (1 - e^{-2at}) - B(t,T) r_t
```
"""
function FinanceCore.discount(m::ShortRate.HullWhite, t, T, r_t)
    a, σ = m.a, m.σ
    B_tT = _decay_integral(a, T - t)
    f0t = Yield.instantaneous_forward(m.curve, t)
    # ln(P(0,T)/P(0,t)) is the curve's log-discount over [t, T], which stays finite where both factors
    # underflow; σ²/(4a)·(1 - e^{-2at}) = σ²/2 · ∫₀ᵗ e^{-2as} ds
    lnA = -Yield.__log_interval(m.curve, t, T) + B_tT * f0t - σ^2 / 2 * B_tT^2 * _decay_integral(2a, t)
    return exp(lnA - B_tT * r_t)
end

# ─── RatePath: a simulated scenario as a yield model ─────────────────────────

"""
    RatePath(interp::DataInterpolations.LinearInterpolation)

A simulated interest-rate path wrapped as an `AbstractYieldModel`. `interp` interpolates the
cumulative integral ∫₀ᵗ r(s) ds linearly over the path's time grid, so
`discount(path, t) = exp(-interp(t))` and the short rate is constant on each step. The grid starts
at `t = 0`, where the integral is 0; otherwise construction throws an `ArgumentError`. The path is
defined where `interp` is; `simulate` builds paths that throw outside their time grid.
"""
struct RatePath{I <: DataInterpolations.LinearInterpolation} <: Yield.AbstractYieldModel
    interp::I
    function RatePath(interp::I) where {I <: DataInterpolations.LinearInterpolation}
        t0, L0 = first(interp.t), first(interp.u)
        (iszero(t0) && iszero(L0)) ||
            throw(ArgumentError("a RatePath's grid must start at t = 0 with value 0, got $L0 at t = $t0"))
        return new{I}(interp)
    end
end

function FinanceCore.discount(p::RatePath, t)
    return exp(-p.interp(t))
end
Yield.__log_discount(p::RatePath, t) = p.interp(t)
Yield.__log_native(::RatePath) = true

# ─── simulate: path generation ───────────────────────────────────────────────

"""
    simulation_steps(horizon, timestep) -> (; nsteps, aligned)

The number of `timestep`-long steps a simulation to `horizon` takes, and whether `horizon` lies on
that grid. With `r = horizon / timestep`, `horizon` is aligned when `r` is within `8eps(r)` of an
integer `n`, and the grid has `nsteps = n` steps: roundoff (`0.07 / 0.01` is `7.000000000000001`)
adds no step. Otherwise `nsteps = ceil(r)`, the first grid point beyond `horizon`, and
`aligned = false`. `horizon` and `timestep` must be finite and positive.

[`simulate`](@ref) ends an aligned grid at `horizon` itself, and covers an unaligned horizon with
the extra step.
"""
function simulation_steps(horizon::Real, timestep::Real)
    (isfinite(horizon) && horizon > 0) || throw(ArgumentError("horizon must be finite and positive, got $horizon"))
    (isfinite(timestep) && timestep > 0) || throw(ArgumentError("timestep must be finite and positive, got $timestep"))
    r = float(horizon) / float(timestep)
    n = round(r)
    abs(r - n) <= 8eps(r) && return (; nsteps = Int(n), aligned = true)
    return (; nsteps = ceil(Int, r), aligned = false)
end

"""
    simulate(model::AbstractStochasticModel;
             n_scenarios=1000, timestep=1/12, horizon=30.0,
             rng=Random.default_rng())

Generate `n_scenarios` interest-rate paths.
Each path is returned as a `RatePath` (an `AbstractYieldModel`) so it plugs
directly into `present_value`, `discount`, etc. A path is defined on its time
grid of [`simulation_steps`](@ref)`(horizon, timestep).nsteps` steps, from 0 to
`horizon` (or the first grid point beyond an unaligned `horizon`); evaluating it
outside that range throws rather than extending the path.

Discretisation schemes:
- **Vasicek / Hull-White**: the exact Gaussian transition density, so the
  simulated short rate has no discretisation bias at any `timestep`.
- **Cox-Ingersoll-Ross**: full truncation (Lord, Koekkoek & Van Dijk, 2010) —
  an auxiliary process may go negative while the observed rate is `max(x, 0)`.
  Weak bias vanishes as `timestep → 0`, but is material for coarse steps when
  the Feller condition is strongly violated.

In all cases the cumulative discount integral ``∫₀ᵗ r(s)\\,ds`` is accumulated
with the trapezoidal rule between grid points, so pathwise discount factors
(and hence `pv_mc`) retain an integration error that grows with `timestep`
even when the short-rate transition itself is exact.
"""
function simulate(
        model::AbstractStochasticModel;
        n_scenarios::Int = 1000,
        timestep::Real = 1 / 12,
        horizon::Real = 30.0,
        rng::Random.AbstractRNG = Random.default_rng()
    )
    grid = simulation_steps(horizon, timestep)
    n_steps = grid.nsteps
    dt = Float64(timestep)
    sqrt_dt = sqrt(dt)

    cache = _sim_cache(model, dt, n_steps)

    # Parameters other than the initial rate may carry an AD number type.
    r0 = _sim_initial_rate(model)
    T = _sim_eltype(model, r0)

    times = Vector{Float64}(undef, n_steps + 1)
    times[1] = 0.0
    for j in 1:n_steps
        times[j + 1] = j * dt
    end
    # An aligned grid ends at the horizon itself: n_steps·dt can round to either side of it (three
    # steps of 0.3 end at 0.8999999999999999, three of 0.1 at 0.30000000000000004), and the path
    # covers exactly the horizon it was asked for. An unaligned grid ends at its extra step.
    grid.aligned && (times[end] = horizon)

    # `map` (rather than filling a Vector{RatePath}, a UnionAll eltype) infers the
    # concrete RatePath{...} element type, so downstream pricing loops dispatch
    # statically instead of dynamically per path
    return map(1:n_scenarios) do _
        cumulative = Vector{T}(undef, n_steps + 1)
        cumulative[1] = zero(T)
        r = r0
        for j in 1:n_steps
            Z = randn(rng)
            t = (j - 1) * dt
            r_new = _step(model, r, dt, sqrt_dt, Z, t, cache, j)
            cumulative[j + 1] = cumulative[j] +
                0.5 * (_observed_rate(model, r) + _observed_rate(model, r_new)) * dt
            r = r_new
        end
        # no extrapolation: the path is defined only on its simulated grid
        interp = DataInterpolations.LinearInterpolation(
            cumulative, times;
            extrapolation = DataInterpolations.ExtrapolationType.None
        )
        RatePath(interp)
    end
end

# Initial rate extractors for simulation
_sim_initial_rate(m::ShortRate.Vasicek) = _initial_rate(m)
_sim_initial_rate(m::ShortRate.CoxIngersollRoss) = _initial_rate(m)
function _sim_initial_rate(m::ShortRate.HullWhite)
    return Yield.instantaneous_forward(m.curve, 0.0)
end

# Include every model parameter that can flow into a simulated state. Falling
# back to the initial-rate type preserves the extension contract for custom
# stochastic models.
_sim_eltype(::AbstractStochasticModel, r0) = typeof(r0)
_sim_eltype(m::ShortRate.Vasicek, r0) =
    promote_type(typeof(r0), typeof(m.a), typeof(m.b), typeof(m.σ))
_sim_eltype(m::ShortRate.CoxIngersollRoss, r0) =
    promote_type(typeof(r0), typeof(m.a), typeof(m.b), typeof(m.σ))
_sim_eltype(m::ShortRate.HullWhite, r0) =
    promote_type(typeof(r0), typeof(m.a), typeof(m.σ))

# Exact Ornstein-Uhlenbeck transition parameters over one step of length dt:
#   x_{t+dt} | x_t ~ Normal(x_t·ϕ, sd²),  ϕ = exp(-a·dt),  sd² = σ²(1-ϕ²)/(2a)
function _ou_step_params(a, σ, dt)
    ϕ = exp(-a * dt)
    var = σ^2 * _decay_integral(2a, dt)
    return (ϕ = ϕ, sd = sqrt(var))
end

# The rate that enters the discount integral: identity except for CIR, where
# the state is the full-truncation auxiliary process and the rate is its
# positive part. The abstract fallback preserves the extension contract for
# custom stochastic models whose simulated state is their observed rate.
_observed_rate(::AbstractStochasticModel, r) = r
_observed_rate(::ShortRate.CoxIngersollRoss, x) = max(x, zero(x))

# Vasicek: exact transition r' = b + (r-b)ϕ + sd·Z (no discretisation bias)
function _step(m::ShortRate.Vasicek, r, dt, sqrt_dt, Z, t, cache, j)
    return m.b + (r - m.b) * cache.ϕ + cache.sd * Z
end

# CIR: full truncation scheme (Lord, Koekkoek & Van Dijk, 2010).
# The state is an auxiliary process x that may go negative; drift and
# diffusion use x⁺ and the observed rate is x⁺ (see `_observed_rate`).
# Flooring x itself at zero (absorption) would bias E[exp(-∫r)] down
# materially when the Feller condition is violated.
function _step(m::ShortRate.CoxIngersollRoss, x, dt, sqrt_dt, Z, t, ::Nothing, j)
    xp = max(x, zero(x))
    return x + m.a * (m.b - xp) * dt + m.σ * sqrt(xp) * sqrt_dt * Z
end

# Hull-White: exact transition via the decomposition r(t) = x(t) + α(t)
# where dx = -a·x dt + σ dW (OU with zero mean, simulated exactly) and
# α(t) = f(0,t) + σ²/(2a²)(1 - exp(-at))² (Brigo & Mercurio 2006, Eq. 3.36).
# This avoids differentiating the forward curve (θ(t) needs f_t(0,t)) and has
# no discretisation bias in r.
function _step(m::ShortRate.HullWhite, r, dt, sqrt_dt, Z, t, cache, j)
    x = r - cache.α[j]  # cache.α[j] = α(t_{j-1})
    return x * cache.ou.ϕ + cache.ou.sd * Z + cache.α[j + 1]
end

function _hw_alpha(m::ShortRate.HullWhite, t)
    a, σ = m.a, m.σ
    f0t = Yield.instantaneous_forward(m.curve, t)
    return f0t + σ^2 / 2 * _decay_integral(a, t)^2
end

# Per-simulation cache: OU transition parameters for the Gaussian models
# (plus the α(t) grid for Hull-White). Models that do not need a cache,
# including CIR and custom AbstractStochasticModel subtypes, receive `nothing`.
_sim_cache(::AbstractStochasticModel, dt, n_steps) = nothing
_sim_cache(m::ShortRate.Vasicek, dt, n_steps) = _ou_step_params(m.a, m.σ, dt)
function _sim_cache(m::ShortRate.HullWhite, dt, n_steps)
    return (
        ou = _ou_step_params(m.a, m.σ, dt),
        α = [_hw_alpha(m, j * dt) for j in 0:n_steps],
    )
end

# ─── pv_mc: Monte Carlo expected present value ───────────────────────────────

"""
    pv_mc(model, contract;
          n_scenarios=1000, timestep=1/12, horizon=nothing,
          rng=Random.default_rng())

Estimate the expected present value of `contract` under the stochastic `model`
by averaging `present_value` across simulated scenarios.

!!! note
    `pv_mc` is designed for fixed-cashflow instruments where each `RatePath` scenario
    provides the discount factors. For floating-rate instruments whose cashflows depend
    on the rate path, project cashflows per scenario using `Projection` instead.

The `horizon` must cover the contract's maturity, since a simulated path throws
beyond its horizon. The default (`maturity + 1`) ensures this.
"""
function pv_mc(
        model::AbstractStochasticModel, contract;
        n_scenarios::Int = 1000,
        timestep::Real = 1 / 12,
        horizon::Union{Nothing, Real} = nothing,
        rng::Random.AbstractRNG = Random.default_rng()
    )
    h = horizon === nothing ? Float64(FinanceModels.maturity(contract)) + 1.0 : Float64(horizon)
    scenarios = simulate(model; n_scenarios, timestep, horizon = h, rng)
    total = sum(FinanceCore.present_value(sc, contract) for sc in scenarios)
    return total / n_scenarios
end

# ─── Hull-White derivative pricing (ZCB options, caps, swaptions) ────────────
#
# Reference: Brigo & Mercurio (2006) "Interest Rate Models", Chapter 3
#            Hull (2018) "Options, Futures, and Other Derivatives", Ch. 32
#            Jamshidian (1989) "An Exact Bond Option Formula"

# Union type for Gaussian (normal) short-rate models that share the same
# ZCB option formula (Black's formula with σ_P from the B(t,T) function).
const _GaussianModel = Union{ShortRate.Vasicek, ShortRate.HullWhite}

"""
    _zcb_option_price(m::Union{ShortRate.Vasicek, ShortRate.HullWhite}, T, S, K)

Closed-form price of a European call and put on a zero-coupon bond
under a Gaussian (Vasicek or Hull-White) one-factor model.

Returns `(call_price, put_price)`.

- `T`: option expiry
- `S`: bond maturity (S > T)
- `K`: strike price

Reference: Brigo & Mercurio (2006), Proposition 3.2.1
"""
function _zcb_option_price(m::_GaussianModel, T, S, K)
    S > T || throw(ArgumentError("Bond maturity S=$S must be greater than option expiry T=$T"))
    K > 0 || throw(ArgumentError("Strike K=$K must be positive"))
    a, σ = m.a, m.σ
    P0T = FinanceCore.discount(m, T)
    P0S = FinanceCore.discount(m, S)

    # σ_P: volatility of the ZCB price at expiry, σ·B(T,S)·sqrt((1 - e^{-2aT})/(2a))
    σ_P = σ * _decay_integral(a, S - T) * sqrt(_decay_integral(2a, T))

    if σ_P < 1.0e-15
        # Degenerate case: no vol → intrinsic value
        call = max(P0S - K * P0T, 0.0)
        put = max(K * P0T - P0S, 0.0)
        return (call, put)
    end

    h = (1 / σ_P) * log(P0S / (K * P0T)) + σ_P / 2

    call = P0S * N(h) - K * P0T * N(h - σ_P)
    put = K * P0T * N(-h + σ_P) - P0S * N(-h)
    return (call, put)
end

# ─── present_value for ZCB options ───────────────────────────────────────────

function FinanceCore.present_value(m::_GaussianModel, c::Option.ZCBCall)
    call, _ = _zcb_option_price(m, c.expiry, c.bond_maturity, c.strike)
    return call
end

function FinanceCore.present_value(m::_GaussianModel, c::Option.ZCBPut)
    _, put = _zcb_option_price(m, c.expiry, c.bond_maturity, c.strike)
    return put
end

_frequency_value(f::FinanceCore.Frequency) = f.frequency
_frequency_value(f::Real) = f

# ─── present_value for Caps and Floors ───────────────────────────────────────
#
# A caplet paying max(L(T_{i-1},T_i) - K, 0)·τ at T_i is equivalent to
# (1 + K·τ) puts on a ZCB with maturity T_i, strike 1/(1+K·τ), expiring at T_{i-1}.
#
# Similarly a floorlet = (1 + K·τ) calls on a ZCB.

function FinanceCore.present_value(m::_GaussianModel, c::Option.Cap)
    K = c.strike
    freq = _frequency_value(c.frequency)
    τ = 1.0 / freq
    # Payment dates: τ, 2τ, ..., maturity
    # Caplet i: reset at T_{i-1}, pays at T_i
    # First caplet (reset at 0, pay at τ) is excluded: its rate is already known
    # at valuation (standard market convention; see Hull 2018, §32.3).
    # For forward-starting caps, adjust the contract maturity accordingly.
    n_periods = _check_integer_periods(c.maturity, freq, "Cap maturity")
    K_bond = 1.0 / (1.0 + K * τ)
    # Without a caplet (a maturity of one period or less) the cap is worth 0. As FinanceCore values
    # an empty collection, that 0 is `zero` of a caplet amount discounted at time zero: the type of a
    # present value under `m`, with no dependence on its value.
    n_periods < 2 && return zero((1.0 + K * τ) * FinanceCore.discount(m, zero(τ)))
    total = 0.0
    for i in 2:n_periods
        T_reset = (i - 1) * τ   # option expiry = reset date
        T_pay = i * τ         # bond maturity = payment date
        _, put = _zcb_option_price(m, T_reset, T_pay, K_bond)
        total += (1.0 + K * τ) * put
    end
    return total
end

function FinanceCore.present_value(m::_GaussianModel, c::Option.Floor)
    K = c.strike
    freq = _frequency_value(c.frequency)
    τ = 1.0 / freq
    n_periods = _check_integer_periods(c.maturity, freq, "Floor maturity")
    K_bond = 1.0 / (1.0 + K * τ)
    # without a floorlet, 0 in the type of a present value under `m` (as for a cap)
    n_periods < 2 && return zero((1.0 + K * τ) * FinanceCore.discount(m, zero(τ)))
    total = 0.0
    for i in 2:n_periods
        T_reset = (i - 1) * τ
        T_pay = i * τ
        call, _ = _zcb_option_price(m, T_reset, T_pay, K_bond)
        total += (1.0 + K * τ) * call
    end
    return total
end

# ─── present_value for European Swaptions (Jamshidian decomposition) ─────────
#
# A payer swaption = right to enter a pay-fixed swap at expiry T₀.
# The underlying swap has payment dates T₁,...,Tₙ with coupon c and frequency f.
# At expiry the swap value (per unit notional) is:
#   V(r) = 1 - P(T₀,Tₙ;r) - c·τ·∑ P(T₀,Tᵢ;r)
# which is positive when rates are high (payer benefits).
#
# Jamshidian (1989): find r* where V(r*)=0, then
#   Payer = ∑ [c·τ·Put(T₀,Tᵢ,Kᵢ)] + Put(T₀,Tₙ,Kₙ)
#   where Kᵢ = P(T₀,Tᵢ;r*)
#
# For a receiver swaption, replace Put with Call.
#
# NOTE: Jamshidian decomposition requires monotonic bond prices in r, which holds
# only for Gaussian models (Vasicek, Hull-White). For CIR, use pv_mc() instead.

function FinanceCore.present_value(m::_GaussianModel, c::Option.Swaption)
    T0 = c.expiry
    freq = _frequency_value(c.frequency)
    τ = 1.0 / freq
    coupon = c.strike

    # Payment dates of the underlying swap
    n_payments = _check_integer_periods(c.swap_maturity - T0, freq, "Swap tenor")
    payment_times = [T0 + i * τ for i in 1:n_payments]

    # Step 1: Find r* such that the swap has zero value at T0
    # Swap value at T0 given r(T0) = r:
    #   V(r) = 1 - P(T0,Tn;r) - c·τ·∑ P(T0,Ti;r)
    # Each bond price is affine in r, P(T0,Ti;r) = Aᵢ e^{-Bᵢ r}, so V is increasing in r
    # with the closed-form slope ∑ wᵢ Bᵢ P(T0,Ti;r) (wᵢ the coupon and principal weights).
    weight(i) = coupon * τ + (i == n_payments ? 1 : 0)
    function swap_value(r)
        total = 1.0
        for (i, Ti) in enumerate(payment_times)
            total -= weight(i) * FinanceCore.discount(m, T0, Ti, r)
        end
        return total
    end
    slope_terms(r) = (
        weight(i) * __primal(_decay_integral(m.a, Ti - T0)) * __primal(FinanceCore.discount(m, T0, Ti, r))
            for (i, Ti) in enumerate(payment_times)
    )

    # r* depends on the model's parameters. The root is solved on primal values and its
    # derivatives come from the implicit function theorem (a Float64 root would silently
    # drop the ∂Kᵢ/∂r*·dr*/dθ terms of every Greek).
    # The slope's own terms measure cancellation (negative strikes give negative coupon weights).
    r_star = __implicit_root(
        swap_value, r -> __primal(swap_value(r)), 0.0;
        who = "Swaption critical rate r*", slope = r -> sum(slope_terms(r)), scale = r -> sum(abs, slope_terms(r))
    )

    # Step 2: Compute strike prices Ki = P(T0, Ti; r*)
    # Step 3: Sum ZCB options, weighted like the swap's payments: a payer swaption is a
    # portfolio of ZCB puts, a receiver swaption of ZCB calls
    price = 0.0
    for (i, Ti) in enumerate(payment_times)
        call, put = _zcb_option_price(m, T0, Ti, FinanceCore.discount(m, T0, Ti, r_star))
        price += weight(i) * (c.payer ? put : call)
    end
    return price
end

# Validate that a value is an integer multiple of the period length
function _check_integer_periods(value, freq, label)
    n = value * freq
    n_int = round(Int, n)
    abs(n - n_int) < 1.0e-8 || throw(
        ArgumentError(
            "$label ($value) must be an integer multiple of the period length (1/$freq)"
        )
    )
    return n_int
end

# ─── short_rate: extract r(t) from a simulated RatePath ──────────────────────

"""
    short_rate(path::RatePath, t)

The instantaneous short rate `r(t)` for a simulated scenario.

`RatePath` stores the cumulative integral `∫₀ᵗ r(s) ds` as a `LinearInterpolation`.
The short rate is the derivative of this cumulative integral.

Because the cumulative integral is built from trapezoidal steps, the
returned rate is piecewise-constant within each timestep — an approximation to the
continuous short-rate process, not the exact value. It is right-continuous: at a grid time
it is the slope of the step that starts there, and at the path's last time the last step's.
It is the path's `Yield.instantaneous_forward`.
"""
short_rate(path::RatePath, t) = Yield.instantaneous_forward(path, t)

# The slope of L = ∫₀ᵗ r over the grid step that starts at `t` (the last step at the path's last time),
# so the rate is right-continuous like a knot curve's forward. Outside its grid the interpolant's own
# extrapolation decides: a simulated path throws, a path built with an extension extends.
function Yield.instantaneous_forward(p::RatePath, t)
    ts, L = p.interp.t, p.interp.u
    first(ts) <= t <= last(ts) || return DataInterpolations.derivative(p.interp, t)
    i = min(searchsortedlast(ts, t), length(ts) - 1)
    return (L[i + 1] - L[i]) / (ts[i + 1] - ts[i])
end
