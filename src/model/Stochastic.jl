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
- `b`: long-term mean rate (continuous compounding). Can be passed as a `Real` or `Continuous(b)`.
- `σ`: volatility
- `initial`: initial short rate `r₀` (a `Rate` object or `Real`)

!!! note
    The Vasicek model allows negative rates. For very negative rates or long horizons,
    discount factors may exceed 1.
"""
struct Vasicek{A,B,S,T} <: AbstractStochasticModel
    a::A
    b::B
    σ::S
    initial::T
end

function Vasicek(a::Real, b::Real, σ::Real, initial::Real)
    return Vasicek(Float64(a), Float64(b), Float64(σ), Continuous(initial))
end

# Accept Continuous(b) for the long-term mean rate parameter
Vasicek(a, b::Rate{<:Any,Continuous}, σ, initial) = Vasicek(a, rate(b), σ, initial)

"""
    CoxIngersollRoss(a, b, σ, initial)

Cox-Ingersoll-Ross (1985) mean-reverting short-rate model:

    dr = a(b - r) dt + σ √r dW

# Arguments
- `a`: speed of mean reversion
- `b`: long-term mean rate (continuous compounding). Can be passed as a `Real` or `Continuous(b)`.
- `σ`: volatility
- `initial`: initial short rate `r₀` (a `Rate` object or `Real`)

!!! note "Feller condition"
    The condition `2ab > σ²` is required for the variance process to stay strictly
    positive. When violated, the short rate can reach zero; simulation uses absorption
    at zero (full truncation scheme) in that case.
"""
struct CoxIngersollRoss{A,B,S,T} <: AbstractStochasticModel
    a::A
    b::B
    σ::S
    initial::T
end

function CoxIngersollRoss(a::Real, b::Real, σ::Real, initial::Real)
    return CoxIngersollRoss(Float64(a), Float64(b), Float64(σ), Continuous(initial))
end

# Accept Continuous(b) for the long-term mean rate parameter
CoxIngersollRoss(a, b::Rate{<:Any,Continuous}, σ, initial) = CoxIngersollRoss(a, rate(b), σ, initial)

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
struct HullWhite{A,S,C} <: AbstractStochasticModel
    a::A
    σ::S
    curve::C
end

end # module ShortRate

# ─── Closed-form discount (zero-coupon bond prices) ──────────────────────────

function _initial_rate(m::ShortRate.Vasicek)
    return rate(Continuous(m.initial))
end

function _initial_rate(m::ShortRate.CoxIngersollRoss)
    return rate(Continuous(m.initial))
end

# Vasicek ZCB price: P = A(τ) exp(-B(τ) r)
# For small |a|, the general formula suffers from catastrophic cancellation
# (two O(σ²τ²/a) terms nearly cancel). Use a Taylor expansion instead.
function _vasicek_zcb(a, b, σ, r, τ)
    if abs(a * τ) < 0.02
        # Taylor expansion in a, avoiding cancellation:
        # B = τ - aτ²/2 + a²τ³/6 - a³τ⁴/24
        # lnA = σ²τ³/6 - a(bτ²/2 + σ²τ⁴/8) + a²(bτ³/6 + 7σ²τ⁵/120)
        B = τ * (1 - a * τ / 2 + (a * τ)^2 / 6 - (a * τ)^3 / 24)
        lnA = σ^2 * τ^3 / 6 -
               a * (b * τ^2 / 2 + σ^2 * τ^4 / 8) +
               a^2 * (b * τ^3 / 6 + 7 * σ^2 * τ^5 / 120)
    else
        B = (1 - exp(-a * τ)) / a
        lnA = (B - τ) * (a^2 * b - 0.5 * σ^2) / a^2 - σ^2 * B^2 / (4a)
    end
    return exp(lnA - B * r)
end

function FinanceCore.discount(m::ShortRate.Vasicek, T)
    return _vasicek_zcb(m.a, m.b, m.σ, _initial_rate(m), T)
end

# CIR ZCB price: P = A(τ) exp(-B(τ) r)
function _cir_zcb(a, b, σ, r, τ)
    if abs(σ) < 1e-15
        # Deterministic limit: dr = a(b-r)dt → r(t) = b + (r0-b)exp(-at)
        # P(0,τ) = exp(-∫₀ᵗ r(s)ds) = exp(-(bτ + (r-b)(1-exp(-aτ))/a))
        if abs(a) < 1e-12
            return exp(-r * τ)
        else
            return exp(-(b * τ + (r - b) * (1 - exp(-a * τ)) / a))
        end
    end
    γ = sqrt(a^2 + 2σ^2)
    expγτ = exp(γ * τ)
    denom = (γ + a) * (expγτ - 1) + 2γ
    B = 2(expγτ - 1) / denom
    A = (2γ * exp((a + γ) * τ / 2) / denom)^(2a * b / σ^2)
    return A * exp(-B * r)
end

function FinanceCore.discount(m::ShortRate.CoxIngersollRoss, T)
    return _cir_zcb(m.a, m.b, m.σ, _initial_rate(m), T)
end

# Hull-White is calibrated to match the initial term structure exactly.
# The model parameters (a, σ) affect derivative pricing and simulation,
# not the initial curve discount factors.
function FinanceCore.discount(m::ShortRate.HullWhite, T)
    return FinanceCore.discount(m.curve, T)
end

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
    P0t = FinanceCore.discount(m.curve, t)
    P0T = FinanceCore.discount(m.curve, T)
    B_tT = _hw_B(a, t, T)
    f0t = _hw_forward_rate(m.curve, t)
    if abs(a) < 1e-12
        lnA = log(P0T / P0t) + B_tT * f0t - 0.5 * σ^2 * t * B_tT^2
    else
        lnA = log(P0T / P0t) + B_tT * f0t - σ^2 / (4a) * B_tT^2 * (1 - exp(-2a * t))
    end
    return exp(lnA - B_tT * r_t)
end

# ─── RatePath: a simulated scenario as a yield model ─────────────────────────

"""
    RatePath(interp)

A simulated interest-rate path wrapped as an `AbstractYieldModel`.
`interp` maps time `t` to the cumulative integral ∫₀ᵗ r(s) ds so that
`discount(path, t) = exp(-interp(t))`.
"""
struct RatePath{I} <: Yield.AbstractYieldModel
    interp::I
end

function FinanceCore.discount(p::RatePath, t)
    return exp(-p.interp(t))
end

# ─── simulate: Euler-Maruyama path generation ────────────────────────────────

"""
    simulate(model::AbstractStochasticModel;
             n_scenarios=1000, timestep=1/12, horizon=30.0,
             rng=Random.default_rng())

Generate `n_scenarios` interest-rate paths via Euler-Maruyama discretisation.
Each path is returned as a `RatePath` (an `AbstractYieldModel`) so it plugs
directly into `present_value`, `discount`, etc.
"""
function simulate(model::AbstractStochasticModel;
                  n_scenarios::Int = 1000,
                  timestep::Real = 1 / 12,
                  horizon::Real = 30.0,
                  rng::Random.AbstractRNG = Random.default_rng())
    timestep > 0 || throw(ArgumentError("timestep must be positive, got $timestep"))
    dt = Float64(timestep)
    n_steps = ceil(Int, horizon / dt)
    sqrt_dt = sqrt(dt)

    drift_cache = _precompute_drift(model, dt, n_steps)

    # Determine element type from initial rate (may be a ForwardDiff.Dual)
    r0 = _sim_initial_rate(model)
    T = typeof(r0)

    paths = Vector{RatePath}(undef, n_scenarios)
    for i in 1:n_scenarios
        times = Vector{Float64}(undef, n_steps + 1)
        cumulative = Vector{T}(undef, n_steps + 1)
        times[1] = 0.0
        cumulative[1] = zero(T)
        r = r0
        for j in 1:n_steps
            Z = randn(rng)
            t = (j - 1) * dt
            r_new = _step(model, r, dt, sqrt_dt, Z, t, drift_cache, j)
            times[j + 1] = j * dt
            cumulative[j + 1] = cumulative[j] + 0.5 * (r + r_new) * dt
            r = r_new
        end
        interp = DataInterpolations.LinearInterpolation(
            cumulative, times;
            extrapolation = DataInterpolations.ExtrapolationType.Extension
        )
        paths[i] = RatePath(interp)
    end
    return paths
end

# Initial rate extractors for simulation
_sim_initial_rate(m::ShortRate.Vasicek) = _initial_rate(m)
_sim_initial_rate(m::ShortRate.CoxIngersollRoss) = _initial_rate(m)
function _sim_initial_rate(m::ShortRate.HullWhite)
    return _hw_forward_rate(m.curve, 0.0)
end

# Euler-Maruyama step for each model
function _step(m::ShortRate.Vasicek, r, dt, sqrt_dt, Z, t, ::Nothing, j)
    return r + m.a * (m.b - r) * dt + m.σ * sqrt_dt * Z
end

function _step(m::ShortRate.CoxIngersollRoss, r, dt, sqrt_dt, Z, t, ::Nothing, j)
    # Full truncation scheme (Lord, Koekkoek & Van Dijk, 2010)
    r_pos = max(r, 0.0)
    return max(r_pos + m.a * (m.b - r_pos) * dt + m.σ * sqrt(r_pos) * sqrt_dt * Z, 0.0)
end

function _step(m::ShortRate.HullWhite, r, dt, sqrt_dt, Z, t, θ_cache, j)
    return r + (θ_cache[j] - m.a * r) * dt + m.σ * sqrt_dt * Z
end

# Compute θ(t) = f_t(0,t) + a·f(0,t) + (σ²/2a)·(1 - exp(-2at))
# where f(0,t) is the instantaneous forward rate from the initial curve
function _hw_theta(m::ShortRate.HullWhite, t)
    a, σ = m.a, m.σ
    f0t = _hw_forward_rate(m.curve, t)
    df_dt = DifferentiationInterface.derivative(
        s -> _hw_forward_rate(m.curve, s),
        AutoForwardDiff(), max(t, 1e-10)
    )
    if abs(a) < 1e-12
        return df_dt + σ^2 * t
    else
        return df_dt + a * f0t + σ^2 / (2a) * (1 - exp(-2a * t))
    end
end

# Pre-compute drift values: nothing for Vasicek/CIR, Vector{Float64} for HW
_precompute_drift(::AbstractStochasticModel, dt, n_steps) = nothing
function _precompute_drift(m::ShortRate.HullWhite, dt, n_steps)
    return [_hw_theta(m, (j - 1) * dt) for j in 1:n_steps]
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

The `horizon` should cover the contract's maturity. The default (`maturity + 1`) ensures this.
"""
function pv_mc(model::AbstractStochasticModel, contract;
               n_scenarios::Int = 1000,
               timestep::Real = 1 / 12,
               horizon::Union{Nothing,Real} = nothing,
               rng::Random.AbstractRNG = Random.default_rng())
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

# Hull-White B(t,T) function
function _hw_B(a, t, T)
    τ = T - t
    if abs(a) < 1e-12
        return τ
    else
        return (1 - exp(-a * τ)) / a
    end
end

# Instantaneous forward rate f(0,t) = -d/dt ln P(0,t) via automatic differentiation
function _hw_forward_rate(curve, t)
    t_eval = max(t, 1e-10)  # avoid exactly zero for AD stability
    return -DifferentiationInterface.derivative(
        s -> log(FinanceCore.discount(curve, s)),
        AutoForwardDiff(), t_eval
    )
end

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

    # σ_P: volatility of the ZCB price at expiry
    B_TS = _hw_B(a, T, S)
    if abs(a) < 1e-12
        σ_P = σ * B_TS * sqrt(T)
    else
        σ_P = σ * B_TS * sqrt((1 - exp(-2a * T)) / (2a))
    end

    if σ_P < 1e-15
        # Degenerate case: no vol → intrinsic value
        call = max(P0S - K * P0T, 0.0)
        put  = max(K * P0T - P0S, 0.0)
        return (call, put)
    end

    h = (1 / σ_P) * log(P0S / (K * P0T)) + σ_P / 2

    call = P0S * N(h) - K * P0T * N(h - σ_P)
    put  = K * P0T * N(-h + σ_P) - P0S * N(-h)
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

# ─── present_value for Caps and Floors ───────────────────────────────────────
#
# A caplet paying max(L(T_{i-1},T_i) - K, 0)·τ at T_i is equivalent to
# (1 + K·τ) puts on a ZCB with maturity T_i, strike 1/(1+K·τ), expiring at T_{i-1}.
#
# Similarly a floorlet = (1 + K·τ) calls on a ZCB.

function FinanceCore.present_value(m::_GaussianModel, c::Option.Cap)
    K = c.strike
    freq = c.frequency isa FinanceCore.Frequency ? c.frequency.frequency : c.frequency
    τ = 1.0 / freq
    # Payment dates: τ, 2τ, ..., maturity
    # Caplet i: reset at T_{i-1}, pays at T_i
    # First caplet (reset at 0, pay at τ) is typically excluded (rate already known)
    n_periods = _check_integer_periods(c.maturity, freq, "Cap maturity")
    K_bond = 1.0 / (1.0 + K * τ)
    total = 0.0
    for i in 2:n_periods
        T_reset = (i - 1) * τ   # option expiry = reset date
        T_pay   = i * τ         # bond maturity = payment date
        _, put = _zcb_option_price(m, T_reset, T_pay, K_bond)
        total += (1.0 + K * τ) * put
    end
    return total
end

function FinanceCore.present_value(m::_GaussianModel, c::Option.Floor)
    K = c.strike
    freq = c.frequency isa FinanceCore.Frequency ? c.frequency.frequency : c.frequency
    τ = 1.0 / freq
    n_periods = _check_integer_periods(c.maturity, freq, "Floor maturity")
    K_bond = 1.0 / (1.0 + K * τ)
    total = 0.0
    for i in 2:n_periods
        T_reset = (i - 1) * τ
        T_pay   = i * τ
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

function FinanceCore.present_value(m::_GaussianModel, c::Option.Swaption)
    T0 = c.expiry
    freq = c.frequency isa FinanceCore.Frequency ? c.frequency.frequency : c.frequency
    τ = 1.0 / freq
    coupon = c.strike

    # Payment dates of the underlying swap
    n_payments = _check_integer_periods(c.swap_maturity - T0, freq, "Swap tenor")
    payment_times = [T0 + i * τ for i in 1:n_payments]

    # Step 1: Find r* such that the swap has zero value at T0
    # Swap value at T0 given r(T0) = r:
    #   V(r) = 1 - P(T0,Tn;r) - c·τ·∑ P(T0,Ti;r)
    function swap_value(r)
        total = 1.0
        for (i, Ti) in enumerate(payment_times)
            P_Ti = FinanceCore.discount(m, T0, Ti, r)
            total -= coupon * τ * P_Ti
            if i == length(payment_times)
                total -= P_Ti  # principal repayment
            end
        end
        return total
    end

    # Bisection to find r*
    r_star = _bisect(swap_value, -1.0, 1.0)

    # Step 2: Compute strike prices Ki = P(T0, Ti; r*)
    # Step 3: Sum ZCB options
    price = 0.0
    for (i, Ti) in enumerate(payment_times)
        Ki = FinanceCore.discount(m, T0, Ti, r_star)
        if c.payer
            # Payer swaption = sum of ZCB puts
            _, put = _zcb_option_price(m, T0, Ti, Ki)
            price += coupon * τ * put
            if i == length(payment_times)
                price += put  # principal
            end
        else
            # Receiver swaption = sum of ZCB calls
            call, _ = _zcb_option_price(m, T0, Ti, Ki)
            price += coupon * τ * call
            if i == length(payment_times)
                price += call  # principal
            end
        end
    end
    return price
end

# Validate that a value is an integer multiple of the period length
function _check_integer_periods(value, freq, label)
    n = value * freq
    n_int = round(Int, n)
    abs(n - n_int) < 1e-8 || throw(ArgumentError(
        "$label ($value) must be an integer multiple of the period length (1/$freq)"))
    return n_int
end

# Simple bisection solver with automatic bracket widening
function _bisect(f, lo, hi; tol = 1e-12, maxiter = 200)
    flo, fhi = f(lo), f(hi)
    if sign(flo) == sign(fhi)
        # Try widening the bounds
        for _ in 1:5
            lo, hi = 2lo, 2hi
            flo, fhi = f(lo), f(hi)
            sign(flo) != sign(fhi) && break
        end
        if sign(flo) == sign(fhi)
            error("Bisection: no sign change found in [$lo, $hi]; f(lo)=$flo, f(hi)=$fhi")
        end
    end
    for _ in 1:maxiter
        mid = (lo + hi) / 2
        fmid = f(mid)
        if abs(fmid) < tol || (hi - lo) / 2 < tol
            return mid
        end
        if sign(fmid) == sign(flo)
            lo = mid
            flo = fmid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

# ─── short_rate: extract r(t) from a simulated RatePath ──────────────────────

"""
    short_rate(path::RatePath, t)

The instantaneous short rate `r(t)` for a simulated scenario.

`RatePath` stores the cumulative integral `∫₀ᵗ r(s) ds` as a `LinearInterpolation`.
The short rate is the derivative of this cumulative integral.
"""
function short_rate(path::RatePath, t)
    return DataInterpolations.derivative(path.interp, t)
end
