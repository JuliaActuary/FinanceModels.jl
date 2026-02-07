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

# Vasicek ZCB price: P(0,T) = A(T) exp(-B(T) r₀)
# For small |a|, the general formula suffers from catastrophic cancellation
# (two O(σ²T²/a) terms nearly cancel). Use a Taylor expansion instead.
function FinanceCore.discount(m::ShortRate.Vasicek, T)
    a, b, σ, r0 = m.a, m.b, m.σ, _initial_rate(m)
    if abs(a * T) < 0.02
        # Taylor expansion in a, avoiding cancellation:
        # B = T - aT²/2 + a²T³/6 - a³T⁴/24
        # lnA = σ²T³/6 - a(bT²/2 + σ²T⁴/8) + a²(bT³/6 + 7σ²T⁵/120)
        B = T * (1 - a * T / 2 + (a * T)^2 / 6 - (a * T)^3 / 24)
        lnA = σ^2 * T^3 / 6 -
               a * (b * T^2 / 2 + σ^2 * T^4 / 8) +
               a^2 * (b * T^3 / 6 + 7 * σ^2 * T^5 / 120)
    else
        B = (1 - exp(-a * T)) / a
        lnA = (B - T) * (a^2 * b - 0.5 * σ^2) / a^2 - σ^2 * B^2 / (4a)
    end
    return exp(lnA - B * r0)
end

# CIR ZCB price: P(0,T) = A(T) exp(-B(T) r₀)
function FinanceCore.discount(m::ShortRate.CoxIngersollRoss, T)
    a, b, σ, r0 = m.a, m.b, m.σ, _initial_rate(m)
    if abs(σ) < 1e-15
        # Deterministic limit: dr = a(b-r)dt → r(t) = b + (r0-b)exp(-at)
        # P(0,T) = exp(-∫₀ᵀ r(s)ds) = exp(-(bT + (r0-b)(1-exp(-aT))/a))
        if abs(a) < 1e-12
            return exp(-r0 * T)
        else
            return exp(-(b * T + (r0 - b) * (1 - exp(-a * T)) / a))
        end
    end
    γ = sqrt(a^2 + 2σ^2)
    expγT = exp(γ * T)
    denom = (γ + a) * (expγT - 1) + 2γ
    B = 2(expγT - 1) / denom
    A = (2γ * exp((a + γ) * T / 2) / denom)^(2a * b / σ^2)
    return A * exp(-B * r0)
end

# Hull-White ZCB price: uses the initial curve + volatility correction
function FinanceCore.discount(m::ShortRate.HullWhite, T)
    a, σ = m.a, m.σ
    P0T = FinanceCore.discount(m.curve, T)
    if abs(a) < 1e-12
        B = T
    else
        B = (1 - exp(-a * T)) / a
    end
    # Hull-White discount is the curve discount (no volatility correction needed
    # for the initial yield curve – the model is calibrated to match it exactly)
    return P0T
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
    dt = Float64(timestep)
    n_steps = ceil(Int, horizon / dt)
    sqrt_dt = sqrt(dt)

    paths = Vector{RatePath}(undef, n_scenarios)
    for i in 1:n_scenarios
        times = Vector{Float64}(undef, n_steps + 1)
        cumulative = Vector{Float64}(undef, n_steps + 1)
        times[1] = 0.0
        cumulative[1] = 0.0
        r = _sim_initial_rate(model)
        for j in 1:n_steps
            Z = randn(rng)
            t = (j - 1) * dt
            r_new = _step(model, r, dt, sqrt_dt, Z, t)
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
    # Instantaneous short rate from curve: -d/dT ln P(0,T) at T→0
    ε = 1e-6
    return -log(FinanceCore.discount(m.curve, ε)) / ε
end

# Euler-Maruyama step for each model
function _step(m::ShortRate.Vasicek, r, dt, sqrt_dt, Z, t)
    return r + m.a * (m.b - r) * dt + m.σ * sqrt_dt * Z
end

function _step(m::ShortRate.CoxIngersollRoss, r, dt, sqrt_dt, Z, t)
    # Full truncation scheme (Lord, Koekkoek & Van Dijk, 2010)
    r_pos = max(r, 0.0)
    return max(r_pos + m.a * (m.b - r_pos) * dt + m.σ * sqrt(r_pos) * sqrt_dt * Z, 0.0)
end

function _step(m::ShortRate.HullWhite, r, dt, sqrt_dt, Z, t)
    # θ(t) = f_t(0,t) + a·f(0,t) + (σ²/2a)·(1 - exp(-2at))
    # where f(0,t) is the instantaneous forward rate from the initial curve
    a, σ = m.a, m.σ
    f0t = _hw_forward_rate(m.curve, t)
    # Use fixed ε for df/dt (not dt, which would make drift timestep-dependent)
    ε_deriv = 1e-4
    f0t_plus = _hw_forward_rate(m.curve, t + ε_deriv)
    df_dt = (f0t_plus - f0t) / ε_deriv
    if abs(a) < 1e-12
        θ = df_dt + σ^2 * t
    else
        θ = df_dt + a * f0t + σ^2 / (2a) * (1 - exp(-2a * t))
    end
    return r + (θ - a * r) * dt + σ * sqrt_dt * Z
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

# Instantaneous forward rate f(0,t) = -d/dt ln P(0,t) from the initial curve
function _hw_forward_rate(curve, t)
    ε = 1e-6
    if t < ε
        # One-sided difference at t ≈ 0
        return -log(FinanceCore.discount(curve, ε)) / ε
    else
        # Central difference
        return -(log(FinanceCore.discount(curve, t + ε)) - log(FinanceCore.discount(curve, t - ε))) / (2ε)
    end
end

# Hull-White ZCB price P(t,T) given r(t), used in Jamshidian decomposition
# ln P(t,T) = ln(P(0,T)/P(0,t)) + B(t,T)·f(0,t) - σ²/(4a)·B(t,T)²·(1-exp(-2at))
function _hw_zcb(m::ShortRate.HullWhite, t, T, r_t)
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

"""
    _hw_zcb_option_price(m::ShortRate.HullWhite, T, S, K)

Closed-form price of a European call and put on a zero-coupon bond
under the Hull-White one-factor model.

Returns `(call_price, put_price)`.

- `T`: option expiry
- `S`: bond maturity (S > T)
- `K`: strike price

Reference: Brigo & Mercurio (2006), Proposition 3.2.1
"""
function _hw_zcb_option_price(m::ShortRate.HullWhite, T, S, K)
    S > T || throw(ArgumentError("Bond maturity S=$S must be greater than option expiry T=$T"))
    a, σ = m.a, m.σ
    P0T = FinanceCore.discount(m.curve, T)
    P0S = FinanceCore.discount(m.curve, S)

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

function FinanceCore.present_value(m::ShortRate.HullWhite, c::Option.ZCBCall)
    call, _ = _hw_zcb_option_price(m, c.expiry, c.bond_maturity, c.strike)
    return call
end

function FinanceCore.present_value(m::ShortRate.HullWhite, c::Option.ZCBPut)
    _, put = _hw_zcb_option_price(m, c.expiry, c.bond_maturity, c.strike)
    return put
end

# ─── present_value for Caps and Floors ───────────────────────────────────────
#
# A caplet paying max(L(T_{i-1},T_i) - K, 0)·τ at T_i is equivalent to
# (1 + K·τ) puts on a ZCB with maturity T_i, strike 1/(1+K·τ), expiring at T_{i-1}.
#
# Similarly a floorlet = (1 + K·τ) calls on a ZCB.

function FinanceCore.present_value(m::ShortRate.HullWhite, c::Option.Cap)
    K = c.strike
    freq = c.frequency isa FinanceCore.Frequency ? c.frequency.frequency : c.frequency
    τ = 1.0 / freq
    # Payment dates: τ, 2τ, ..., maturity
    # Caplet i: reset at T_{i-1}, pays at T_i
    # First caplet (reset at 0, pay at τ) is typically excluded (rate already known)
    n_periods = round(Int, c.maturity * freq)
    K_bond = 1.0 / (1.0 + K * τ)
    total = 0.0
    for i in 2:n_periods
        T_reset = (i - 1) * τ   # option expiry = reset date
        T_pay   = i * τ         # bond maturity = payment date
        _, put = _hw_zcb_option_price(m, T_reset, T_pay, K_bond)
        total += (1.0 + K * τ) * put
    end
    return total
end

function FinanceCore.present_value(m::ShortRate.HullWhite, c::Option.Floor)
    K = c.strike
    freq = c.frequency isa FinanceCore.Frequency ? c.frequency.frequency : c.frequency
    τ = 1.0 / freq
    n_periods = round(Int, c.maturity * freq)
    K_bond = 1.0 / (1.0 + K * τ)
    total = 0.0
    for i in 2:n_periods
        T_reset = (i - 1) * τ
        T_pay   = i * τ
        call, _ = _hw_zcb_option_price(m, T_reset, T_pay, K_bond)
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

function FinanceCore.present_value(m::ShortRate.HullWhite, c::Option.Swaption)
    T0 = c.expiry
    freq = c.frequency isa FinanceCore.Frequency ? c.frequency.frequency : c.frequency
    τ = 1.0 / freq
    coupon = c.strike

    # Payment dates of the underlying swap
    n_payments = round(Int, (c.swap_maturity - T0) * freq)
    payment_times = [T0 + i * τ for i in 1:n_payments]

    # Step 1: Find r* such that the swap has zero value at T0
    # Swap value at T0 given r(T0) = r:
    #   V(r) = 1 - P(T0,Tn;r) - c·τ·∑ P(T0,Ti;r)
    function swap_value(r)
        total = 1.0
        for (i, Ti) in enumerate(payment_times)
            P_Ti = _hw_zcb(m, T0, Ti, r)
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
        Ki = _hw_zcb(m, T0, Ti, r_star)
        if c.payer
            # Payer swaption = sum of ZCB puts
            _, put = _hw_zcb_option_price(m, T0, Ti, Ki)
            price += coupon * τ * put
            if i == length(payment_times)
                price += put  # principal
            end
        else
            # Receiver swaption = sum of ZCB calls
            call, _ = _hw_zcb_option_price(m, T0, Ti, Ki)
            price += coupon * τ * call
            if i == length(payment_times)
                price += call  # principal
            end
        end
    end
    return price
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
