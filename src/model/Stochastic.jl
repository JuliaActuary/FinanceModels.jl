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
using ..FinanceCore: Continuous, rate

"""
    Vasicek(a, b, σ, initial)

Vasicek (1977) mean-reverting short-rate model:

    dr = a(b - r) dt + σ dW

# Arguments
- `a`: speed of mean reversion
- `b`: long-term mean rate (continuous)
- `σ`: volatility
- `initial`: initial short rate `r₀` (a `Rate` object or `Real`)
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

"""
    CoxIngersollRoss(a, b, σ, initial)

Cox-Ingersoll-Ross (1985) mean-reverting short-rate model:

    dr = a(b - r) dt + σ √r dW

# Arguments
- `a`: speed of mean reversion
- `b`: long-term mean rate (continuous)
- `σ`: volatility
- `initial`: initial short rate `r₀` (a `Rate` object or `Real`)
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
function FinanceCore.discount(m::ShortRate.Vasicek, T)
    a, b, σ, r0 = m.a, m.b, m.σ, _initial_rate(m)
    if abs(a) < 1e-12
        B = T
        lnA = -0.5 * σ^2 * T^3 / 3
    else
        B = (1 - exp(-a * T)) / a
        lnA = (B - T) * (a^2 * b - 0.5 * σ^2) / a^2 - σ^2 * B^2 / (4a)
    end
    return exp(lnA - B * r0)
end

# CIR ZCB price: P(0,T) = A(T) exp(-B(T) r₀)
function FinanceCore.discount(m::ShortRate.CoxIngersollRoss, T)
    a, b, σ, r0 = m.a, m.b, m.σ, _initial_rate(m)
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
    n_steps = round(Int, horizon / dt)
    sqrt_dt = sqrt(dt)

    paths = Vector{RatePath}(undef, n_scenarios)
    for i in 1:n_scenarios
        times = Vector{Float64}(undef, n_steps + 1)
        cumulative = Vector{Float64}(undef, n_steps + 1)
        times[1] = 0.0
        cumulative[1] = 0.0
        r = _sim_initial_rate(model)
        for j in 1:n_steps
            Z = _randn(rng)
            r_new = _step(model, r, dt, sqrt_dt, Z)
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

# Standard normal via inverse CDF using erfinv (avoids Distributions.jl)
function _randn(rng::Random.AbstractRNG)
    u = rand(rng)
    # Clamp to avoid ±Inf
    u = clamp(u, 1e-10, 1 - 1e-10)
    return sqrt(2.0) * SpecialFunctions.erfinv(2u - 1)
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
function _step(m::ShortRate.Vasicek, r, dt, sqrt_dt, Z)
    return r + m.a * (m.b - r) * dt + m.σ * sqrt_dt * Z
end

function _step(m::ShortRate.CoxIngersollRoss, r, dt, sqrt_dt, Z)
    r_pos = max(r, 0.0)
    return r + m.a * (m.b - r) * dt + m.σ * sqrt(r_pos) * sqrt_dt * Z
end

function _step(m::ShortRate.HullWhite, r, dt, sqrt_dt, Z)
    # θ(t) ≈ f'(0,t) + a f(0,t) + σ²/(2a)(1 - exp(-2at))
    # For simulation from the initial curve we use a simplified drift
    a, σ = m.a, m.σ
    return r + a * (_sim_initial_rate(m) - r) * dt + σ * sqrt_dt * Z
end

# ─── pv_mc: Monte Carlo expected present value ───────────────────────────────

"""
    pv_mc(model, contract;
          n_scenarios=1000, timestep=1/12, horizon=nothing,
          rng=Random.default_rng())

Estimate the expected present value of `contract` under the stochastic `model`
by averaging `present_value` across simulated scenarios.
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
