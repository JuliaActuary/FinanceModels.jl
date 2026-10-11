# ─── The shared simulation loop ──────────────────────────────────────────────
#
# `simulate` steps every model of a scenario together on one grid. Each model has a kernel, built
# once per call: a short-rate model's from its simulation hooks (`__sim_cache`,
# `__sim_initial_rate`, `__sim_eltype`, `__step`, `__observed_rate`), Black–Scholes–Merton's from
# its drift and volatility. A scenario draws one standard normal per model per step, in the order
# scenario, step, model, and feeds each model its shock. A single model is a one-model tuple, so it
# draws exactly the shocks it drew before the loop was shared.

# The time grid of a simulation: `simulation_steps(horizon, timestep).nsteps` steps of `dt`.
struct __SimGrid
    times::Vector{Float64}
    dt::Float64
    sqrt_dt::Float64
    nsteps::Int
end

function __sim_grid(horizon, timestep)
    grid = simulation_steps(horizon, timestep)
    n = grid.nsteps
    dt = Float64(timestep)
    times = Vector{Float64}(undef, n + 1)
    times[1] = 0.0
    for j in 1:n
        times[j + 1] = j * dt
    end
    # An aligned grid ends at the horizon itself: n·dt can round to either side of it (three
    # steps of 0.3 end at 0.8999999999999999, three of 0.1 at 0.30000000000000004), and the path
    # covers exactly the horizon it was asked for. An unaligned grid ends at its extra step.
    grid.aligned && (times[end] = horizon)
    return __SimGrid(times, dt, sqrt(dt), n)
end

# The element type of a path's discount integral: the model's own (`__sim_eltype`), widened only by a
# shock type that carries derivatives, so Float64 shocks leave a Float32 model's integral Float32.
__integral_eltype(T, ::Type{W}) where {W} = __ad_depth(W) > 0 ? promote_type(T, W) : T

# Independent shocks: one standard normal per model per step, drawn in model order.
struct __IndependentShocks{d} end
@inline __draw(::__IndependentShocks{d}, rng) where {d} = ntuple(_ -> randn(rng), Val(d))

# `n_scenarios` scenarios on grid `g`, each the tuple of the kernels' paths passed through `finish`
# (`only` for a single model).
function __simulate_paths(finish, kernels::Tuple, g::__SimGrid, n_scenarios::Int, rng, shocks)
    # `map` (rather than filling a vector of an abstract path type) infers the concrete element
    # type, so downstream pricing loops dispatch statically instead of dynamically per path
    return map(1:n_scenarios) do _
        started = map(k -> __path_start(k, g), kernels)
        buffers, states0 = map(first, started), map(last, started)
        # The first step is taken outside the loop: a state can change type there (a Float64
        # initial rate becomes a dual number when the parameters carry derivatives, a Float32 one
        # becomes Float64), and inside the loop it then keeps one type.
        states = __step_all(kernels, buffers, states0, 1, g, __draw(shocks, rng))
        for j in 2:g.nsteps
            states = __step_all(kernels, buffers, states, j, g, __draw(shocks, rng))
        end
        finish(map((k, b) -> __path_finish(k, b, g), kernels, buffers))
    end
end

# One step of every model, recursing over the tuples: a `map` with a closure over four tuples was
# not inlined for Hull–White under ForwardDiff and took about 1.8 times as long.
@inline __step_all(::Tuple{}, ::Tuple{}, ::Tuple{}, j, g, ::Tuple{}) = ()
@inline __step_all(kernels::Tuple, buffers::Tuple, states::Tuple, j, g, Z::Tuple) = (
    __path_step!(first(kernels), first(buffers), first(states), j, g, first(Z)),
    __step_all(Base.tail(kernels), Base.tail(buffers), Base.tail(states), j, g, Base.tail(Z))...,
)

# ─── Short-rate kernels ──────────────────────────────────────────────────────

struct __RateKernel{M, C, R, T, S}
    model::M
    cache::C
    r0::R
end

function __sim_kernel(m::AbstractStochasticModel, g::__SimGrid, ::Type{W}) where {W}
    cache = __sim_cache(m, g.dt, g.nsteps)
    # Parameters other than the initial rate may carry an AD number type.
    r0 = __sim_initial_rate(m)
    T = __integral_eltype(__sim_eltype(m, r0), W)
    # The recorded rates keep the precision of the states the integral is accumulated from: a Float64
    # shock makes a Float32 model's states after the first step Float64.
    S = promote_type(T, W)
    return __RateKernel{typeof(m), typeof(cache), typeof(r0), T, S}(m, cache, r0)
end

function __path_start(k::__RateKernel{M, C, R, T, S}, g::__SimGrid) where {M, C, R, T, S}
    cumulative = Vector{T}(undef, g.nsteps + 1)
    cumulative[1] = zero(T)
    rates = Vector{S}(undef, g.nsteps + 1)
    rates[1] = __observed_rate(k.model, k.r0)
    return (cumulative, rates), k.r0
end

@inline function __path_step!(k::__RateKernel, (cumulative, rates), r, j, g::__SimGrid, Z)
    m = k.model
    r_new = __step(m, r, g.dt, g.sqrt_dt, Z, (j - 1) * g.dt, k.cache, j)
    # the integral is accumulated from the working states, never from the recorded rates
    cumulative[j + 1] = cumulative[j] + 0.5 * (__observed_rate(m, r) + __observed_rate(m, r_new)) * g.dt
    rates[j + 1] = __observed_rate(m, r_new)
    return r_new
end

function __path_finish(::__RateKernel, (cumulative, rates), g::__SimGrid)
    # no extrapolation: the path is defined only on its simulated grid
    interp = DataInterpolations.LinearInterpolation(
        cumulative, g.times;
        extrapolation = DataInterpolations.ExtrapolationType.None
    )
    return RatePath(interp, rates)
end

# ─── Black–Scholes–Merton price paths ────────────────────────────────────────

"""
    PricePath

A simulated equity price path from [`simulate`](@ref) with an `Equity.BlackScholesMerton` model:
the ex-dividend price relative to its start, so `p(0) == 1`.

`p(t)` is the simulated price at grid time `t`, one of
[`FinanceModels.simulation_times`](@ref)`(p)`. A time within roundoff of one is read as that time, as
for [`short_rate`](@ref). Any other time throws a `DomainError`: the simulation has no price between
or beyond its grid times, so a payoff at 5.03 needs a grid that contains 5.03.

For prices in currency, multiply by the starting spot: `S0 .* p.(FinanceModels.simulation_times(p))`.
One share with its dividends reinvested is worth `exp(q * t) * p(t)`.

A `PricePath` is not a yield model: it has no `discount`. Only `simulate` builds one.
"""
struct PricePath{T <: Real} <: AbstractModel
    _times::Vector{Float64}
    _prices::Vector{T}
    # simulate's constructor: the grid and the price at each of its times
    function PricePath(times::Vector{Float64}, prices::Vector{T}) where {T <: Real}
        length(times) == length(prices) ||
            throw(DimensionMismatch("a PricePath needs one price per grid time, got $(length(prices)) for $(length(times))"))
        return new{T}(times, prices)
    end
end

(p::PricePath)(t) = p._prices[__grid_index(p._times, t)]
simulation_times(p::PricePath) = copy(p._times)
Base.show(io::IO, p::PricePath) =
    print(io, "PricePath(", length(p._times), " grid times from ", first(p._times), " to ", last(p._times), ")")

struct __PriceKernel{D, V, E}
    drift::D  # (r − q − σ²/2)·dt: the log price's change over a step, without its shock
    vol::V    # σ·√dt
end

# The volatility a simulation uses: a number, or a `Volatility.Constant`'s value (zero for the
# unfitted `Volatility.Constant()`). A volatility that depends on strike or maturity has no single
# path volatility, so other types are a `MethodError`; `Equity.volatility` would read them at the
# money.
__path_vol(σ::Real) = σ
__path_vol(v::Volatility.Constant) = v.σ

function __sim_kernel(m::Equity.BlackScholesMerton, g::__SimGrid, ::Type{W}) where {W}
    r, q, σ = m.r, m.q, __path_vol(m.σ)
    # Infinite rates give prices of 0 or Inf from the first step on. A negative σ gives the paths
    # of |σ| with the shocks' signs flipped, which reverses their correlation with other models.
    (isfinite(r) && isfinite(q)) ||
        throw(ArgumentError("a simulated BlackScholesMerton needs finite rates, got r = $r and q = $q"))
    (isfinite(σ) && σ >= 0) ||
        throw(ArgumentError("a simulated BlackScholesMerton needs a finite, non-negative volatility, got σ = $σ"))
    drift = (r - q - σ^2 / 2) * g.dt
    vol = σ * g.sqrt_dt
    E = promote_type(typeof(drift), typeof(vol), W)
    return __PriceKernel{typeof(drift), typeof(vol), E}(drift, vol)
end

function __path_start(::__PriceKernel{D, V, E}, g::__SimGrid) where {D, V, E}
    prices = Vector{E}(undef, g.nsteps + 1)
    prices[1] = one(E)
    return prices, zero(E)
end

# The exact lognormal transition: the log price takes the step's drift and its shock.
@inline function __path_step!(k::__PriceKernel, prices, logS, j, g::__SimGrid, Z)
    logS_new = logS + k.drift + k.vol * Z
    prices[j + 1] = exp(logS_new)
    return logS_new
end

__path_finish(::__PriceKernel, prices, g::__SimGrid) = PricePath(g.times, prices)

"""
    simulate(model::Equity.BlackScholesMerton;
             n_scenarios=1000, timestep=1/12, horizon=30.0,
             rng=Random.default_rng())

Generate `n_scenarios` Black–Scholes–Merton price paths, as [`PricePath`](@ref)s that start at 1. A
path is the ex-dividend price relative to its start, on the grid of
[`FinanceModels.simulation_steps`](@ref)`(horizon, timestep).nsteps` steps that `simulate` uses for
rate models.

Each step is the exact lognormal transition
``S(t + Δt) = S(t)\\exp\\big((r - q - σ^2/2)Δt + σ\\sqrt{Δt}\\,Z\\big)``, so the prices at the grid
times have no discretisation bias at any `timestep`.

- `σ` is the model's number or its `Volatility.Constant`'s value. `Volatility.Constant()` is zero
  volatility, which gives the deterministic path `exp((r - q)t)`. Other volatility models are a
  `MethodError`.
- The drift is the model's constant `r - q`: risk-neutral when `r` is the risk-free rate, real-world
  when `r` is an expected return. It is never a simulated short rate.
- Rates that are not finite, or a volatility that is negative or not finite, throw an
  `ArgumentError`.

Paths are repeatable as for rate models: the same `rng` state, arguments and package versions give
the same paths, and a model whose parameters carry ForwardDiff numbers draws the same shocks as its
`Float64` counterpart.
"""
function simulate(
        model::Equity.BlackScholesMerton;
        n_scenarios::Int = 1000,
        timestep::Real = 1 / 12,
        horizon::Real = 30.0,
        rng::Random.AbstractRNG = Random.default_rng()
    )
    grid = __sim_grid(horizon, timestep)
    kernel = __sim_kernel(model, grid, Float64)
    return __simulate_paths(only, (kernel,), grid, n_scenarios, rng, __IndependentShocks{1}())
end
