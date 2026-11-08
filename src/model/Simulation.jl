# ─── The shared simulation loop ──────────────────────────────────────────────
#
# `simulate` steps every model of a scenario together on one grid. Each model has a kernel, built
# once per call: a short-rate model's from its simulation hooks (`__sim_cache`,
# `__sim_initial_rate`, `__sim_eltype`, `__step`, `__observed_rate`), Black–Scholes–Merton's from
# its drift and volatility. A scenario draws one standard normal per model per step, in the order
# scenario, step, model, and feeds each model its shock. A single model is a one-model tuple, so it
# draws exactly the shocks it drew before the loop was shared.

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
    simulation_times(path)

A copy of the time grid of `path`, a [`RatePath`](@ref) or a [`PricePath`](@ref), from 0 to its
last time. For a path from [`simulate`](@ref) it is a new `Vector{Float64}` of multiples of the
timestep, ending at an aligned horizon or at the extra step that covers an unaligned one (see
[`FinanceModels.simulation_steps`](@ref)). [`short_rate`](@ref) and a `PricePath`'s prices answer at
these times only.
"""
simulation_times(path::RatePath) = copy(path.interp.t)

# The index of the grid time that `t` names. The nearest grid time is found first and accepted when
# `t` is within 8eps of it, at the scale of that time or of the first step, whichever is larger: a
# time computed as `k * dt` or read from a range may differ from the grid's by a few ulps, and
# `t = 0` gets the same tolerance as the first step. Any other `t` throws.
function __grid_index(ts::AbstractVector, t)
    n = length(ts)
    i = clamp(searchsortedfirst(ts, t), 1, n)
    k = (i > 1 && abs(t - ts[i - 1]) <= abs(t - ts[i])) ? i - 1 : i
    scale = n > 1 ? max(abs(ts[k]), ts[2] - ts[1]) : abs(ts[k])
    abs(t - ts[k]) <= 8eps(float(scale)) && return k
    where = if isnan(t)
        "not a number"
    elseif t > last(ts)
        "beyond the path's last simulated time $(last(ts))"
    elseif t < first(ts)
        "before the path's first simulated time $(first(ts))"
    else
        j = searchsortedlast(ts, t)
        "between the simulated times $(ts[j]) and $(ts[j + 1])"
    end
    throw(DomainError(t, "the time is $where; the simulation has no state there. `FinanceModels.simulation_times(path)` lists the path's times"))
end

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
__shock_eltype(::__IndependentShocks) = Float64

# Every dependence adapter supplies a shock source and its numeric type. Preparation belongs here,
# so model buffers and the innovations they receive always use the same type contract.
function __simulate(finish::F, models::Tuple, shocks, n_scenarios, timestep, horizon, rng) where {F}
    grid = __sim_grid(horizon, timestep)
    W = __shock_eltype(shocks)
    kernels = map(m -> __sim_kernel(m, grid, W), models)
    return __simulate_paths(finish, kernels, grid, n_scenarios, rng, shocks)
end

# `n_scenarios` scenarios on grid `g`, each the tuple of the kernels' paths passed through `finish`
# (`only` for a single model).
function __simulate_paths(finish::F, kernels::Tuple, g::__SimGrid, n_scenarios::Int, rng, shocks) where {F}
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

# An inline callable keeps the per-model update visible to the compiler without recursively
# splatting the state tuple. Recursive reconstruction allocates heavily for larger model tuples.
struct __PathStep{G}
    j::Int
    grid::G
end
@inline (step::__PathStep)(kernel, buffer, state, Z) =
    __path_step!(kernel, buffer, state, step.j, step.grid, Z)
@inline __step_all(kernels, buffers, states, j, grid, Z) =
    map(__PathStep(j, grid), kernels, buffers, states, Z)

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

function __sim_kernel(m::Equity.BlackScholesMerton, g::__SimGrid, ::Type{W}) where {W}
    r, q, σ = m.r, m.q, Equity.__path_vol(m.σ)
    # Infinite rates give prices of 0 or Inf from the first step on. A negative σ gives the paths
    # of |σ| with the shocks' signs flipped, which reverses their correlation with other models.
    # The sign is read from σ's value: under ForwardDiff 1, a dual number at 0 with a negative
    # partial compares below 0, which would reject a derivative through σ = 0.
    (isfinite(r) && isfinite(q)) ||
        throw(ArgumentError("a simulated BlackScholesMerton needs finite rates, got r = $r and q = $q"))
    (isfinite(σ) && __primal(σ) >= 0) ||
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
    return __simulate(only, (model,), __IndependentShocks{1}(), n_scenarios, timestep, horizon, rng)
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

Each path also records the simulated short rate at its grid times
([`FinanceModels.simulation_times`](@ref)), read with [`short_rate`](@ref).

Paths are repeatable: the same `rng` state, arguments and package versions give the same paths.
A model whose parameters carry ForwardDiff numbers draws the same shocks as its `Float64`
counterpart, so paired sensitivity runs see identical scenarios. A change to which draws feed
which step is recorded in the NEWS.
"""
function simulate(
        model::AbstractStochasticModel;
        n_scenarios::Int = 1000,
        timestep::Real = 1 / 12,
        horizon::Real = 30.0,
        rng::Random.AbstractRNG = Random.default_rng()
    )
    return __simulate(only, (model,), __IndependentShocks{1}(), n_scenarios, timestep, horizon, rng)
end

# ─── Correlated simulation ───────────────────────────────────────────────────

# Shocks correlated through a factor F of the models' correlation matrix (F Fᵀ = ρ, d × rank): each
# step draws d standard normals z, whatever the rank, and model i takes row i of F times z.
struct __CorrelatedShocks{d, M <: AbstractMatrix}
    F::M
end
function __CorrelatedShocks(ρ::AbstractMatrix, d::Int)
    F = __correlation_factor(ρ, d)
    return __CorrelatedShocks{d, typeof(F)}(F)
end
__shock_eltype(s::__CorrelatedShocks) = promote_type(eltype(s.F), Float64)
@inline function __draw(s::__CorrelatedShocks{d}, rng) where {d}
    z = ntuple(_ -> randn(rng), Val(d))
    return ntuple(i -> __factor_row(s.F, i, z), Val(d))
end
# Row i of F times z, starting from F[i, 1]·z[1]: for ρ = [1;;] the shock is z[1] itself.
@inline function __factor_row(F, i, z)
    acc = F[i, 1] * z[1]
    for l in 2:size(F, 2)
        acc += F[i, l] * z[l]
    end
    return acc
end

# A factor F of the correlation matrix ρ of `d` models, with F Fᵀ = ρ. The checks guard against
# silent wrong numbers: a matrix of the wrong size, a covariance matrix (its diagonal is not 1) and
# an asymmetric matrix (Cholesky reads one triangle) would otherwise simulate something else.
function __correlation_factor(ρ::AbstractMatrix, d::Int)
    size(ρ) == (d, d) ||
        throw(DimensionMismatch("$d models need a $d×$d correlation matrix, got $(join(size(ρ), '×'))"))
    tol = 8eps(float(typeof(__primal(zero(eltype(ρ))))))
    for i in 1:d
        abs(ρ[i, i] - 1) <= tol || throw(
            ArgumentError("a correlation matrix has a unit diagonal, got $(ρ[i, i]) at [$i, $i]; scale a covariance matrix to correlations")
        )
        for j in 1:(i - 1)
            abs(ρ[i, j] - ρ[j, i]) <= tol || throw(
                ArgumentError("a correlation matrix is symmetric, got $(ρ[i, j]) at [$i, $j] and $(ρ[j, i]) at [$j, $i]")
            )
        end
    end
    # The loop consumes a dense factor. Normalize storage here too: sparse Cholesky otherwise
    # selects CHOLMOD, whose L is a backend-specific component with a fill-reducing permutation.
    A = float.(Matrix(ρ))
    has_derivatives = __ad_depth(eltype(A)) > 0
    values = has_derivatives ? __primal.(A) : A
    C = cholesky(Symmetric(values); check = false)
    if issuccess(C)
        return has_derivatives ? Matrix(cholesky(Symmetric(A)).L) : Matrix(C.L)
    end
    # At a singular correlation the dual factor can report success with infinite derivatives.
    # Check the values first, before any draw, and differentiate only in the PD interior.
    has_derivatives && throw(
        ArgumentError("correlation sensitivities need a positive-definite correlation matrix; this one is singular or not positive semidefinite")
    )
    # Pivoting gives Pᵀ L Lᵀ P = A: restore the original model order and retain the rank columns.
    # Julia supplies generic pivoted Cholesky for additional numeric types on newer releases.
    P = cholesky(Symmetric(A), RowMaximum(); check = false)
    F = Matrix(P.L)[invperm(P.p), 1:P.rank]
    maximum(abs, F * F' - A) <= 8 * d * eps(eltype(A)) ||
        throw(ArgumentError("a correlation matrix must be positive semidefinite; this one is not"))
    return F
end

"""
    simulate(models::Tuple, ρ::AbstractMatrix;
             n_scenarios=1000, timestep=1/12, horizon=30.0,
             rng=Random.default_rng())
    simulate(models::NamedTuple, ρ::AbstractMatrix; ...)

Simulate several models together with correlated shocks. Each scenario is a tuple of paths, one
per model in order (a `NamedTuple` with the same names for a `NamedTuple` of models): a
[`RatePath`](@ref) for a short-rate model and a [`PricePath`](@ref) for an
`Equity.BlackScholesMerton` model, all on the grid a single model's `simulate` uses. Each model
steps with its own transition: exact for Vasicek, Hull–White and Black–Scholes–Merton, full
truncation for Cox–Ingersoll–Ross.

`ρ` is the correlation between the standardized innovations the models receive at each step: a
step draws `d` standard normals `z` and gives model `i` the `i`-th entry of `F z`, where
`F Fᵀ = ρ`. It is not, in general, the Brownian-driver correlation of an exact continuous-time
joint model. An exact transition weights a step's Brownian increments by when they occur: mean
reversion discounts the earlier ones, at its own speed, and a model without mean reversion weights
them equally. Two models that weight them the same way, such as two Vasicek models with the same
speed, have innovations with their drivers' correlation. Models that weight them differently have
innovations correlated less strongly than their drivers, by an amount that grows with the
difference in weighting and with the timestep. The dependence of the models' values at a later time
follows from the models and the timestep.

This is a scenario generator, not a hybrid pricing model: a Black–Scholes–Merton model keeps its
constant-rate drift, so a `RatePath`'s discount factor times a `PricePath` is not a martingale.

`ρ` must be `d × d`, finite, with a unit diagonal and symmetric (within `8eps`), or a
`DimensionMismatch` or `ArgumentError` is thrown. A positive-definite `ρ` is factored by Cholesky.
A singular positive-semidefinite one, with correlations of ±1, is factored by pivoted Cholesky
and accepted when `F Fᵀ` reproduces it within `8 d eps`; any other matrix throws. Supported element
types for pivoted Cholesky follow Julia's implementation (`Float32`/`Float64`, plus generic types
such as `BigFloat` on Julia 1.12 and later). A `ρ` carrying ForwardDiff derivatives must be positive
definite, because the derivatives at a correlation of ±1 don't exist.

A scenario draws `d` standard normals per step, in the order scenario, step, model. A single model
with `ρ = [1;;]` gives exactly the paths of `simulate(model; rng)`.

With Copulas.jl loaded, a `d`-dimensional copula can take the place of `ρ`,
`simulate(models, copula; ...)`. Each step draws one sample of uniforms from the copula and maps
each through the standard normal quantile to a model's innovation, so the copula is the dependence
of each step's innovations, not of the models' later values. A `GaussianCopula` is simulated as its
correlation matrix and an `IndependentCopula` as the identity, exactly as the matrix method. A
copula sample of exactly 0 or 1 throws a `DomainError` rather than giving an infinite shock.
For a custom copula, its declared `eltype` determines the path buffers; a sampler that returns
ForwardDiff numbers can pass pathwise parameter derivatives through the paths.
"""
function simulate(
        models::Tuple, ρ::AbstractMatrix;
        n_scenarios::Int = 1000,
        timestep::Real = 1 / 12,
        horizon::Real = 30.0,
        rng::Random.AbstractRNG = Random.default_rng()
    )
    return __simulate(identity, models, __CorrelatedShocks(ρ, length(models)), n_scenarios, timestep, horizon, rng)
end

# Shape adaptation is shared by every dependence source, including package extensions.
simulate(models::NamedTuple{names}, dependence; kw...) where {names} =
    map(NamedTuple{names}, simulate(Tuple(models), dependence; kw...))
