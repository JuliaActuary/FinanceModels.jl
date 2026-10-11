module FinanceModelsCopulasExt
# Correlated simulation with a Copulas.jl copula in place of a correlation matrix: each step draws one
# sample of uniforms from the copula and maps each through the standard normal quantile to a model's
# standardized innovation.
import FinanceModels
import Copulas
import Random
using LinearAlgebra: I

const FM = FinanceModels

struct CopulaShocks{d, C}
    copula::C
end
@inline function FM.__draw(s::CopulaShocks{d}, rng) where {d}
    u = rand(rng, s.copula)
    return ntuple(i -> FM.__normal_quantile(u[i]), Val(d))
end

function FM.simulate(
        models::Tuple, C::Copulas.Copula;
        n_scenarios::Int = 1000,
        timestep::Real = 1 / 12,
        horizon::Real = 30.0,
        rng::Random.AbstractRNG = Random.default_rng()
    )
    return copula_simulate(identity, models, C, n_scenarios, timestep, horizon, rng)
end

function FM.simulate(
        models::NamedTuple{names}, C::Copulas.Copula;
        n_scenarios::Int = 1000,
        timestep::Real = 1 / 12,
        horizon::Real = 30.0,
        rng::Random.AbstractRNG = Random.default_rng()
    ) where {names}
    return copula_simulate(NamedTuple{names}, Tuple(models), C, n_scenarios, timestep, horizon, rng)
end

# A Gaussian copula is its correlation matrix, and the independence copula the identity: both take
# the matrix method, which needs no uniforms (a Gaussian sample's uniforms can round to exactly 0 or
# 1). `GaussianCopula(I)` is an `IndependentCopula` in some Copulas versions.
copula_simulate(finish, models, C::Copulas.GaussianCopula, args...) =
    FM.__simulate_correlated(finish, models, C.Σ, args...)
copula_simulate(finish, models, C::Copulas.IndependentCopula, args...) =
    FM.__simulate_correlated(finish, models, Matrix{Float64}(I, length(C), length(C)), args...)
function copula_simulate(finish, models, C::Copulas.Copula, n_scenarios, timestep, horizon, rng)
    d = length(models)
    length(C) == d || throw(DimensionMismatch("$d models need a $d-dimensional copula, got $(length(C))"))
    grid = FM.__sim_grid(horizon, timestep)
    kernels = map(m -> FM.__sim_kernel(m, grid, Float64), models)
    return FM.__simulate_paths(finish, kernels, grid, n_scenarios, rng, CopulaShocks{d, typeof(C)}(C))
end

end
