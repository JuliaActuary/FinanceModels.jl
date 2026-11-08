module FinanceModelsCopulasExt
# Correlated simulation with a Copulas.jl copula in place of a correlation matrix: each step draws one
# sample of uniforms from the copula and maps each through the standard normal quantile to a model's
# standardized innovation.
import FinanceModels
import Copulas
import Random

const FM = FinanceModels

struct CopulaShocks{d, C}
    copula::C
end
# Copulas' eltype is its declared sample type; the quantile also uses Float64 constants.
FM.__shock_eltype(s::CopulaShocks) = promote_type(eltype(s.copula), Float64)
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
    d = length(models)
    length(C) == d || throw(DimensionMismatch("$d models need a $d-dimensional copula, got $(length(C))"))
    return FM.__simulate(identity, models, copula_shocks(C, d), n_scenarios, timestep, horizon, rng)
end

# Gaussian and independent copulas need no uniforms (a Gaussian sample's uniforms can round to 0
# or 1). Independence also needs no identity matrix or factorization. All adapters use core setup.
copula_shocks(C::Copulas.GaussianCopula, d) = FM.__CorrelatedShocks(C.Σ, d)
copula_shocks(::Copulas.IndependentCopula, d) = FM.__IndependentShocks{d}()
copula_shocks(C::Copulas.Copula, d) = CopulaShocks{d, typeof(C)}(C)

end
