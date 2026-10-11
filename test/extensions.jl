@testset "Makie Extension" begin
    using CairoMakie
    cfs = [Cashflow(t + randn(), t) for t in 1:10]

    @test stem(cfs) isa Makie.FigureAxisPlot

    b = Bond.Fixed(0.1, Periodic(2), 20)
    proj = Projection(b, NullModel(), CashflowProjection())
    # a stem plot:
    @test stem(proj) isa Makie.FigureAxisPlot

    @test stem(b) isa Makie.FigureAxisPlot
end

@testset "UnicodePlots Extension" begin
    using UnicodePlots
    c = Yield.Constant(0.04)
    # rich display renders a plot...
    plain = sprint(show, MIME("text/plain"), c)
    @test occursin("Yield Curve", plain)
    # ...but plain `print`/string interpolation must NOT (it used to render a
    # 60-char-wide plot into every log line and interpolated string)
    @test !occursin("Yield Curve", string(c))
    @test !occursin("Yield Curve", sprint(show, c))
    # a simulated path is plotted up to its horizon, beyond which it throws (its display threw)
    v = ShortRate.Vasicek(0.1, 0.03, 0.01, Continuous(0.03))
    for horizon in (5.0, 0.5)
        path = only(simulate(v; n_scenarios = 1, timestep = 0.1, horizon))
        @test occursin("Yield Curve", sprint(show, MIME("text/plain"), path))
    end
end

import Copulas
import Random
# a copula whose samples include an exact 0: its infinite normal quantile must not become a shock
struct EdgeCopula <: Copulas.Copula{2} end
Base.rand(::Random.AbstractRNG, ::EdgeCopula) = [0.0, 0.5]

@testset "Copulas Extension" begin
    using Copulas
    m = Equity.BlackScholesMerton(0.03, 0.01, 0.2)
    v = ShortRate.Vasicek(0.5, 0.03, 0.01, Continuous(0.02))
    dt = 0.25
    same(a, b) = all(i -> a[i][1].interp.u == b[i][1].interp.u && a[i][2](2.0) == b[i][2](2.0), eachindex(a, b))

    # a Gaussian copula is its correlation matrix and the independence copula the identity, in the
    # Tuple and NamedTuple forms; `GaussianCopula(I)` is an `IndependentCopula` in some versions
    Σ = [1.0 0.4; 0.4 1.0]
    Id = [1.0 0.0; 0.0 1.0]
    for (C, ρ) in ((GaussianCopula(copy(Σ)), Σ), (IndependentCopula(2), Id), (GaussianCopula(copy(Id)), Id))
        kw = (; n_scenarios = 3, timestep = dt, horizon = 2.0)
        b = simulate((v, m), ρ; kw..., rng = Random.Xoshiro(1))
        @test same(simulate((v, m), C; kw..., rng = Random.Xoshiro(1)), b)
        nt = simulate((rate = v, stock = m), C; kw..., rng = Random.Xoshiro(1))
        @test keys(first(nt)) == (:rate, :stock)
        @test same([Tuple(s) for s in nt], b)
    end

    # Clayton(2, θ): Kendall's τ of the innovations is θ/(θ + 2), whatever the marginals
    θ = 7.0
    p = only(simulate((m, m), ClaytonCopula(2, θ); n_scenarios = 1, timestep = dt, horizon = 1000.0, rng = Random.Xoshiro(3)))
    function innovations(path)
        x = log.(path.(FinanceModels.simulation_times(path)))
        return (diff(x) .- (0.03 - 0.01 - 0.2^2 / 2) * dt) ./ (0.2 * sqrt(dt))
    end
    z1, z2 = innovations(p[1]), innovations(p[2])
    n = length(z1)
    @test n == 4000
    concordance = 0.0
    for i in 1:n, j in (i + 1):n
        concordance += sign(z1[i] - z1[j]) * sign(z2[i] - z2[j])
    end
    @test abs(concordance / (n * (n - 1) / 2) - θ / (θ + 2)) < 0.03

    @test_throws DimensionMismatch simulate((m, m), ClaytonCopula(3, 2.0); n_scenarios = 1)
    @test_throws DimensionMismatch simulate((m, m), IndependentCopula(3); n_scenarios = 1)

    # the standard normal quantile, accurate in both tails, and never infinite
    @test FinanceModels.__normal_quantile(0.5) == 0.0
    @test FinanceModels.__normal_quantile(0.975) ≈ 1.959963984540054 rtol = 1.0e-14
    for u in (1.0e-300, 1.0e-10, 0.01, 0.3, 0.7, 0.99, 1 - 1.0e-10)
        @test FinanceModels.__ϕ(FinanceModels.__normal_quantile(u)) ≈ u rtol = 1.0e-12
    end
    @test_throws DomainError FinanceModels.__normal_quantile(0.0)
    @test_throws DomainError FinanceModels.__normal_quantile(1.0)
    # a copula sample of exactly 0 throws rather than becoming an infinite shock
    @test_throws DomainError simulate((m, m), EdgeCopula(); n_scenarios = 1, timestep = dt, horizon = 1.0)
end
