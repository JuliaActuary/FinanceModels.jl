using Random
using ForwardDiff

# sample moments, without a Statistics dependency
sample_mean(x) = sum(x) / length(x)
sample_var(x) = (μ = sample_mean(x); sum(abs2, x .- μ) / (length(x) - 1))
sample_std(x) = sqrt(sample_var(x))

# a volatility that depends on strike or maturity has no single path volatility
struct SmileVolatility <: FinanceModels.Volatility.AbstractVolatilityModel end

@testset "Black–Scholes–Merton price paths" begin
    r, q, σ = 0.03, 0.01, 0.2
    m = Equity.BlackScholesMerton(r, q, σ)

    @testset "exact lognormal steps from a unit spot" begin
        p = only(simulate(m; n_scenarios = 1, timestep = 0.25, horizon = 2.0, rng = Random.Xoshiro(1)))
        ts = FinanceModels.simulation_times(p)
        @test ts == collect(0.0:0.25:2.0)
        @test p(0.0) == 1.0
        # a replay of the transition with the same draws
        rng = Random.Xoshiro(1)
        drift, vol = (r - q - σ^2 / 2) * 0.25, σ * sqrt(0.25)
        logS = 0.0
        for t in ts[2:end]
            logS = logS + drift + vol * randn(rng)
            @test p(t) == exp(logS)
        end
        # zero volatility is the deterministic path exp((r - q)t), whether a number or the unfitted
        # Volatility.Constant()
        for σ0 in (0.0, Volatility.Constant())
            p0 = only(simulate(Equity.BlackScholesMerton(r, q, σ0); n_scenarios = 1, timestep = 0.25, horizon = 2.0))
            @test all(t -> p0(t) ≈ exp((r - q) * t), ts)
        end
        # a number and a Volatility.Constant of it give the same paths
        pc = only(simulate(Equity.BlackScholesMerton(r, q, Volatility.Constant(σ)); n_scenarios = 1, timestep = 0.25, horizon = 2.0, rng = Random.Xoshiro(1)))
        @test all(t -> pc(t) == p(t), ts)
        # ESG's daily case: 252 steps cover one year
        pd = only(simulate(m; n_scenarios = 1, timestep = 1 / 252, horizon = 1.0))
        @test length(FinanceModels.simulation_times(pd)) == 253
    end

    @testset "distribution at one coarse step" begin
        n, T = 200_000, 2.0
        paths = simulate(m; n_scenarios = n, timestep = T, horizon = T, rng = Random.Xoshiro(7))
        @test isconcretetype(eltype(paths))
        ST = [p(T) for p in paths]
        logs = log.(ST)
        # log S_T ~ Normal((r - q - σ²/2)T, σ²T), within 4 standard errors
        @test abs(sample_mean(logs) - (r - q - σ^2 / 2) * T) < 4 * σ * sqrt(T / n)
        @test abs(sample_var(logs) - σ^2 * T) < 4 * σ^2 * T * sqrt(2 / (n - 1))
        # an ex-dividend share: E[e^{-rT} S_T] = e^{-qT}; with its dividends reinvested, 1
        x = exp(-r * T) .* ST
        @test abs(sample_mean(x) - exp(-q * T)) < 4 * sample_std(x) / sqrt(n)
        y = exp(-(r - q) * T) .* ST
        @test abs(sample_mean(y) - 1) < 4 * sample_std(y) / sqrt(n)
        # the Monte Carlo call price agrees with the closed form for a unit stock
        K = 1.05
        payoff = exp(-r * T) .* max.(ST .- K, 0)
        @test abs(sample_mean(payoff) - pv(m, Option.EuroCall(CommonEquity(), K, T))) < 4 * sample_std(payoff) / sqrt(n)
    end

    @testset "prices at grid times only" begin
        p = only(simulate(m; n_scenarios = 1, timestep = 0.1, horizon = 0.95, rng = Random.Xoshiro(2)))
        ts = FinanceModels.simulation_times(p)
        # an unaligned horizon is covered by the step to 1.0
        @test last(ts) ≈ 1.0
        for t in ts
            @test p(nextfloat(t)) == p(t)
        end
        for t in (0.05, 0.95, last(ts) + 0.01, -0.1, NaN)
            @test_throws DomainError p(t)
        end
        ts[1] = 99.0
        @test FinanceModels.simulation_times(p)[1] == 0.0
        # not a yield model
        @test !(p isa Yield.AbstractYieldModel)
        @test_throws MethodError discount(p, 0.5)
        @test repr(p) == "PricePath(11 grid times from 0.0 to 1.0)"
    end

    @testset "invalid inputs" begin
        for (rr, qq, σσ) in ((r, q, -0.2), (r, q, NaN), (r, q, Inf), (-Inf, q, σ), (r, Inf, σ), (NaN, q, σ))
            @test_throws ArgumentError simulate(Equity.BlackScholesMerton(rr, qq, σσ); n_scenarios = 1)
        end
        @test_throws MethodError simulate(Equity.BlackScholesMerton(r, q, SmileVolatility()); n_scenarios = 1)
    end

    @testset "derivatives" begin
        # pathwise ∂S_T/∂r = T·S_T
        ST(rr) = only(simulate(Equity.BlackScholesMerton(rr, q, σ); n_scenarios = 1, timestep = 0.25, horizon = 2.0, rng = Random.Xoshiro(3)))(2.0)
        @test ForwardDiff.derivative(ST, r) ≈ 2.0 * ST(r) rtol = 1.0e-12
        # a model with derivative parameters draws the same shocks
        rng_f, rng_d = Random.Xoshiro(5), Random.Xoshiro(5)
        pf = only(simulate(m; n_scenarios = 1, timestep = 0.25, horizon = 2.0, rng = rng_f))
        pd = only(simulate(Equity.BlackScholesMerton(ForwardDiff.Dual{Nothing}(r, 1.0), q, σ); n_scenarios = 1, timestep = 0.25, horizon = 2.0, rng = rng_d))
        @test rng_f == rng_d
        @test all(t -> ForwardDiff.value(pd(t)) == pf(t), FinanceModels.simulation_times(pf))
        # a derivative through σ = 0: σ = (x - y)² is never negative, but its dual numbers at
        # x = y = 0 carry signed partials, which the sign check must not read
        Sσ(s) = only(simulate(Equity.BlackScholesMerton(r, q, s); n_scenarios = 1, timestep = 0.25, horizon = 2.0, rng = Random.Xoshiro(3)))(2.0)
        mixed = ForwardDiff.derivative(x -> ForwardDiff.derivative(y -> Sσ((x - y)^2), 0.0), 0.0)
        # ∂²/∂x∂y S((x - y)²) at 0 is -2S'(0); its central difference is 2(S(0) - S(h²))/h²
        h = 1.0e-4
        @test mixed ≈ 2 * (Sσ(0.0) - Sσ(h^2)) / h^2 rtol = 1.0e-3
        # a negative seed at σ = 0
        @test ForwardDiff.derivative(s -> Sσ(-s), 0.0) == -ForwardDiff.derivative(Sσ, 0.0)
    end
end

# a custom short-rate model: a constant rate
struct ConstantRateModel <: AbstractStochasticModel
    r::Float64
end
FinanceModels.__sim_initial_rate(m::ConstantRateModel) = m.r
FinanceModels.__step(::ConstantRateModel, r, dt, sqrt_dt, Z, t, ::Nothing, j) = r

@testset "correlated simulation" begin
    v = ShortRate.Vasicek(0.5, 0.03, 0.01, Continuous(0.02))
    hw = ShortRate.HullWhite(0.3, 0.012, Yield.Constant(Continuous(0.03)))
    m = Equity.BlackScholesMerton(0.03, 0.01, 0.2)
    D = ForwardDiff.Dual{Nothing}
    dt = 0.25

    # the standardized innovations each model received, recovered from its recorded path
    function innovations(model::ShortRate.Vasicek, p)
        ou = FinanceModels.__ou_step_params(model.a, model.σ, dt)
        r = p._rates
        return [(r[j + 1] - model.b - (r[j] - model.b) * ou.ϕ) / ou.sd for j in 1:(length(r) - 1)]
    end
    function innovations(model::ShortRate.HullWhite, p)
        r = p._rates
        c = FinanceModels.__sim_cache(model, dt, length(r) - 1)
        return [(r[j + 1] - c.α[j + 1] - (r[j] - c.α[j]) * c.ou.ϕ) / c.ou.sd for j in 1:(length(r) - 1)]
    end
    function innovations(model::Equity.BlackScholesMerton, p)
        x = log.(p.(FinanceModels.simulation_times(p)))
        drift, vol = (model.r - model.q - model.σ^2 / 2) * dt, model.σ * sqrt(dt)
        return [(x[j + 1] - x[j] - drift) / vol for j in 1:(length(x) - 1)]
    end
    sample_cor(x, y) = (μx = sample_mean(x); μy = sample_mean(y); sum((x .- μx) .* (y .- μy)) / sqrt(sum(abs2, x .- μx) * sum(abs2, y .- μy)))

    @testset "one model with ρ = [1;;] is the single-model simulation" begin
        for model in (v, hw, m)
            rng_a, rng_b = Random.Xoshiro(1), Random.Xoshiro(1)
            a = simulate((model,), [1.0;;]; n_scenarios = 3, timestep = dt, horizon = 2.0, rng = rng_a)
            b = simulate(model; n_scenarios = 3, timestep = dt, horizon = 2.0, rng = rng_b)
            @test rng_a == rng_b
            for (x, y) in zip(a, b)
                ts = FinanceModels.simulation_times(y)
                @test model isa Equity.BlackScholesMerton ? all(t -> only(x)(t) == y(t), ts) :
                    (only(x).interp.u == y.interp.u && only(x)._rates == y._rates)
            end
        end
    end

    @testset "innovations are F z" begin
        ρ = [1 0.6 -0.3; 0.6 1 0.2; -0.3 0.2 1]
        models = (v, hw, m)
        F = FinanceModels.__correlation_factor(ρ, 3)
        @test F * F' ≈ ρ
        s = only(simulate(models, ρ; n_scenarios = 1, timestep = dt, horizon = 5.0, rng = Random.Xoshiro(3)))
        rec = map(innovations, models, s)
        rng = Random.Xoshiro(3)
        for j in eachindex(rec[1])
            z = [randn(rng) for _ in 1:3]
            Zc = F * z
            @test all(i -> isapprox(rec[i][j], Zc[i]; atol = 1.0e-8), 1:3)
        end
        # over many steps their sample correlation is ρ
        paths = simulate(models, ρ; n_scenarios = 2000, timestep = dt, horizon = 25.0, rng = Random.Xoshiro(4))
        recs = [reduce(vcat, [innovations(models[i], s[i]) for s in paths]) for i in 1:3]
        n = length(recs[1])
        @test n == 200_000
        for i in 1:3, k in 1:(i - 1)
            @test abs(sample_cor(recs[i], recs[k]) - ρ[i, k]) < 4 * (1 - ρ[i, k]^2) / sqrt(n)
        end
    end

    @testset "terminal covariance follows the discrete scheme" begin
        # Strong mean reversion and a coarse step: the innovation correlation differs visibly from
        # a Brownian-driver correlation of the same value.
        a, σr, σS, c, T, h = 2.0, 0.01, 0.2, 0.6, 5.0, 1.0
        vr = ShortRate.Vasicek(a, 0.03, σr, Continuous(0.03))
        bs = Equity.BlackScholesMerton(0.03, 0.0, σS)
        n = 100_000
        paths = simulate((vr, bs), [1 c; c 1]; n_scenarios = n, timestep = h, horizon = T, rng = Random.Xoshiro(5))
        rT = [rate(short_rate(s[1], T)) for s in paths]
        lT = [log(s[2](T)) for s in paths]
        prods = (rT .- sample_mean(rT)) .* (lT .- sample_mean(lT))
        cov, se = sample_mean(prods), sample_std(prods) / sqrt(n)
        ϕ = exp(-a * h)
        sd = σr * sqrt((1 - ϕ^2) / (2a))
        discrete = c * sd * σS * sqrt(h) * (1 - ϕ^(T / h)) / (1 - ϕ)
        brownian = c * σr * σS * (1 - exp(-a * T)) / a
        @test abs(cov - discrete) < 4se
        @test abs(cov - brownian) > 4se
    end

    @testset "derivatives in ρ" begin
        payoff(c) = sample_mean(
            [
                discount(s[1], 5.0) * s[2](5.0)
                    for s in simulate((v, m), [1 c; c 1]; n_scenarios = 200, timestep = dt, horizon = 5.0, rng = Random.Xoshiro(6))
            ]
        )
        d1 = ForwardDiff.derivative(payoff, 0.3)
        @test isfinite(d1)
        @test d1 ≈ (payoff(0.3 + 1.0e-6) - payoff(0.3 - 1.0e-6)) / 2.0e-6 rtol = 1.0e-5
        d2 = ForwardDiff.derivative(x -> ForwardDiff.derivative(payoff, x), 0.3)
        @test d2 ≈ (ForwardDiff.derivative(payoff, 0.3 + 1.0e-5) - ForwardDiff.derivative(payoff, 0.3 - 1.0e-5)) / 2.0e-5 rtol = 1.0e-4
        # a correlation carrying derivatives widens the discount integral, nested derivatives too
        cc = D(0.3, 1.0)
        s = only(simulate((v, m), [1 cc; cc 1]; n_scenarios = 1, timestep = dt, horizon = 1.0))
        @test eltype(s[1].interp.u) <: ForwardDiff.Dual
        DD = ForwardDiff.Dual{:outer}
        s2 = only(simulate((v, m), [1 DD(cc, cc); DD(cc, cc) 1]; n_scenarios = 1, timestep = dt, horizon = 1.0))
        @test FinanceModels.__ad_depth(eltype(s2[1].interp.u)) == 2
        # at a correlation of ±1 the derivative doesn't exist: it throws before any draw
        for seed in (1.0, -1.0), sgn in (1, -1)
            rng = Random.Xoshiro(9)
            before = copy(rng)
            ρ1 = [1 D(sgn * 1.0, seed); D(sgn * 1.0, seed) 1]
            @test_throws ArgumentError simulate((m, m), ρ1; n_scenarios = 1, rng)
            @test rng == before
        end
    end

    @testset "factorization" begin
        # correlations of ±1 are valid: identical or mirrored innovations
        p = only(simulate((m, m), [1 1; 1 1]; n_scenarios = 1, timestep = dt, horizon = 2.0, rng = Random.Xoshiro(7)))
        @test all(t -> p[1](t) == p[2](t), FinanceModels.simulation_times(p[1]))
        p = only(simulate((m, m), [1 -1; -1 1]; n_scenarios = 1, timestep = dt, horizon = 2.0, rng = Random.Xoshiro(7)))
        @test innovations(m, p[1]) ≈ -innovations(m, p[2]) atol = 1.0e-12
        # a nearly singular positive-definite matrix takes the plain Cholesky factor
        @test size(FinanceModels.__correlation_factor([1 1 - 1.0e-12; 1 - 1.0e-12 1], 2)) == (2, 2)
        # pivoting reorders the variables: the factor puts its rows back
        ρp = [1.0 1.0 0.0; 1.0 1.0 0.0; 0.0 0.0 1.0]
        F = FinanceModels.__correlation_factor(ρp, 3)
        @test size(F) == (3, 2)
        @test maximum(abs, F * F' - ρp) < 1.0e-12
        p = only(simulate((m, m, m), ρp; n_scenarios = 1, timestep = dt, horizon = 2.0, rng = Random.Xoshiro(8)))
        @test all(t -> p[1](t) == p[2](t), FinanceModels.simulation_times(p[1]))
        @test F[1, :]' * F[3, :] == 0
        # not positive semidefinite
        @test_throws ArgumentError simulate((m, m, m), [1 0.9 -0.9; 0.9 1 0.9; -0.9 0.9 1]; n_scenarios = 1)
        # a singular matrix of another element type has no pivoted factor
        @test_throws ArgumentError simulate((m, m), BigFloat[1 1; 1 1]; n_scenarios = 1)
        # a Float64 correlation leaves a Float32 model's integral Float32
        v32 = ShortRate.Vasicek(0.5f0, 0.03f0, 0.01f0, Continuous(0.02f0))
        s = only(simulate((v32, m), [1 0.5; 0.5 1]; n_scenarios = 1, timestep = dt, horizon = 1.0))
        @test eltype(s[1].interp.u) == Float32
    end

    @testset "validation" begin
        @test_throws DimensionMismatch simulate((v, m), [1 0 0; 0 1 0; 0 0 1]; n_scenarios = 1)
        @test_throws ArgumentError simulate((v, m), [1 NaN; NaN 1]; n_scenarios = 1)
        # a covariance matrix
        @test_throws ArgumentError simulate((v, m), [2 0.5; 0.5 1]; n_scenarios = 1)
        @test_throws ArgumentError simulate((v, m), [1 0.5; 0.4 1]; n_scenarios = 1)
        @test_throws MethodError simulate((v, 0.03), [1 0; 0 1]; n_scenarios = 1)
    end

    @testset "tuples, names and custom models" begin
        ρ = [1 0.4 0; 0.4 1 0; 0 0 1]
        t = simulate((v, m, ConstantRateModel(0.03)), ρ; n_scenarios = 4, timestep = dt, horizon = 1.0, rng = Random.Xoshiro(10))
        @test isconcretetype(eltype(t))
        @test all(s -> discount(s[3], 1.0) ≈ exp(-0.03), t)
        nt = simulate((rate = v, stock = m, flat = ConstantRateModel(0.03)), ρ; n_scenarios = 4, timestep = dt, horizon = 1.0, rng = Random.Xoshiro(10))
        @test isconcretetype(eltype(nt))
        @test keys(first(nt)) == (:rate, :stock, :flat)
        @test all(i -> nt[i].rate.interp.u == t[i][1].interp.u && nt[i].stock(1.0) == t[i][2](1.0), 1:4)
    end

    @testset "CIR in a correlated pair" begin
        # CIR's recorded max(x, 0) loses information where it clips, so the latent state is replayed
        cir = ShortRate.CoxIngersollRoss(0.1, 0.1, 0.5, Continuous(0.05))
        ρ = [1 0.5; 0.5 1]
        s = only(simulate((cir, m), ρ; n_scenarios = 1, timestep = dt, horizon = 5.0, rng = Random.Xoshiro(11)))
        F = FinanceModels.__correlation_factor(ρ, 2)
        rng = Random.Xoshiro(11)
        x = FinanceModels.__sim_initial_rate(cir)
        rates, integral = [max(x, 0)], [0.0]
        for j in 1:20
            Z = F * [randn(rng), randn(rng)]
            x_new = FinanceModels.__step(cir, x, dt, sqrt(dt), Z[1], (j - 1) * dt, nothing, j)
            push!(integral, integral[end] + 0.5 * (max(x, 0) + max(x_new, 0)) * dt)
            push!(rates, max(x_new, 0))
            x = x_new
        end
        @test s[1]._rates == rates
        @test s[1].interp.u == integral
        @test any(iszero, rates)
    end
end
