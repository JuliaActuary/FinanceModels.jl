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
    end
end
