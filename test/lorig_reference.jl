using ForwardDiff

# M. Lorig, An Introduction to Financial Derivatives: Open Lecture Notes (2026),
# arXiv:2611.01355v1. Cite sections and exercises to avoid differences between
# equation numbering in the HTML and PDF versions.
@testset "Lorig: option time derivatives" begin
    @testset "Section 7.6 and Exercise 6.4: the dividend-adjusted Black-Scholes PDE" begin
        # Exercise 6.4 supplies the dividend adjustment to Section 7.6's PDE:
        # theta + (r-q)*S*delta + sigma^2*S^2*gamma/2 - r*V = 0.
        # Calendar theta is minus the derivative in time to maturity. A nonzero
        # dividend yield checks its contribution through the public pricing API.
        # https://arxiv.org/html/2611.01355v1#Ch7.S6
        # https://arxiv.org/html/2611.01355v1#Ch6.Thmexercise4
        r, q, σ, τ = 0.05, 0.02, 0.2, 0.75
        for price in (FinanceModels.eurocall, FinanceModels.europut), s in (0.8, 1.0, 1.2)
            value(spot, time = τ) = price(; S = spot, K = 1.0, τ = time, r, q, σ)
            delta = ForwardDiff.derivative(value, s)
            gamma = ForwardDiff.derivative(x -> ForwardDiff.derivative(value, x), s)
            theta = -ForwardDiff.derivative(time -> value(s, time), τ)
            @test theta + (r - q) * s * delta + σ^2 * s^2 * gamma / 2 - r * value(s) ≈ 0 atol = 1.0e-13
        end
    end

    @testset "Exercise 9.1: local variance from a volatility term structure" begin
        # With zero interest and dividend rates and sigma(T) independent of
        # strike, Dupire's formula gives f^2 = sigma(T)^2 + 2(T-t)*sigma(T)*sigma'(T).
        # Check time and second strike derivatives of the public pricing API.
        # https://arxiv.org/html/2611.01355v1#Ch9.Thmexercise1
        σ(τ) = 0.2 + 0.05τ
        price(τ, k) = FinanceModels.eurocall(; S = 1.0, K = k, τ, r = 0.0, q = 0.0, σ = σ(τ))
        for τ in (0.25, 0.75), k in (0.8, 1.0, 1.2)
            time_slope = ForwardDiff.derivative(time -> price(time, k), τ)
            strike_curvature = ForwardDiff.derivative(strike -> ForwardDiff.derivative(x -> price(τ, x), strike), k)
            @test time_slope / (k^2 * strike_curvature / 2) ≈ σ(τ)^2 + 2τ * σ(τ) * 0.05 rtol = 1.0e-12
        end
    end
end
