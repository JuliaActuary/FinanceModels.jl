using ForwardDiff

# M. Lorig, An Introduction to Financial Derivatives: Open Lecture Notes (2026),
# arXiv:2611.01355v1. Each block identifies its example, exercise or equation.
@testset "Lorig: option prices and Greeks" begin
    @testset "Example 2.6.2: a forward after inception" begin
        # https://arxiv.org/html/2611.01355v1#Ch2.S6.Thmtheorem2
        # The paper's S0=100, St=50, r=0.10, T=1 and t=1/2 give -55.13.
        # Replicate the forward using put-call parity (Section 2.10).
        s0, s, r, maturity, t = 100.0, 50.0, 0.10, 1.0, 0.5
        delivery = s0 * exp(r * maturity)
        expected = s - s0 * exp(r * t)
        @test expected ≈ -55.13 atol = 0.005
        for σ in (0.1, 0.4, 0.8)
            call = FinanceModels.eurocall(; S = s, K = delivery, τ = maturity - t, r, q = 0.0, σ)
            put = FinanceModels.europut(; S = s, K = delivery, τ = maturity - t, r, q = 0.0, σ)
            @test call - put ≈ expected rtol = 1.0e-13
            # Contract strikes are normalized by spot in BlackScholesMerton.
            model = Equity.BlackScholesMerton(r, 0.0, σ)
            c = Option.EuroCall(CommonEquity(), delivery / s, maturity - t)
            p = Option.EuroPut(CommonEquity(), delivery / s, maturity - t)
            @test s * (pv(model, c) - pv(model, p)) ≈ expected rtol = 1.0e-13
        end
    end

    @testset "Exercise 2.6: convexity in strike" begin
        # Eqs. (2.76)-(2.79); unequal spacing requires a non-half weight.
        # https://arxiv.org/html/2611.01355v1#Ch2.Thmexercise6
        k1, k2, k3 = 0.7, 1.0, 1.4
        λ = (k3 - k2) / (k3 - k1)
        for r in (0.0, 0.05), σ in (0.1, 0.4), τ in (0.25, 1.0),
                option in (Option.EuroCall, Option.EuroPut)
            model = Equity.BlackScholesMerton(r, 0.0, σ)
            prices = [pv(model, option(CommonEquity(), k, τ)) for k in (k1, k2, k3)]
            @test prices[2] <= λ * prices[1] + (1 - λ) * prices[3]
        end
    end

    @testset "Exercise 7.1: Greeks at three times" begin
        # r=0, sigma=0.2, K=T=1 as specified in the exercise. The ATM price
        # references evaluate Eq. (7.53) as erf(sigma*sqrt(T-t)/(2sqrt(2))).
        # https://arxiv.org/html/2611.01355v1#Ch7.Thmexercise1
        for (t, expected) in ((0.0, 0.07965567455405796), (0.5, 0.05637197779701664), (0.75, 0.039877611676744924))
            τ, σ = 1 - t, 0.2
            call(s, vol = σ) = FinanceModels.eurocall(; S = s, K = 1.0, τ, r = 0.0, q = 0.0, σ = vol)
            @test call(1.0) ≈ expected rtol = 1.0e-13
            for s in (0.8, 1.0, 1.2)
                delta = ForwardDiff.derivative(call, s)
                gamma = ForwardDiff.derivative(x -> ForwardDiff.derivative(call, x), s)
                vega = ForwardDiff.derivative(vol -> call(s, vol), σ)
                @test 0 < delta < 1
                @test gamma > 0
                @test vega ≈ τ * σ * s^2 * gamma rtol = 1.0e-12
            end
        end
    end

    @testset "Section 7.6: the Black-Scholes PDE" begin
        # Eq. (7.48), tested for both prices with a nonzero rate.
        # Calendar theta is minus the derivative in time to maturity.
        # https://arxiv.org/html/2611.01355v1#Ch7.S6
        r, σ, τ = 0.05, 0.2, 0.75
        for price in (FinanceModels.eurocall, FinanceModels.europut), s in (0.8, 1.0, 1.2)
            value(spot, time = τ) = price(; S = spot, K = 1.0, τ = time, r, q = 0.0, σ)
            delta = ForwardDiff.derivative(value, s)
            gamma = ForwardDiff.derivative(x -> ForwardDiff.derivative(value, x), s)
            theta = -ForwardDiff.derivative(time -> value(s, time), τ)
            @test theta + r * s * delta + σ^2 * s^2 * gamma / 2 - r * value(s) ≈ 0 atol = 1.0e-13
        end
    end

    @testset "Exercise 7.2: Greeks through a volatility smile" begin
        # Exercise 7.2 / Eqs. (A.331), (A.335): implied volatility is
        # I(tau, log(K/S)), so its spot and time derivatives enter the Greeks.
        # I below is a chosen smooth smile, with nonzero skew and term slope.
        # https://arxiv.org/html/2611.01355v1#Ch7.Thmexercise2
        smile(τ, x) = 0.2 + 0.03τ + 0.04x + 0.02x^2
        τ, k = 0.75, 1.0
        for s in (0.8, 1.0, 1.2)
            x = log(k / s)
            σ = smile(τ, x)
            value(spot, time, vol) = FinanceModels.eurocall(; S = spot, K = k, τ = time, r = 0.0, q = 0.0, σ = vol)
            delta = ForwardDiff.derivative(spot -> value(spot, τ, σ), s)
            theta = -ForwardDiff.derivative(time -> value(s, time, σ), τ)
            vega = ForwardDiff.derivative(vol -> value(s, τ, vol), σ)
            smile_delta = ForwardDiff.derivative(spot -> value(spot, τ, smile(τ, log(k / spot))), s)
            smile_theta = -ForwardDiff.derivative(time -> value(s, time, smile(time, x)), τ)
            @test smile_delta ≈ delta - vega * (0.04 + 0.04x) / s rtol = 1.0e-12
            @test smile_theta ≈ theta - vega * 0.03 rtol = 1.0e-12
        end
    end

    @testset "Exercise 9.1: local variance from a volatility term structure" begin
        # Dupire's formula, Eqs. (9.41), (A.373): with sigma(T) independent
        # of strike, f^2 = sigma(T)^2 + 2(T-t)*sigma(T)*sigma'(T).
        # This checks time and second strike derivatives of the public price API.
        # https://arxiv.org/html/2611.01355v1#Ch9.Thmexercise1
        σ(τ) = 0.2 + 0.05τ
        price(τ, k) = FinanceModels.eurocall(; S = 1.0, K = k, τ, r = 0.0, q = 0.0, σ = σ(τ))
        for τ in (0.25, 0.75), k in (0.8, 1.0, 1.2)
            time_slope = ForwardDiff.derivative(time -> price(time, k), τ)
            strike_curvature = ForwardDiff.derivative(strike -> ForwardDiff.derivative(x -> price(τ, x), strike), k)
            @test time_slope / (k^2 * strike_curvature / 2) ≈ σ(τ)^2 + 2τ * σ(τ) * 0.05 rtol = 1.0e-12
        end
    end

    @testset "Exercise 11.1: put-call symmetry" begin
        # Eq. (11.29), with G(x)=(x-K)^+ and S=1, gives C(1,K)=K*P(1,1/K).
        # This symmetry requires zero interest and dividend rates.
        # https://arxiv.org/html/2611.01355v1#Ch11.Thmexercise1
        for k in (0.7, 1.0, 1.4), τ in (0.25, 1.0), σ in (0.2, 0.6)
            model = Equity.BlackScholesMerton(0.0, 0.0, σ)
            call = Option.EuroCall(CommonEquity(), k, τ)
            put = Option.EuroPut(CommonEquity(), 1 / k, τ)
            @test pv(model, call) ≈ k * pv(model, put) rtol = 1.0e-12
        end
    end
end
