using ForwardDiff

@testset "Derivatives" begin

    @testset "Euro Options" begin
        # tested against https://option-price.com/index.php
        params = (S = 1.0, K = 1.0, τ = 1, r = 0.05, σ = 0.25, q = 0.0)

        @test FinanceModels.eurocall(; params...) ≈ 0.12336 atol = 1.0e-5
        @test FinanceModels.europut(; params...) ≈ 0.07459 atol = 1.0e-5

        params = (S = 1.0, K = 1.0, τ = 1, r = 0.05, σ = 0.25, q = 0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.105493 atol = 1.0e-5
        @test FinanceModels.europut(; params...) ≈ 0.086277 atol = 1.0e-5

        params = (S = 1.0, K = 0.5, τ = 1, r = 0.05, σ = 0.25, q = 0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.49494 atol = 1.0e-5
        @test FinanceModels.europut(; params...) ≈ 0.00011 atol = 1.0e-5

        params = (S = 1.0, K = 0.5, τ = 1, r = 0.05, σ = 0.25, q = 0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.49494 atol = 1.0e-5
        @test FinanceModels.europut(; params...) ≈ 0.00011 atol = 1.0e-5

        params = (S = 1.0, K = 0.5, τ = 0, r = 0.05, σ = 0.25, q = 0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.5 atol = 1.0e-5
        @test FinanceModels.europut(; params...) ≈ 0.0 atol = 1.0e-5

        params = (S = 1.0, K = 1.5, τ = 0, r = 0.05, σ = 0.25, q = 0.03)
        @test FinanceModels.eurocall(; params...) ≈ 0.0 atol = 1.0e-5
        @test FinanceModels.europut(; params...) ≈ 0.5 atol = 1.0e-5

        # `r` and `q` may be `Rate`s, converted to their continuous values
        base = (S = 1.0, K = 1.1, τ = 2, σ = 0.25)
        r, q = Periodic(0.05, 1), Periodic(0.01, 4)
        rc, qc = rate(Continuous(r)), rate(Continuous(q))
        for f in (FinanceModels.eurocall, FinanceModels.europut)
            @test f(; base..., r, q) == f(; base..., r = rc, q = qc)
            @test f(; base..., r = Continuous(0.05)) == f(; base..., r = 0.05)
        end
    end
end


@testset "Equity Models" begin
    m = Equity.BlackScholesMerton(0.01, 0.02, 0.15)

    a = Option.EuroCall(CommonEquity(), 1.0, 1.0)
    b = Option.EuroCall(CommonEquity(), 1.0, 2.0)

    @test pv(m, a) ≈ 0.05410094201902403

    ap = Option.EuroPut(CommonEquity(), 1.0, 1.0)
    @test maturity(a) == 1.0
    @test maturity(b) == 2.0
    @test maturity(ap) == 1.0

    # put call parity with continuous dividend yield
    # https://web.ma.utexas.edu/users/mcudina/m339d_lecture6.pdf, pg 2
    @test pv(m, a) - pv(m, ap) ≈ 1.0 * exp(-0.02) - exp(-0.01) * 1.0


    # test fitting volatility
    qs = [
        Quote(0.0541, a),
        Quote(0.072636, b),
    ]
    m = Equity.BlackScholesMerton(0.01, 0.02, Volatility.Constant())
    fit(m, qs)
    @test fit(m, qs).σ ≈ 0.15 atol = 1.0e-4

    # test constructor with Rates
    @test Equity.BlackScholesMerton(Continuous(0.01), Continuous(0.02), 0.15) isa Equity.BlackScholesMerton
    m = Equity.BlackScholesMerton(Periodic(0.01, 1), Periodic(0.02, 1), 0.15)
    # periodic rates should convert to Continuous equivalent that's higher than the Continuous rate
    @test m.r < 0.01
    @test m.q < 0.02


end

# Black's formula in 256-bit arithmetic, from the given (Float64) discounted forward F, discounted
# strike K and volatility v > 0: the prices and the Greeks (∂F, ∂K, ∂v) of the call and the put
function black_reference(F, K, v)
    erfc = FinanceModels.SpecialFunctions.erfc
    return setprecision(BigFloat, 256) do
        Φ(x) = erfc(-x / sqrt(big(2))) / 2
        φ(x) = exp(-x^2 / 2) / sqrt(2 * big(π))
        F, K, v = big(F), big(K), big(v)
        d₊ = log(F / K) / v + v / 2
        d₋ = d₊ - v
        (
            call = F * Φ(d₊) - K * Φ(d₋), put = K * Φ(-d₋) - F * Φ(-d₊),
            Δcall = [Φ(d₊), -Φ(d₋), F * φ(d₊)], Δput = [-Φ(-d₊), Φ(-d₋), F * φ(d₊)],
        )
    end
end

@testset "Black's formula at zero and tiny volatility (#291)" begin
    black = FinanceModels.__black
    D(f, x) = ForwardDiff.derivative(f, x)
    greeks(price, F, K, v) = ForwardDiff.gradient(x -> price(black(x...)), [F, K, v])
    F = 0.97
    # at, near, out of and in the money, and far from it
    strikes = (F, F * (1 + 1.0e-12), F * (1 - 1.0e-12), F * (1 + 1.0e-6), F * (1 - 1.0e-6), F * 1.05, F / 1.05, F * 1.0e-20, F * 1.0e20)
    # 7.0.0's bond options switched to intrinsic values at a volatility of 1e-15, and differentiating
    # Black's formula overflows below about 1e-155
    vols = (nextfloat(0.0), 1.0e-300, 1.0e-200, 1.0e-160, 1.0e-150, 1.0e-20, prevfloat(1.0e-15), 1.0e-15, nextfloat(1.0e-15), 1.0e-10, 1.0e-4, 0.2, 2.0)

    @testset "prices and Greeks against 256-bit references" begin
        for K in strikes, v in vols
            ref = black_reference(F, K, v)
            call, put = black(F, K, v)
            scale = max(F, K)
            @test abs(call - ref.call) <= 4eps() * scale
            @test abs(put - ref.put) <= 4eps() * scale
            gc, gp = greeks(first, F, K, v), greeks(last, F, K, v)
            # the deltas are probabilities and the vega is at most F/√(2π)
            @test maximum(abs, gc - ref.Δcall) <= 1.0e-13
            @test maximum(abs, gp - ref.Δput) <= 1.0e-13
            # put-call parity, for the prices and their derivatives
            @test abs((call - put) - (F - K)) <= 4eps() * scale
            @test maximum(abs, (gc - gp) - [1, -1, 0]) <= 2eps()
        end
    end

    @testset "zero volatility: intrinsic values and one-sided Greeks" begin
        # the limits as v → 0⁺ at fixed F and K: (Φ(d₊), Φ(d₋)) is (1, 1) in the money, (0, 0) out of
        # it and (1/2, 1/2) at the money, with the vega F/√(2π) at the money only
        for (K, Φ₀, vega) in ((F / 1.05, 1.0, 0.0), (F, 0.5, F / sqrt(2π)), (F * 1.05, 0.0, 0.0))
            call, put = black(F, K, 0.0)
            @test call == max(F - K, 0) && put == max(K - F, 0)
            gc, gp = greeks(first, F, K, 0.0), greeks(last, F, K, 0.0)
            @test gc ≈ [Φ₀, -Φ₀, vega] atol = 1.0e-16
            @test gp ≈ [Φ₀ - 1, 1 - Φ₀, vega] atol = 1.0e-16
            # they are the limits of the Greeks at positive volatilities
            @test gc ≈ greeks(first, F, K, 1.0e-300) atol = 1.0e-16
            @test gp ≈ greeks(last, F, K, 1.0e-300) atol = 1.0e-16
        end
        # a forward a rounding error from the strike is out of the money, with no vega
        @test D(v -> first(black(F, nextfloat(F), v)), 0.0) == 0
    end

    @testset "signed and nested derivatives" begin
        for K in (F, F * (1 + 1.0e-12), F * 1.05, F / 1.05), v in (0.0, 1.0e-300, 0.2)
            price(w) = first(black(F, K, w))
            vega = D(price, v)
            @test D(s -> price(v - s), 0.0) == -vega
            # (x - y)² is never negative, and ∂²/∂x∂y price(v + (x - y)²) at x = y = 0 is -2·vega(v): the
            # term in the second derivative of the price is weighted by zero
            @test D(x -> D(y -> price(v + (x - y)^2), 0.0), 0.0) == -2vega
        end
        # second derivatives at v = 0: zero away from the money, and the limits in v alone at it, where
        # gamma is infinite
        gamma(K) = D(f -> D(g -> first(black(g, K, 0.0)), f), F)
        volga(K) = D(v -> D(w -> first(black(F, K, w)), v), 0.0)
        @test gamma(F * 1.05) == 0 && gamma(F / 1.05) == 0
        @test volga(F * 1.05) == 0 && volga(F / 1.05) == 0 && volga(F) == 0
        @test gamma(F) == Inf
        # at ordinary volatilities the Greeks' derivatives are the formula's second derivatives, here
        # from differentiating it in 256-bit arithmetic
        erfc = FinanceModels.SpecialFunctions.erfc
        for (K, v) in ((1.0, 0.2), (0.97, 0.05), (1.3, 1.0), (0.5, 0.3)), (i, price) in enumerate((first, last))
            H = ForwardDiff.hessian(x -> price(black(x...)), [F, K, v])
            ref = setprecision(BigFloat, 256) do
                Φ(x) = erfc(-x / sqrt(big(2))) / 2
                function formula(x)
                    d₊ = log(x[1] / x[2]) / x[3] + x[3] / 2
                    d₋ = d₊ - x[3]
                    return i == 1 ? x[1] * Φ(d₊) - x[2] * Φ(d₋) : x[2] * Φ(-d₋) - x[1] * Φ(-d₊)
                end
                ForwardDiff.hessian(formula, big.([F, K, v]))
            end
            @test maximum(abs, H - ref) <= 4eps() * maximum(abs, ref)
        end
    end

    @testset "number types and domain" begin
        # Float32, BigFloat, integer and rational inputs
        @test first(black(0.97f0, 1.0f0, 0.2f0)) ≈ first(black(0.97, 1.0, 0.2)) rtol = 1.0e-6
        @test black(1, 1 // 2, 0) == (0.5, 0.0)
        @test all(isfinite, black(0.97f0, 0.97f0, 0.0f0))
        ref = black_reference(0.97, 1.0, 0.2)
        @test first(black(big(0.97), big(1.0), big(0.2))) ≈ ref.call rtol = 1.0e-15
        # a negative volatility is not |v|; a NaN propagates
        @test_throws DomainError black(1.0, 1.0, -1.0e-300)
        @test all(isnan, black(1.0, 1.0, NaN))
    end
end

@testset "Black-Scholes at zero and tiny volatility (#291)" begin
    D(f, x) = ForwardDiff.derivative(f, x)
    r, q, τ, S = 0.03, 0.01, 2.0, 1.0
    Fd, Kd = S * exp(-q * τ), K -> K * exp(-r * τ)
    # at the money for the forward S·e^{(r - q)τ}, at the spot, near, in and out of the money
    for K in (S * exp((r - q) * τ), S, S * exp((r - q) * τ) * (1 + 1.0e-12), 0.9, 1.2),
            σ in (0.0, 1.0e-300, 1.0e-200, 1.0e-15, 1.0e-8, 0.25)
        call(s) = FinanceModels.eurocall(; S, K, τ, r, σ = s, q)
        put(s) = FinanceModels.europut(; S, K, τ, r, σ = s, q)
        # the helper's reference, from the discounted forward and strike as eurocall forms them
        if σ > 0
            ref = black_reference(Fd, Kd(K), σ * sqrt(τ))
            @test abs(call(σ) - ref.call) <= 4eps()
            @test abs(put(σ) - ref.put) <= 4eps()
            @test abs(D(call, σ) - ref.Δcall[3] * sqrt(τ)) <= 1.0e-13
            @test abs(D(put, σ) - ref.Δput[3] * sqrt(τ)) <= 1.0e-13
        else
            @test call(σ) == max(Fd - Kd(K), 0) && put(σ) == max(Kd(K) - Fd, 0)
            @test D(call, σ) == (Fd == Kd(K) ? Fd * sqrt(τ) / sqrt(2π) : 0)
        end
        @test abs((call(σ) - put(σ)) - (Fd - Kd(K))) <= 4eps()
        @test D(call, σ) == D(put, σ)
        # Equity.BlackScholesMerton values options with eurocall and europut
        m = Equity.BlackScholesMerton(r, q, Volatility.Constant(σ))
        @test pv(m, Option.EuroCall(CommonEquity(), K, τ)) == call(σ)
        @test pv(m, Option.EuroPut(CommonEquity(), K, τ)) == put(σ)
    end
    # an exact-zero log-moneyness (S = K, r = q) gave NaN at σ = 0
    @test FinanceModels.eurocall(S = 1.0, K = 1.0, τ = 1.0, r = 0.03, q = 0.03, σ = 0.0) == 0
    @test D(s -> FinanceModels.eurocall(S = 1.0, K = 1.0, τ = 1.0, r = 0.03, q = 0.03, σ = s), 0.0) ≈ exp(-0.03) / sqrt(2π)
    # Volatility.Constant() is zero volatility, with the one-sided vega at the money
    m(s) = Equity.BlackScholesMerton(0.02, 0.02, Volatility.Constant(s))
    @test pv(m(0.0), Option.EuroCall(CommonEquity(), 1.0, 1.0)) == pv(Equity.BlackScholesMerton(0.02, 0.02, Volatility.Constant()), Option.EuroCall(CommonEquity(), 1.0, 1.0)) == 0
    @test D(s -> pv(m(s), Option.EuroCall(CommonEquity(), 1.0, 1.0)), 0.0) ≈ exp(-0.02) / sqrt(2π)
    # at τ = 0 a derivative in τ is the forward's drift in the money, zero out of it, and infinite at it
    θ(K) = D(t -> FinanceModels.eurocall(; S, K, τ = t, r, q, σ = 0.2), 0.0)
    @test θ(0.9) ≈ r * 0.9 - q * S
    @test θ(1.2) == 0
    @test θ(S) == Inf
    # a negative volatility gave a negative price at the money
    @test_throws DomainError FinanceModels.eurocall(S = 1.0, K = exp(0.02), τ = 1.0, r = 0.03, q = 0.01, σ = -0.2)
    @test_throws DomainError pv(Equity.BlackScholesMerton(r, q, -0.2), Option.EuroPut(CommonEquity(), 1.0, 1.0))
end
