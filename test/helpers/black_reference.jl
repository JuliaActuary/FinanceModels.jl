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
