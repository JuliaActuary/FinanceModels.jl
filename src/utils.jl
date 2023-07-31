N(x) = Distributions.cdf(Distributions.Normal(), x)

function d1(S, K, τ, r, σ, q)
    return (log(S / K) + (r - q + σ^2 / 2) * τ) / (σ * √(τ))
end

function d2(S, K, τ, r, σ, q)
    return d1(S, K, τ, r, σ, q) - σ * √(τ)
end

"""
    eurocall(;S=1.,K=1.,τ=1,r,σ,q=0.)

Calculate the Black-Scholes implied option price for a european call, where:

- `S` is the current asset price
- `K` is the strike or exercise price
- `τ` is the time remaining to maturity (can be typed with \\tau[tab])
- `r` is the continuously compounded risk free rate
- `σ` is the (implied) volatility (can be typed with \\sigma[tab])
- `q` is the continuously paid dividend rate

Rates should be input as rates (not percentages), e.g.: `0.05` instead of `5` for a rate of five percent.

!!! Experimental: this function is well-tested, but the derivatives functionality (API) may change in a future version of ActuaryUtilities.

# Extended Help

This is the same as the formulation presented in the [dividend extension of the BS model in Wikipedia](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model#Black%E2%80%93Scholes_equation).

## Other general comments:

- Swap/OIS curves are generally better sources for `r` than government debt (e.g. US Treasury) due to the collateralized nature of swap instruments.
- (Implied) volatility is characterized by a curve that is a function of the strike price (among other things), so take care when using 
- Yields.jl can assist with converting rates to continuously compounded if you need to perform conversions.

"""
function eurocall(; S=1.0, K=1.0, τ=1, r, σ, q=0.0)
    iszero(τ) && return max(zero(S), S - K)
    d₁ = d1(S, K, τ, r, σ, q)
    d₂ = d2(S, K, τ, r, σ, q)
    return (N(d₁) * S * exp(τ * (r - q)) - N(d₂) * K) * exp(-r * τ)
end

"""
    europut(;S=1.,K=1.,τ=1,r,σ,q=0.)

Calculate the Black-Scholes implied option price for a european call, where:

- `S` is the current asset price
- `K` is the strike or exercise price
- `τ` is the time remaining to maturity (can be typed with \\tau[tab])
- `r` is the continuously compounded risk free rate
- `σ` is the (implied) volatility (can be typed with \\sigma[tab])
- `q` is the continuously paid dividend rate

Rates should be input as rates (not percentages), e.g.: `0.05` instead of `5` for a rate of five percent.


!!! Experimental: this function is well-tested, but the derivatives functionality (API) may change in a future version of ActuaryUtilities.

# Extended Help

This is the same as the formulation presented in the [dividend extension of the BS model in Wikipedia](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model#Black%E2%80%93Scholes_equation).

## Other general comments:

- Swap/OIS curves are generally better sources for `r` than government debt (e.g. US Treasury) due to the collateralized nature of swap instruments.
- (Implied) volatility is characterized by a curve that is a function of the strike price (among other things), so take care when using 
- Yields.jl can assist with converting rates to continuously compounded if you need to perform conversions.

"""
function europut(; S=1.0, K=1.0, τ=1, r, σ, q=0.0)
    iszero(τ) && return max(zero(S), K - S)
    d₁ = d1(S, K, τ, r, σ, q)
    d₂ = d2(S, K, τ, r, σ, q)
    return (N(-d₂) * K - N(-d₁) * S * exp(τ * (r - q))) * exp(-r * τ)
end