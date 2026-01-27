# NS and NSS
## Originally developed by leeyuntien <leeyuntien@gmail.com>

"""
    NelsonSiegel(β₀, β₁, β₂, τ₁)
    NelsonSiegel(τ₁=1.0) # used in fitting


A Nelson-Siegel yield curve model 
Parameters of Nelson and Siegel (1987) parametric model, along with default parameter ranges used in the fitting:

- β₀ represents a long-term interest rate: `-10.0 .. 10.0`
- β₁ represents a time-decay component: `-10.0 .. 10.0`
- β₂ represents a hump: `-10.0 .. 10.0`
- τ₁ controls the location of the hump: `0.0 .. 100.0`

# Examples

```julia-repl
julia> β₀, β₁, β₂, τ₁ = 0.6, -1.2, -1.9, 3.0
julia> nsm = Yields.NelsonSiegel.(β₀, β₁, β₂, τ₁)

# Extended Help

NelsonSiegel has generally been replaced by NelsonSiegelSvensson, which is a more flexible model.

## References
- https://onriskandreturn.com/2019/12/01/nelson-siegel-yield-curve-model/
- https://www.bis.org/publ/bppdf/bispap25.pdf

```
"""
struct NelsonSiegel{T} <: AbstractYieldModel
    τ₁::T
    β₀::T
    β₁::T
    β₂::T

    function NelsonSiegel(τ₁::T, β₀::T, β₁::T, β₂::T) where {T}
        (τ₁ <= 0) && throw(DomainError("Wrong tau parameter ranges (must be positive)"))
        return new{T}(τ₁, β₀, β₁, β₂)
    end
end

function NelsonSiegel(τ₁, β₀, β₁, β₂)
    T = promote_type(typeof(τ₁), typeof(β₀), typeof(β₁), typeof(β₂))
    return NelsonSiegel(convert(T, τ₁), convert(T, β₀), convert(T, β₁), convert(T, β₂))
end

function NelsonSiegel(τ₁ = 1.0)
    return NelsonSiegel(τ₁, 1.0, 0.0, 0.0)
end

function Base.zero(ns::NelsonSiegel, t)
    if iszero(t)
        # zero rate is undefined for t = 0
        t += eps()
    end
    return Continuous.(ns.β₀ .+ ns.β₁ .* (1.0 .- exp.(-t ./ ns.τ₁)) ./ (t ./ ns.τ₁) .+ ns.β₂ .* ((1.0 .- exp.(-t ./ ns.τ₁)) ./ (t ./ ns.τ₁) .- exp.(-t ./ ns.τ₁)))
end
FinanceCore.discount(ns::NelsonSiegel, t) = discount.(zero.(ns, t), t)

"""
    NelsonSiegelSvensson(τ₁, τ₂, β₀, β₁, β₂, β₃)
    NelsonSiegelSvensson(τ₁=1.0, τ₂=1.0)

Return the NelsonSiegelSvensson yield curve. The rates should be continuous zero spot rates. If `rates` are not `Rate`s, then they will be interpreted as `Continuous` `Rate`s.

Parameters of Svensson (1994) parametric model, along with the default parameter bounds used in the fit routine:

- τ₁ controls the location of the hump: `0.0 .. 100.0`
- τ₁ controls the location of the second hump: `0.0 .. 100.0`
- β₀ represents a long-term interest rate: `-10.0 .. 10.0`
- β₁ represents a time-decay component: `-10.0 .. 10.0`
- β₂ represents a hump: `-10.0 .. 10.0`
- β₃ represents a second hump: `-10.0 .. 10.0`

# Examples

```julia-repl
julia> β₀, β₁, β₂, β₃, τ₁, τ₂ = 0.6, -1.2, -2.1, 3.0, 1.5
julia> nssm = NelsonSiegelSvensson.NelsonSiegelSvensson.(β₀, β₁, β₂, β₃, τ₁, τ₂)

# Extended Help

Nelson-Siegel-Svensson Pros:

- Simplicity: With only six parameters, the model is quite parsimonious and easy to estimate. It's also easier to interpret and communicate than more complex models.
- Economic Interpretability: Each of the model's components can be given an economic interpretation, with parameters representing long term rate, short term rate, the rates of decay towards the long term rate, and humps in the yield curve.

Nelson-Siegel-Svensson Cons:

- Unusual Curves: NSS makes some assumptions about the shape of the yield curve (e.g. generally has a hump in short to medium term maturities). It might not be the best choice for fitting unusual curves.
- Arbitrage Opportunities: The NSS model does not guarantee absence of arbitrage opportunities. More sophisticated models, like the ones based on no-arbitrage conditions, might provide better pricing accuracy in some contexts.
- Sensitivity: Similar inputs may produce different parameters due to the highly convex, non-linear region to solve for the parameters. Entities like the ECB will partially mitigate this by using the prior business day's parameters as the starting point for the current day's yield curve.

## References
- https://onriskandreturn.com/2019/12/01/nelson-siegel-yield-curve-model/
- https://www.bis.org/publ/bppdf/bispap25.pdf

```
"""
struct NelsonSiegelSvensson{T} <: AbstractYieldModel
    τ₁::T
    τ₂::T
    β₀::T
    β₁::T
    β₂::T
    β₃::T

    function NelsonSiegelSvensson(τ₁::T, τ₂::T, β₀::T, β₁::T, β₂::T, β₃::T) where {T}
        (τ₁ <= 0 || τ₂ <= 0) && throw(DomainError("Wrong tau parameter ranges (must be positive)"))
        return new{T}(τ₁, τ₂, β₀, β₁, β₂, β₃)
    end
end

function NelsonSiegelSvensson(τ₁, τ₂, β₀, β₁, β₂, β₃)
    T = promote_type(typeof(τ₁), typeof(τ₂), typeof(β₀), typeof(β₁), typeof(β₂), typeof(β₃))
    return NelsonSiegelSvensson(convert(T, τ₁), convert(T, τ₂), convert(T, β₀), convert(T, β₁), convert(T, β₂), convert(T, β₃))
end

NelsonSiegelSvensson(τ₁ = 1.0, τ₂ = 1.0) = NelsonSiegelSvensson(τ₁, τ₂, 0.0, 0.0, 0.0, 0.0)

function Base.zero(nss::NelsonSiegelSvensson, t)
    if iszero(t)
        # zero rate is undefined for t = 0
        t += eps()
    end
    return Continuous.(nss.β₀ .+ nss.β₁ .* (1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .+ nss.β₂ .* ((1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .- exp.(-t ./ nss.τ₁)) .+ nss.β₃ .* ((1.0 .- exp.(-t ./ nss.τ₂)) ./ (t ./ nss.τ₂) .- exp.(-t ./ nss.τ₂)))
end
FinanceCore.discount(nss::NelsonSiegelSvensson, t) = discount.(zero.(nss, t), t)
