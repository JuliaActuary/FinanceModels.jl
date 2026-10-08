# NS and NSS
## Originally developed by leeyuntien <leeyuntien@gmail.com>

"""
    NelsonSiegel(τ₁, β₀, β₁, β₂)
    NelsonSiegel(τ₁=1.0) # used in fitting


A Nelson-Siegel yield curve model
Parameters of Nelson and Siegel (1987) parametric model, along with default parameter ranges used in the fitting:

- τ₁ controls the location of the hump: `0.0 .. 100.0`
- β₀ represents a long-term interest rate: `-10.0 .. 10.0`
- β₁ represents a time-decay component: `-10.0 .. 10.0`
- β₂ represents a hump: `-10.0 .. 10.0`

# Examples

```julia-repl
julia> τ₁, β₀, β₁, β₂ = 3.0, 0.6, -1.2, -1.9;

julia> nsm = Yield.NelsonSiegel(τ₁, β₀, β₁, β₂);
```

# Extended Help

NelsonSiegel has generally been replaced by NelsonSiegelSvensson, which is a more flexible model.

## References
- https://onriskandreturn.com/2019/12/01/nelson-siegel-yield-curve-model/
- https://www.bis.org/publ/bppdf/bispap25.pdf
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

# Promote mixed argument types to a common type so the inner constructor's
# `where {T}` constraint is satisfied. This is needed for ForwardDiff, which
# passes Dual numbers for the parameters being differentiated while the
# remaining parameters stay as Float64.
function NelsonSiegel(τ₁, β₀, β₁, β₂)
    T = promote_type(typeof(τ₁), typeof(β₀), typeof(β₁), typeof(β₂))
    return NelsonSiegel(convert(T, τ₁), convert(T, β₀), convert(T, β₁), convert(T, β₂))
end

function NelsonSiegel(τ₁ = 1.0)
    return NelsonSiegel(τ₁, 1.0, 0.0, 0.0)
end

# The loadings near t = 0, the decay (1 - e^{-q})/q and the hump decay - e^{-q}. Their closed forms
# cancel there: at q = 5e-17, e^{-q} rounds to 1, so the decay is 0 instead of 1 (the zero rate was off
# by β₁ + β₂). For |q| < 0.1 they come from their Taylor series instead,
#     decay = Σⱼ (-q)ʲ/(j + 1)!,   hump = q·Σⱼ (-q)ʲ (j + 1)/(j + 2)!,
# whose omitted terms are below Float64 rounding there, and which are exact at q = 0 in value and in
# their derivatives (a time derivative at t = 0). Out of line, so that `zero` stays small enough to
# inline; beyond the band the closed form is unchanged.
const __NS_SERIES_Q = 0.1
const __NS_DECAY_COEFFS = ntuple(k -> 1 // factorial(k), 12)
const __NS_HUMP_COEFFS = ntuple(k -> k // factorial(k + 1), 12)
__ns_short(q) = abs(q) < __NS_SERIES_Q
function __ns_loadings(q)
    __ns_short(q) && return (__evalpoly_exact(-q, __NS_DECAY_COEFFS), q * __evalpoly_exact(-q, __NS_HUMP_COEFFS))
    d = -expm1(-q) / q
    return (d, d - exp(-q))
end
@noinline function __ns_zero_short(ns, q)
    d, h = __ns_loadings(q)
    return Continuous(ns.β₀ + ns.β₁ * d + ns.β₂ * h)
end
@noinline function __nss_zero_short(nss, q₁, q₂)
    (d₁, h₁), (_, h₂) = __ns_loadings(q₁), __ns_loadings(q₂)
    return Continuous(nss.β₀ + nss.β₁ * d₁ + nss.β₂ * h₁ + nss.β₃ * h₂)
end

function Base.zero(ns::NelsonSiegel, t)
    q = t / ns.τ₁
    __ns_short(q) && return __ns_zero_short(ns, q)
    # Bind leaf subexpressions (q, e) only — do NOT combine into `decay = (1-e)/q` and
    # write `β·decay`: that reassociates `(β·(1-e))/q → β·((1-e)/q)`, and the sub-ULP
    # gradient shift tips the (documented, highly sensitive) NSS calibration into NaN.
    e = exp(-q)
    return Continuous(ns.β₀ + ns.β₁ * (1 - e) / q + ns.β₂ * ((1 - e) / q - e))
end
# f = (z·t)′ = β₀ + β₁e^{-q} + β₂qe^{-q}: no cancellation, and β₀ at t = ∞
function instantaneous_forward(ns::NelsonSiegel, t)
    __at_infinity(t) && return ns.β₀ + zero(t)
    q = t / ns.τ₁
    e = exp(-q)
    return ns.β₀ + ns.β₁ * e + ns.β₂ * q * e
end
FinanceCore.discount(ns::NelsonSiegel, t) = _discount_from_zero(ns, t)
__log_discount(ns::NelsonSiegel, t) = __zero_log_discount(ns, t)
__log_native(::NelsonSiegel) = true
Base.zero(ns::NelsonSiegel, ts::AbstractArray) = zero.(Ref(ns), ts)
FinanceCore.discount(ns::NelsonSiegel, ts::AbstractArray) = discount.(Ref(ns), ts)

"""
    NelsonSiegelSvensson(τ₁, τ₂, β₀, β₁, β₂, β₃)
    NelsonSiegelSvensson(τ₁=1.0, τ₂=1.0)

Return the NelsonSiegelSvensson yield curve.

Parameters of Svensson (1994) parametric model, along with the default parameter bounds used in the fit routine:

- τ₁ controls the location of the hump: `0.0 .. 100.0`
- τ₂ controls the location of the second hump: `0.0 .. 100.0`
- β₀ represents a long-term interest rate: `-10.0 .. 10.0`
- β₁ represents a time-decay component: `-10.0 .. 10.0`
- β₂ represents a hump: `-10.0 .. 10.0`
- β₃ represents a second hump: `-10.0 .. 10.0`

# Examples

```julia-repl
julia> τ₁, τ₂, β₀, β₁, β₂, β₃ = 1.5, 3.0, 0.6, -1.2, -2.1, 3.0;

julia> nssm = Yield.NelsonSiegelSvensson(τ₁, τ₂, β₀, β₁, β₂, β₃);
```

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

# See NelsonSiegel promotion comment above.
function NelsonSiegelSvensson(τ₁, τ₂, β₀, β₁, β₂, β₃)
    T = promote_type(typeof(τ₁), typeof(τ₂), typeof(β₀), typeof(β₁), typeof(β₂), typeof(β₃))
    return NelsonSiegelSvensson(convert(T, τ₁), convert(T, τ₂), convert(T, β₀), convert(T, β₁), convert(T, β₂), convert(T, β₃))
end

NelsonSiegelSvensson(τ₁ = 1.0, τ₂ = 1.0) = NelsonSiegelSvensson(τ₁, τ₂, 0.0, 0.0, 0.0, 0.0)

function Base.zero(nss::NelsonSiegelSvensson, t)
    q₁ = t / nss.τ₁
    q₂ = t / nss.τ₂
    (__ns_short(q₁) || __ns_short(q₂)) && return __nss_zero_short(nss, q₁, q₂)
    # Bind leaf subexpressions (q, e) only — see the NelsonSiegel `zero` above. Do NOT
    # combine into shared `decay = (1-e)/q` intermediates: reassociating the `β·(1-e)/q`
    # products shifts ForwardDiff gradients enough to tip the sensitive NSS fit into NaN.
    e₁ = exp(-q₁)
    e₂ = exp(-q₂)
    return Continuous(nss.β₀ + nss.β₁ * (1 - e₁) / q₁ + nss.β₂ * ((1 - e₁) / q₁ - e₁) + nss.β₃ * ((1 - e₂) / q₂ - e₂))
end
FinanceCore.discount(nss::NelsonSiegelSvensson, t) = _discount_from_zero(nss, t)
__log_discount(nss::NelsonSiegelSvensson, t) = __zero_log_discount(nss, t)
__log_native(::NelsonSiegelSvensson) = true
function instantaneous_forward(nss::NelsonSiegelSvensson, t)
    __at_infinity(t) && return nss.β₀ + zero(t)
    q₁, q₂ = t / nss.τ₁, t / nss.τ₂
    e₁, e₂ = exp(-q₁), exp(-q₂)
    return nss.β₀ + nss.β₁ * e₁ + nss.β₂ * q₁ * e₁ + nss.β₃ * q₂ * e₂
end
Base.zero(nss::NelsonSiegelSvensson, ts::AbstractArray) = zero.(Ref(nss), ts)
FinanceCore.discount(nss::NelsonSiegelSvensson, ts::AbstractArray) = discount.(Ref(nss), ts)
