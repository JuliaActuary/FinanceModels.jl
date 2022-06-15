abstract type ParametricModel <: AbstractYieldCurve end
Base.Broadcast.broadcastable(x::T) where {T<:ParametricModel} = Ref(x)
__ratetype(::Type{T}) where {T<:ParametricModel}= Yields.Rate{Float64, typeof(DEFAULT_COMPOUNDING)}


"""
    NelsonSiegel(rates::AbstractVector, maturities::AbstractVector; τ_initial=1.0)

Return the NelsonSiegel fitted parameters. The rates should be zero spot rates. If `rates` are not `Rate`s, then they will be interpreted as `Continuous` `Rate`s.

    NelsonSiegel(β₀, β₁, β₂, τ₁)

Parameters of Nelson and Siegel (1987) parametric model:

- β₀ represents a long-term interest rate
- β₁ represents a time-decay component
- β₂ represents a hump
- τ₁ controls the location of the hump

# Examples

```julia-repl
julia> β₀, β₁, β₂, τ₁ = 0.6, -1.2, -1.9, 3.0
julia> nsm = Yields.NelsonSiegel.(β₀, β₁, β₂, τ₁)

# Extend Help

## References
- https://onriskandreturn.com/2019/12/01/nelson-siegel-yield-curve-model/
- https://www.bis.org/publ/bppdf/bispap25.pdf

```
"""
struct NelsonSiegelCurve{T} <: ParametricModel
    β₀::T
    β₁::T
    β₂::T
    τ₁::T

    function NelsonSiegelCurve(β₀::T, β₁::T, β₂::T, τ₁::T) where {T<:Real}
        (τ₁ <= 0) && throw(DomainError("Wrong parameter ranges"))
        return new{T}(β₀, β₁, β₂, τ₁)
    end
end


"""
    NelsonSiegel(τ_initial)
    NelsonSiegel() # defaults to τ_initial=1.0
    
This parameter set is used to fit the Nelson-Siegel parametric model to given rates. `τ_initial` should be a scalar and is used as the starting τ value in the optimization. The default value for `τ_initial` is 1.0.

When fitting rates using this `YieldCurveFitParameters` object, the Nelson-Siegel model is used. If constructing curves and the rates are not `Rate`s (ie you pass a `Vector{Float64}`), then they will be interpreted as `Continuous` `Rate`s.

See for more:

- [`Zero`](@ref)
- [`Forward`](@ref)
- [`Par`](@ref)
- [`CMT`](@ref)
- [`OIS`](@ref)
"""
struct NelsonSiegel{T} <: YieldCurveFitParameters
    τ_initial::T
end
NelsonSiegel() = NelsonSiegel(1.0)
__default_rate_interpretation(ns::NelsonSiegel,r) = Continuous(r)

function Base.zero(ns::NelsonSiegelCurve, t)
    if iszero(t)
        # zero rate is undefined for t = 0
        t += eps()
    end
    Continuous.(ns.β₀ .+ ns.β₁ .* (1.0 .- exp.(-t ./ ns.τ₁)) ./ (t ./ ns.τ₁) .+ ns.β₂ .* ((1.0 .- exp.(-t ./ ns.τ₁)) ./ (t ./ ns.τ₁) .- exp.(-t ./ ns.τ₁)))
end
discount(ns::NelsonSiegelCurve, t) = discount.(zero.(ns,t),t)


function fit_β(ns::NelsonSiegel,func,yields,maturities,τ) 
    Δₘ = vcat([maturities[1]], diff(maturities))
    param₀ = [1.0, 0.0, 0.0]
    _rate(m, p) = rate.(func.(NelsonSiegelCurve(p[1], p[2], p[3],only(τ)), m))
    
    return LsqFit.curve_fit(_rate, maturities, yields, Δₘ,param₀)
end

function __fit_NS(ns::NelsonSiegel,func,yields,maturities,τ)
    f(τ) = β_sum_sq_resid(ns,func,yields,maturities,τ)
    r = Optim.optimize(f, [ns.τ_initial])

    τ = only(Optim.minimizer(r))

    return τ, fit_β(ns,func,yields,maturities,τ) 
end

function β_sum_sq_resid(ns,func,yields,maturities,τ)
    result = fit_β(ns,func,yields,maturities,τ) 
    return sum(r^2 for r in result.resid)
end

function Zero(ns::NelsonSiegel,yields::T, maturities::U=eachindex(yields))  where {T<:AbstractArray,U<:AbstractVector}
    yields = rate.(__default_rate_interpretation.(ns,yields))
    func = zero
    τ, result = __fit_NS(ns,func,yields,maturities,ns.τ_initial) 
    return NelsonSiegelCurve(result.param[1], result.param[2], result.param[3], τ)
end

function Par(ns::NelsonSiegel,yields::T, maturities::U=eachindex(yields)) where {T<:AbstractArray,U<:AbstractVector}
    yields = rate.(__default_rate_interpretation.(ns,yields))
    func = par
    τ, result = __fit_NS(ns,func,yields,maturities,ns.τ_initial)
    return NelsonSiegelCurve(result.param[1], result.param[2], result.param[3],result.param[4],  first(τ), last(τ))
end

function Forward(ns::NelsonSiegel,yields::T, maturities::U=eachindex(yields)) where {T<:AbstractArray,U<:AbstractVector}
    yields = rate.(__default_rate_interpretation.(ns,yields))
    func = forward
    τ, result = __fit_NS(ns,func,yields,maturities,ns.τ_initial)
    return NelsonSiegelCurve(result.param[1], result.param[2], result.param[3],result.param[4],  first(τ), last(τ))
end


"""
    NelsonSiegelSvensson(yields::AbstractVector, maturities::AbstractVector; τ_initial=[1.0,1.0])

Return the NelsonSiegelSvensson fitted parameters. The rates should be continuous zero spot rates. If `rates` are not `Rate`s, then they will be interpreted as `Continuous` `Rate`s.

When fitting rates using this `YieldCurveFitParameters` object, the Nelson-Siegel model is used. If constructing curves and the rates are not `Rate`s (ie you pass a `Vector{Float64}`), then they will be interpreted as `Continuous` `Rate`s.

See for more:

    - [`Zero`](@ref)
    - [`Forward`](@ref)
    - [`Par`](@ref)
    - [`CMT`](@ref)
    - [`OIS`](@ref)

    NelsonSiegelSvensson(β₀, β₁, β₂, β₃, τ₁, τ₂)

Parameters of Svensson (1994) parametric model:

- β₀ represents a long-term interest rate
- β₁ represents a time-decay component
- β₂ represents a hump
- β₃ represents a second hum
- τ₁ controls the location of the hump 
- τ₁ controls the location of the second hump 


# Examples

```julia-repl
julia> β₀, β₁, β₂, β₃, τ₁, τ₂ = 0.6, -1.2, -2.1, 3.0, 1.5
julia> nssm = NelsonSiegelSvensson.NelsonSiegelSvensson.(β₀, β₁, β₂, β₃, τ₁, τ₂)

## References
- https://onriskandreturn.com/2019/12/01/nelson-siegel-yield-curve-model/
- https://www.bis.org/publ/bppdf/bispap25.pdf

```
"""
struct NelsonSiegelSvenssonCurve{T} <: ParametricModel
    β₀::T
    β₁::T
    β₂::T
    β₃::T
    τ₁::T
    τ₂::T

    function NelsonSiegelSvenssonCurve(β₀::T, β₁::T, β₂::T, β₃::T, τ₁::T, τ₂::T) where {T<:Real}
        (τ₁ <= 0 || τ₂ <= 0) && throw(DomainError("Wrong parameter ranges"))
        return new{T}(β₀, β₁, β₂, β₃, τ₁, τ₂)
    end
end

"""
    NelsonSiegelSvensson(τ_initial) 
    NelsonSiegelSvensson() # defaults to τ_initial=[1.0,1.0]
        
This parameter set is used to fit the Nelson-Siegel parametric model to given rates. `τ_initial` should be a two element vector and is used as the starting τ value in the optimization. The default value for `τ_initial` is [1.0,1.0].

See for more:

- [`Zero`](@ref)
- [`Forward`](@ref)
- [`Par`](@ref)
- [`CMT`](@ref)
- [`OIS`](@ref)
"""
struct NelsonSiegelSvensson{T} <: YieldCurveFitParameters
    τ_initial::T
end
NelsonSiegelSvensson() = NelsonSiegelSvensson([1.0,1.0])
__default_rate_interpretation(ns::NelsonSiegelSvensson,r) = Continuous(r)

function fit_β(ns::NelsonSiegelSvensson,func,yields,maturities,τ) 
    Δₘ = vcat([maturities[1]], diff(maturities))
    param₀ = [1.0, 0.0, 0.0, 0.0]
    _rate(m, p) = rate.(func.(NelsonSiegelSvenssonCurve(p[1], p[2], p[3],p[4],first(τ),last(τ)), m))
    return LsqFit.curve_fit(_rate, maturities, yields, Δₘ,param₀)
end

function __fit_NS(ns::NelsonSiegelSvensson,func,yields,maturities,τ)
    f(τ) = β_sum_sq_resid(ns,func,yields,maturities,τ)
    r = Optim.optimize(f, ns.τ_initial)

    τ = Optim.minimizer(r)[[1,2]]

    return τ, fit_β(ns,func,yields,maturities,τ) 
end


function Zero(ns::NelsonSiegelSvensson,yields::T, maturities::U=eachindex(yields)) where {T<:AbstractVector,U<:AbstractVector}
    yields = rate.(__default_rate_interpretation.(ns,yields))
    func = zero
    τ, result = __fit_NS(ns,func,yields,maturities,ns.τ_initial)
    return NelsonSiegelSvenssonCurve(result.param[1], result.param[2], result.param[3],result.param[4],  first(τ), last(τ))
end

function Par(ns::NelsonSiegelSvensson,yields::T, maturities::U=eachindex(yields)) where {T<:AbstractVector,U<:AbstractVector}
    yields = rate.(__default_rate_interpretation.(ns,yields))
    func = par
    τ, result = __fit_NS(ns,func,yields,maturities,ns.τ_initial)
    return NelsonSiegelSvenssonCurve(result.param[1], result.param[2], result.param[3],result.param[4],  first(τ), last(τ))
end

function Forward(ns::NelsonSiegelSvensson,yields::T, maturities::U=eachindex(yields)) where {T<:AbstractVector,U<:AbstractVector}
    yields = rate.(__default_rate_interpretation.(ns,yields))
    func = forward
    τ, result = __fit_NS(ns,func,yields,maturities,ns.τ_initial)
    return NelsonSiegelSvenssonCurve(result.param[1], result.param[2], result.param[3],result.param[4],  first(τ), last(τ))
end

function Base.zero(nss::NelsonSiegelSvenssonCurve, t)
    if iszero(t)
        # zero rate is undefined for t = 0
        t += eps()
    end
    Continuous.(nss.β₀ .+ nss.β₁ .* (1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .+ nss.β₂ .* ((1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .- exp.(-t ./ nss.τ₁)) .+ nss.β₃ .* ((1.0 .- exp.(-t ./ nss.τ₂)) ./ (t ./ nss.τ₂) .- exp.(-t ./ nss.τ₂)))
end
discount(nss::NelsonSiegelSvenssonCurve, t) = discount.(zero.(nss,t),t)
