abstract type ParametricModel <: AbstractYieldCurve end

"""
    NelsonSiegel(rates::AbstractVector, maturities::AbstractVector; τ_init=1.0)

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
struct NelsonSiegel <: ParametricModel
    β₀
    β₁
    β₂
    τ₁

    function NelsonSiegel(β₀, β₁, β₂, τ₁)
        (τ₁ <= 0) && throw(DomainError("Wrong parameter ranges"))
        return new(β₀, β₁, β₂, τ₁)
    end
end

function Base.zero(ns::NelsonSiegel, t)
    if iszero(t)
        # zero rate is undefined for t = 0
        t += eps()
    end
    Continuous.(ns.β₀ .+ ns.β₁ .* (1.0 .- exp.(-t ./ ns.τ₁)) ./ (t ./ ns.τ₁) .+ ns.β₂ .* ((1.0 .- exp.(-t ./ ns.τ₁)) ./ (t ./ ns.τ₁) .- exp.(-t ./ ns.τ₁)))
end
discount(ns::NelsonSiegel, t) = discount.(zero.(ns,t),t)

function Par
function NelsonSiegel(yields::Vector{T}, maturities::Vector{U}; τ_init=1.0)  where {T<:Real,U<:Real}
    function fit_β(yields,maturities,τ) 
        Δₘ = vcat([maturities[1]], diff(maturities))
        param₀ = [1.0, 0.0, 0.0]
        spot(m, p) = rate.(zero.(NelsonSiegel(p[1], p[2], p[3],only(τ)), m))
        
        return LsqFit.curve_fit(spot, maturities, yields, Δₘ,param₀)
    end

    function β_sum_sq_resid(τ)
        result = fit_β(yields,maturities,τ) 
        return sum(r^2 for r in result.resid)
    end
    r = Optim.optimize(β_sum_sq_resid, [τ_init])

    τ = only(Optim.minimizer(r))

    result = fit_β(yields,maturities,τ) 
    return NelsonSiegel(result.param[1], result.param[2], result.param[3], τ)
end

function NelsonSiegel(yields::Vector{T}, maturities::Vector{U}; τ_init=1.0) where {T<:Rate,U<:Real}
    cont = [convert(Continuous,r) for r in yields]
    return NelsonSiegelSvensson(cont, maturities; τ_init)
end

"""
    NelsonSiegelSvensson(yields::AbstractVector, maturities::AbstractVector; τ_init=[1.0,1.0])

Return the NelsonSiegelSvensson fitted parameters. The rates should be continuous zero spot rates. If `rates` are not `Rate`s, then they will be interpreted as `Continuous` `Rate`s.

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
struct NelsonSiegelSvensson <: ParametricModel
    β₀
    β₁
    β₂
    β₃
    τ₁
    τ₂

    function NelsonSiegelSvensson(β₀, β₁, β₂, β₃, τ₁, τ₂)
        (τ₁ <= 0 || τ₂ <= 0) && throw(DomainError("Wrong parameter ranges"))
        return new(β₀, β₁, β₂, β₃, τ₁, τ₂)
    end
end

function NelsonSiegelSvensson(yields::Vector{T}, maturities::Vector{U}; τ_init=[1.0,1.0]) where {T<:Real,U<:Real}
    function fit_β(yields,maturities,τ) 
        Δₘ = vcat([maturities[1]], diff(maturities))
        param₀ = [1.0, 0.0, 0.0, 0.0]
        spot(m, p) = rate.(zero.(NelsonSiegelSvensson(p[1], p[2], p[3],p[4],first(τ),last(τ)), m))
        
        return LsqFit.curve_fit(spot, maturities, yields, Δₘ,param₀)
    end

    function β_sum_sq_resid(τ)
        result = fit_β(yields,maturities,τ) 
        return sum(r^2 for r in result.resid)
    end

    r = Optim.optimize(β_sum_sq_resid, τ_init)

    τ = Optim.minimizer(r)[[1,2]]

    result = fit_β(yields,maturities,τ) 
    return NelsonSiegelSvensson(result.param[1], result.param[2], result.param[3],result.param[4],  first(τ), last(τ))
end

function NelsonSiegelSvensson(yields::Vector{T}, maturities::Vector{U}; τ_init=[1.0,1.0]) where {T<:Rate,U<:Real}
    cont = [convert(Continuous,r) for r in yields]
    return NelsonSiegelSvensson(cont, maturities; τ_init)
end



function Base.zero(nss::NelsonSiegelSvensson, t)
    if iszero(t)
        # zero rate is undefined for t = 0
        t += eps()
    end
    Continuous.(nss.β₀ .+ nss.β₁ .* (1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .+ nss.β₂ .* ((1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .- exp.(-t ./ nss.τ₁)) .+ nss.β₃ .* ((1.0 .- exp.(-t ./ nss.τ₂)) ./ (t ./ nss.τ₂) .- exp.(-t ./ nss.τ₂)))
end
discount(nss::NelsonSiegelSvensson, t) = discount.(zero.(nss,t),t)
