abstract type ParametricModel <: AbstractYield end

"""
    NelsonSiegel(yields::AbstractVector, maturities::AbstractVector, τₐₗₗ::Array{Float64, 1}=[0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10])

Return the NelsonSiegel fitted parameters. Please note there must be no 0's in maturities.

    NelsonSiegel(β₀, β₁, β₂, τ₁)

Parameters of Nelson and Siegel (1987) parametric model. 

# Examples

```julia-repl
julia> β₀, β₁, β₂, τ₁ = 0.6, -1.2, -1.9, 3.0
julia> nsm = NelsonSiegelSvensson.NelsonSiegel.(β₀, β₁, β₂, τ₁)
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


function NelsonSiegel(yields::AbstractVector, maturities::AbstractVector, τₐₗₗ::Array{Float64, 1} = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10])
    Δₘ = vcat([maturities[1]], diff(maturities))
    total_resid = Inf
    ns_param = NelsonSiegel(1.0, 0.0, 0.0, 1)

    for τ in τₐₗₗ
        spot(m, param) = rate.(zero.(NelsonSiegel(param[1], param[2], param[3], τ), m))
        param₀ = [1.0, 0.0, 0.0, 1.0]
        res = LsqFit.curve_fit(spot, maturities, yields, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            ns_param = NelsonSiegel(res.param[1], res.param[2], res.param[3], τ)
        end
    end

    return ns_param
end

"""
    NelsonSiegelSvensson(yields::AbstractVector, maturities::AbstractVector, τₐₗₗ::Array{Float64, 1}=[0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10])

Return the NelsonSiegelSvensson fitted parameters. Please note there must be no 0's in maturities.

    NelsonSiegelSvensson(β₀, β₁, β₂, β₃, τ₁, τ₂)

Parameters of Svensson (1994) parametric model. 


# Examples

```julia-repl
julia> β₀, β₁, β₂, β₃, τ₁, τ₂ = 0.6, -1.2, -2.1, 3.0, 1.5
julia> nssm = NelsonSiegelSvensson.NelsonSiegelSvensson.(β₀, β₁, β₂, β₃, τ₁, τ₂)
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

function NelsonSiegelSvensson(yields::AbstractVector, maturities::AbstractVector, τₐₗₗ::Array{Float64, 1} = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10])
    Δₘ = vcat([maturities[1]], diff(maturities))
    total_resid = Inf
    nss_param = NelsonSiegelSvensson(1.0, 0.0, 0.0, 0.0, 1, 1)

    for τ₁ in τₐₗₗ, τ₂ in τₐₗₗ
        spot(m, param) = rate.(zero.(NelsonSiegelSvensson(param[1], param[2], param[3], param[4], τ₁, τ₂), m))
        param₀ = [1.0, 0.0, 0.0, 0.0, 1.0, 1.0]
        res = LsqFit.curve_fit(spot, maturities, yields, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            nss_param = NelsonSiegelSvensson(res.param[1], res.param[2], res.param[3], res.param[4], τ₁, τ₂)
        end
    end

    return nss_param
end


function Base.zero(nss::NelsonSiegelSvensson, t)
    if iszero(t)
        # zero rate is undefined for t = 0
        t += eps()
    end
    Continuous.(nss.β₀ .+ nss.β₁ .* (1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .+ nss.β₂ .* ((1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .- exp.(-t ./ nss.τ₁)) .+ nss.β₃ .* ((1.0 .- exp.(-t ./ nss.τ₂)) ./ (t ./ nss.τ₂) .- exp.(-t ./ nss.τ₂)))
end
discount(nss::NelsonSiegelSvensson, t) = discount.(zero.(nss,t),t)
