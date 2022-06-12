abstract type ParametricModel <: AbstractYield end

"""
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
    #r² no need

    #function NelsonSiegel(β₀, β₁, β₂, τ₁, r²)
    function NelsonSiegel(β₀, β₁, β₂, τ₁)
        #(τ₁ <= 0 || r² < 0 || r² > 1) && throw(DomainError("Wrong parameter ranges"))
        (τ₁ <= 0) && throw(DomainError("Wrong parameter ranges"))
        #return new(β₀, β₁, β₂, τ₁, r²)
        return new(β₀, β₁, β₂, τ₁)
    end
end

"""
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
    #r² no need

    #function NelsonSiegelSvensson(β₀, β₁, β₂, β₃, τ₁, τ₂, r²)
    function NelsonSiegelSvensson(β₀, β₁, β₂, β₃, τ₁, τ₂)
        #(τ₁ <= 0 || τ₂ <= 0 || r² < 0 || r² > 1) && throw(DomainError("Wrong parameter ranges"))
        (τ₁ <= 0 || τ₂ <= 0) && throw(DomainError("Wrong parameter ranges"))
        #return new(β₀, β₁, β₂, β₃, τ₁, τ₂, r²)
        return new(β₀, β₁, β₂, β₃, τ₁, τ₂)
    end
end

Base.zero(ns::NelsonSiegel, t) = Continuous.(ns.β₀ .+ ns.β₁ .* (1.0 .- exp.(-t ./ ns.τ₁)) ./ (t ./ ns.τ₁) .+ ns.β₂ .* ((1.0 .- exp.(-t ./ ns.τ₁)) ./ (t ./ ns.τ₁) .- exp.(-t ./ ns.τ₁)))
discount(ns::NelsonSiegel, t) = discount.(zero.(ns,t),t)

Base.zero(nss::NelsonSiegelSvensson, t) = Continuous.(nss.β₀ .+ nss.β₁ .* (1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .+ nss.β₂ .* ((1.0 .- exp.(-t ./ nss.τ₁)) ./ (t ./ nss.τ₁) .- exp.(-t ./ nss.τ₁)) .+ nss.β₃ .* ((1.0 .- exp.(-t ./ nss.τ₂)) ./ (t ./ nss.τ₂) .- exp.(-t ./ nss.τ₂)))
discount(nss::NelsonSiegelSvensson, t) = discount.(zero.(nss,t),t)

#=""" 
    est_ns_params(swq::Vector{SwapQuote}, τₐₗₗ::Array{Float, 1})

Return the NelsonSiegel fitted parameters. Please note there must be no 0's in maturities.
"""
function est_ns_params(swq::Vector{SwapQuote}, τₐₗₗ::Array{Float64, 1} = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10]) where {Q<:ObservableQuote}
    # τₐₗₗ = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10] move to parameter to allow more freedom
    yields = [q.yield / q.frequency for q in swq] # assume simple interests
    maturities = [q.maturity for q in swq]
    Δₘ = vcat([maturities[1]], diff(maturities))
    #sp = sum(prices .* prices)
    total_resid = Inf
    ns_param = NelsonSiegel(1.0, 0.0, 0.0, 1)

    for τ in τₐₗₗ
        spot(m, param) = zero_disc(NelsonSiegel(param[1], param[2], param[3], param[4]), m)
        param₀ = [1.0, 0.0, 0.0, 1.0]
        res = LsqFit.curve_fit(spot, maturities, yields, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            #r² = 1.0 - sr / sp
            ns_param = NelsonSiegel(res.param[1], res.param[2], res.param[3], τ)
        end
    end

    return ns_param
end=#

""" 
    est_ns_params(yields::AbstractVector, maturities::AbstractVector, τₐₗₗ::Array{Float64, 1})

Return the NelsonSiegel fitted parameters. Please note there must be no 0's in maturities.
"""
function est_ns_params(yields::AbstractVector, maturities::AbstractVector, τₐₗₗ::Array{Float64, 1} = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10])
    # τₐₗₗ = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10] move to parameter to allow more freedom
    Δₘ = vcat([maturities[1]], diff(maturities))
    #sp = sum(prices .* prices)
    total_resid = Inf
    ns_param = NelsonSiegel(1.0, 0.0, 0.0, 1)

    for τ in τₐₗₗ
        spot(m, param) = rate.(zero.(NelsonSiegel(param[1], param[2], param[3], τ), m))
        param₀ = [1.0, 0.0, 0.0, 1.0]
        res = LsqFit.curve_fit(spot, maturities, yields, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            #r² = 1.0 - sr / sp
            ns_param = NelsonSiegel(res.param[1], res.param[2], res.param[3], τ)
        end
    end

    return ns_param
end

#=""" 
    est_nss_params(swq::Vector{SwapQuote}, τₐₗₗ::Array{Float64, 1})

Return the NelsonSiegelSvensson fitted parameters. Please note there must be no 0's in maturities.
"""
function est_nss_params(swq::Vector{SwapQuote}, τₐₗₗ::Array{Float64, 1} = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10]) where {Q<:ObservableQuote}
    # τₐₗₗ = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10] move to parameter to allow more freedom
    yields = [q.yield / q.frequency for q in swq] # assume simple interests
    maturities = [q.maturity for q in swq]
    Δₘ = vcat([maturities[1]], diff(maturities))
    #sp = sum(prices .* prices)
    total_resid = Inf
    nss_param = NelsonSiegelSvensson(1.0, 0.0, 0.0, 0.0, 1, 1)

    for τ₁ in τₐₗₗ, τ₂ in τₐₗₗ
        spot(m, param) = zero_disc(NelsonSiegelSvensson(param[1], param[2], param[3], param[4], τ₁, τ₂), m)
        param₀ = [1.0, 0.0, 0.0, 0.0, 1.0, 1.0]
        res = LsqFit.curve_fit(spot, maturities, yields, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            #r² = 1.0 - sr / sp
            nss_param = NelsonSiegelSvensson(res.param[1], res.param[2], res.param[3], res.param[4], τ₁, τ₂)
        end
    end

    return nss_param
end=#

""" 
    est_nss_params(yields::AbstractVector, maturities::AbstractVector, τₐₗₗ::Array{Float64, 1})

Return the NelsonSiegelSvensson fitted parameters. Please note there must be no 0's in maturities.
"""
function est_nss_params(yields::AbstractVector, maturities::AbstractVector, τₐₗₗ::Array{Float64, 1} = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10])
    # τₐₗₗ = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10] more to parameter to allow more freedom
    Δₘ = vcat([maturities[1]], diff(maturities))
    #sp = sum(prices .* prices)
    total_resid = Inf
    nss_param = NelsonSiegelSvensson(1.0, 0.0, 0.0, 0.0, 1, 1)

    for τ₁ in τₐₗₗ, τ₂ in τₐₗₗ
        spot(m, param) = rate.(zero.(NelsonSiegelSvensson(param[1], param[2], param[3], param[4], τ₁, τ₂), m))
        param₀ = [1.0, 0.0, 0.0, 0.0, 1.0, 1.0]
        res = LsqFit.curve_fit(spot, maturities, yields, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            #r² = 1.0 - sr / sp
            nss_param = NelsonSiegelSvensson(res.param[1], res.param[2], res.param[3], res.param[4], τ₁, τ₂)
        end
    end

    return nss_param
end