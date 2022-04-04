abstract type ParametricModel end

"""
    NelsonSiegel(β₀, β₁, β₂, τ₁, r²)
Parameters of Nelson and Siegel (1987) parametric model. 
# Examples
```julia-repl
julia> β₀, β₁, β₂, τ₁, r² = 0.6, -1.2, -1.9, 3.0, 1.0
julia> nsm = NelsonSiegelSvensson.NelsonSiegel.(β₀, β₁, β₂, τ₁, r²)
```
"""
struct NelsonSiegel <: ParametricModel
    β₀
    β₁
    β₂
    τ₁
    r²

    function NelsonSiegel(β₀, β₁, β₂, τ₁, r²)
        (τ₁ <= 0 || r² < 0 || r² > 1) && throw(DomainError("Wrong parameter ranges"))
        return new(β₀, β₁, β₂, τ₁, r²)
    end
end

"""
    NelsonSiegelSvensson(β₀, β₁, β₂, β₃, τ₁, τ₂, r²)
Parameters of Svensson (1994) parametric model. 
# Examples
```julia-repl
julia> β₀, β₁, β₂, β₃, τ₁, τ₂, r² = 0.6, -1.2, -2.1, 3.0, 1.5, 1.0
julia> nssm = NelsonSiegelSvensson.NelsonSiegelSvensson.(β₀, β₁, β₂, β₃, τ₁, τ₂, r²)
```
"""
struct NelsonSiegelSvensson <: ParametricModel
    β₀
    β₁
    β₂
    β₃
    τ₁
    τ₂
    r²

    function NelsonSiegelSvensson(β₀, β₁, β₂, β₃, τ₁, τ₂, r²)
        (τ₁ <= 0 || τ₂ <= 0 || r² < 0 || r² > 1) && throw(DomainError("Wrong parameter ranges"))
        return new(β₀, β₁, β₂, β₃, τ₁, τ₂, r²)
    end
end

using LsqFit

""" 
    est_ns_params(zcq::Vector{ZeroCouponQuote})
Return the NelsonSiegel fitted parameters.
"""
function est_ns_params(qs::Vector{Q}) where {Q<:ObservableQuote}
    τₐₗₗ = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10]
    prices = [q.price for q in qs]
    maturities = [q.maturity for q in qs]
    Δₘ = vcat([maturities[1]], diff(maturities))
    sp = sum(prices .* prices)
    total_resid = Inf
    ns_param = NelsonSiegel(1.0, 0.0, 0.0, 1, 0.01)

    for τ in τₐₗₗ
        spot(m, param) = param[1] .+ param[2] .* (1.0 .- exp.(-m ./ τ)) ./ (m ./ τ) .+ param[3] .* ((1.0 .- exp.(-m ./ τ)) ./ (m ./ τ) .- exp.(-m ./ τ))
        param₀ = [1.0, 0.0, 0.0]
        res = LsqFit.curve_fit(spot, maturities, prices, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            ns_param = NelsonSiegel(res.param[1], res.param[2], res.param[3], τ, 1.0 - sr / sp)
        end
    end

    return ns_param
end

""" 
    est_ns_params(prices::AbstractVector, maturities::AbstractVector)
Return the NelsonSiegel fitted parameters.
"""
function est_ns_params(prices::AbstractVector, maturities::AbstractVector)
    τₐₗₗ = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10]
    Δₘ = vcat([maturities[1]], diff(maturities))
    sp = sum(prices .* prices)
    total_resid = Inf
    ns_param = NelsonSiegel(1.0, 0.0, 0.0, 1, 0.01)

    for τ in τₐₗₗ
        spot(m, param) = param[1] .+ param[2] .* (1.0 .- exp.(-m ./ τ)) ./ (m ./ τ) .+ param[3] .* ((1.0 .- exp.(-m ./ τ)) ./ (m ./ τ) .- exp.(-m ./ τ))
        param₀ = [1.0, 0.0, 0.0]
        res = LsqFit.curve_fit(spot, maturities, prices, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            ns_param = NelsonSiegel(res.param[1], res.param[2], res.param[3], τ, 1.0 - sr / sp)
        end
    end

    return ns_param
end

""" 
    est_nss_params(zcq::Vector{ZeroCouponQuote})
Return the NelsonSiegelSvensson fitted parameters.
"""
function est_nss_params(qs::Vector{Q}) where {Q<:ObservableQuote}
    τₐₗₗ = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10]
    prices = [q.price for q in qs]
    maturities = [q.maturity for q in qs]
    Δₘ = vcat([maturities[1]], diff(maturities))
    sp = sum(prices .* prices)
    total_resid = Inf
    nss_param = NelsonSiegelSvensson(1.0, 0.0, 0.0, 0.0, 1, 1, 0.01)

    for τ₁ in τₐₗₗ, τ₂ in τₐₗₗ
        spot(m, param) = param[1] .+ param[2] .* (1.0 .- exp.(-m ./ τ₁)) ./ (m ./ τ₁) .+ param[3] .* ((1.0 .- exp.(-m ./ τ₁)) ./ (m ./ τ₁) .- exp.(-m ./ τ₁)) .+ param[4] .* ((1.0 .- exp.(-m ./ τ₂)) ./ (m ./ τ₂) .- exp.(-m ./ τ₂))
        param₀ = [1.0, 0.0, 0.0, 0.0]
        res = LsqFit.curve_fit(spot, maturities, prices, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            nss_param = NelsonSiegelSvensson(res.param[1], res.param[2], res.param[3], res.param[4], τ₁, τ₂, 1.0 - sr / sp)
        end
    end

    return nss_param
end

""" 
    est_nss_params(prices::AbstractVector, maturities::AbstractVector)
Return the NelsonSiegelSvensson fitted parameters.
"""
function est_nss_params(prices::AbstractVector, maturities::AbstractVector)
    τₐₗₗ = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10]
    Δₘ = vcat([maturities[1]], diff(maturities))
    sp = sum(prices .* prices)
    total_resid = Inf
    nss_param = NelsonSiegelSvensson(1.0, 0.0, 0.0, 0.0, 1, 1, 0.01)

    for τ₁ in τₐₗₗ, τ₂ in τₐₗₗ
        spot(m, param) = param[1] .+ param[2] .* (1.0 .- exp.(-m ./ τ₁)) ./ (m ./ τ₁) .+ param[3] .* ((1.0 .- exp.(-m ./ τ₁)) ./ (m ./ τ₁) .- exp.(-m ./ τ₁)) .+ param[4] .* ((1.0 .- exp.(-m ./ τ₂)) ./ (m ./ τ₂) .- exp.(-m ./ τ₂))
        param₀ = [1.0, 0.0, 0.0, 0.0]
        res = LsqFit.curve_fit(spot, maturities, prices, Δₘ, param₀)
        sr = sum(res.resid .* res.resid)
        if sr < total_resid # take the smallest sum of squares of residuals
            total_resid = sr
            nss_param = NelsonSiegelSvensson(res.param[1], res.param[2], res.param[3], res.param[4], τ₁, τ₂, 1.0 - sr / sp)
        end
    end

    return nss_param
end