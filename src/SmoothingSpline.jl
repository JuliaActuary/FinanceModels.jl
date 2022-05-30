"""
    SmoothingSpline(df, β, xs)
Parameters of a smoothing spline model. 
# Examples
```julia-repl
julia> df, β, xs = 3, zeros(df), [1, 2, 3]
julia> ssm = SmoothingSpline.SmoothingSpline.(df, β, xs)
```
"""
mutable struct SmoothingSpline
    df
    β
    xs

    function SmoothingSpline(df, β, xs)
        (df <= 0) && throw(DomainError("Wrong parameter ranges"))
        return new(df, β, xs)
    end
end

""" 
    createBasisMatrix(maturities::AbstractVector, df)
Create the basis function matrix with a degree of freedom.
"""
function createBasisMatrix(maturities::AbstractVector, df)
    @assert df <= length(maturities)

    m = zeros(length(maturities), length(maturities))
    for i in 1:length(maturities)
        m[i, 1] = 1
        for j in 2:df
            m[i, j] = m[i, j - 1] * maturities[i]
        end
        for j in (df + 1):length(maturities)
            m[i, j] = max((maturities[i] - maturities[j - df]) ^ df, 0)
        end
    end

    return m
end

""" 
    est_ss_params(yields::AbstractVector, maturities::AbstractVector, df=3:3)
Return the SmoothingSpline fitted parameters. Reference: https://www.stat.cmu.edu/~ryantibs/advmethods/notes/smoothspline.pdf.
"""
function est_ss_params(yields::AbstractVector, maturities::AbstractVector, df=3:3)
    dis = Inf
    ssm = SmoothingSpline(3, zeros(3), maturities)
    for d in df
        x = createBasisMatrix(maturities, d)
        svd_fac = svd(x)
        gt = svd_fac.V * Diagonal(svd_fac.S)
        for λ in 0:10:1000 # ignore not so important penalty matrix
            β = svd_fac.V * (gt * gt' .+ λ) * Diagonal(svd_fac.S) * svd_fac.U' * yields
            e = yields .- x * β
            t = sum(e .* e) + λ * dot(β, β)
            if t < dis
                dis = t
                ssm.β = β
                ssm.df = d
            end
        end
    end
    return ssm
end

function generateSplinePoints(ssm::SmoothingSpline, data)
    step = (max(data) - min(data)) / 1000
    xs = [min(data) + i * step for i in 1:1000]
    ys = zeros(1000)
    for (i, x) in enumerate(xs)
        ys = ssm.β[1]
        xt = x
        for j = 2:ssm.df
            ys += ssm.β[j] * xt
            xt *= x
        end
        for j = (ssm.df + 1):length(ssm.β)
            ys += ssm.β[j] * max((x - ssm.xs[j - ssm.df]) ^ df, 0)
        end
    end
end