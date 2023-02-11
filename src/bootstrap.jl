# bootstrapped class of curve methods
"""
    Bootstrap(;interpolation=QuadraticSpline)
    
This `CurveMethod` object defines the interpolation method to use when bootstrapping the curve. Provided options are `QuadraticSpline()` (the default) and `LinearSpline()`. You may also pass a custom interpolation method with the function signature of `f(xs, ys) -> f(x) -> y`.

If constructing curves and the rates are not `Rate`s (ie you pass a `Vector{Float64}`), then they will be interpreted as `Periodic(1)` `Rate`s, except the [`Par`](@ref) curve, which is interpreted as `Periodic(2)` `Rate`s. [`CMT`](@ref) and [`OIS`](@ref) FinanceCore.CompoundingFrequency assumption depends on the corresponding maturity.

See for more:

- [`Zero`](@ref)
- [`Forward`](@ref)
- [`Par`](@ref)
- [`CMT`](@ref)
- [`OIS`](@ref)
"""
Base.@kwdef struct Bootstrap{T} <: CurveMethod
    interpolation::T = LinearSpline()
end

__default_rate_interpretation(ns,r::T) where {T<:Rate} = r
__default_rate_interpretation(::Type{Bootstrap{T}},r::U) where {T,U<:Real} = Periodic(r,1)

struct BootstrapCurve{T,U,V} <: AbstractYieldCurve
    rates::T
    maturities::U
    zero::V # function time -> continuous zero rate
end
FinanceCore.discount(yc::T, time) where {T<:BootstrapCurve} = exp(-yc.zero(time) * time)

__ratetype(::Type{BootstrapCurve{T,U,V}}) where {T,U,V}= Yields.Rate{Float64, typeof(DEFAULT_COMPOUNDING)}

function (b::Bootstrap)(quotes::Vector{Quote{T,I}}) where {I<:Cashflow}
    continuous_zeros = [rate(-log(q.price)/q.instrument.time) for q in quotes]
    times = [q.instrument.time for q in quotes]
    intp = b.interpolation([0.0;times],[first(continuous_zeros);continuous_zeros])
    return BootstrapCurve(continuous_zeros, times, intp)
end


function Par(b::Bootstrap,rates, maturities)
    rates = __coerce_rate.(rates,Periodic(2))
    return BootstrapCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise semi-annual
        # per Hull 4.7
        bootstrap(rates, maturities, [m <= 1 ? nothing : 1 / r.compounding.frequency for (r, m) in zip(rates, maturities)], b.interpolation)
    )
end

function Forward(b::Bootstrap,rates, maturities)
    rates = __default_rate_interpretation.(typeof(b),rates)
    # convert to zeros and pass to Zero
    disc_v = Vector{Float64}(undef, length(rates))

    v = 1.0

    for (i,r) = enumerate(rates)
        Δt = maturities[i] - (i == 1 ? 0 : maturities[i-1])
        v *= FinanceCore.discount(r, Δt)
        disc_v[i] = v
    end

    z = (1.0 ./ disc_v) .^ (1 ./ maturities) .- 1 # convert disc_v to zero
    return Zero(b,z, maturities)
end

function CMT(b::Bootstrap,rates, maturities)
    rs = map(zip(rates, maturities)) do (r, m)
        if m <= 1
            Rate(r, Periodic(1 / m))
        else
            Rate(r, Periodic(2))
        end
    end

    CMT(b,rs, maturities)
end

function CMT(b::Bootstrap,rates::Vector{T}, maturities) where {T<:Rate}
    return BootstrapCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise semi-annual
        # per Hull 4.7
        bootstrap(rates, maturities, [m <= 1 ? nothing : 0.5 for m in maturities], b.interpolation)
    )
end


function OIS(b::Bootstrap,rates, maturities)
    rs = map(zip(rates, maturities)) do (r, m)
        if m <= 1
            Rate(r, Periodic(1 / m))
        else
            Rate(r, Periodic(4))
        end
    end

    return OIS(b,rs, maturities)
end
function OIS(b::Bootstrap,rates::Vector{<:Rate}, maturities)
    return BootstrapCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise quarterly
        # per Hull 4.7
        bootstrap(rates, maturities, [m <= 1 ? nothing : 1 / 4 for m in maturities], b.interpolation)
    )
end
