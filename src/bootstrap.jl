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

function (b::Bootstrap)(quotes::Vector{Quote{T,I}}) where {T,I<:Cashflow}
    continuous_zeros = [-log(q.price)/q.instrument.time for q in quotes]
    times = [q.instrument.time for q in quotes]
    intp = b.interpolation([0.0;times],[first(continuous_zeros);continuous_zeros])
    return BootstrapCurve(continuous_zeros, times, intp)
end

# tried to just use the version without the `guess` argument, but in the solving,
# it was getting stuck where the ForwardDiff.Dual type was different than the Float64 `T`
# and not dispatching
function (b::Bootstrap)(quotes::Vector{Quote{T,I}},guess) where {T,I<:Cashflow}
    continuous_zeros = [-log(q.price)/q.instrument.time for q in quotes]
    continuous_zeros = vcat(continuous_zeros,-log(guess.price)/guess.instrument.time)
    times = [q.instrument.time for q in quotes]
    times = vcat(times,guess.instrument.time)
    intp = b.interpolation([0.0;times],[first(continuous_zeros);continuous_zeros])
    return BootstrapCurve(continuous_zeros, times, intp)
end

function (b::Bootstrap)(quotes::Vector{Quote{T,I}}) where {T,I<:Bond}
    _bootstrap_instrument(b,quotes)
end

function _bootstrap_instrument(bs::Bootstrap,quotes::Vector{Quote{P,I}}) where {P,I<:Bond}
    # use first coupon rate as the initial guess
    maturities = [q.instrument.maturity for q in quotes]
    z = ZCBYield.(zeros(length(quotes)), maturities)

    ## TODO: Optim is looking for minimum, not root. Need to get root finding to work or switch Optim algo
    
    # we have to take the first rate as the starting point
    for (i,q) in enumerate(quotes)
        # it's a discount/premium bond so we need to solve for the rate that works
        b = q.instrument
        
        function root_func(v_guess)
            # z_inner = vcat(z[1:i-1],ZCBYield(v_guess[1],maturities[i]))
            # @show z_inner
            c = bs(z[1:i-1],ZCBYield(v_guess[1],maturities[i])) 
            _pv(c,b) - q.price
        end
        root_func′(v_guess) = ForwardDiff.derivative(root_func, v_guess)

        z[i] = ZCBYield(solve(root_func, root_func′, q.instrument.coupon_rate),maturities[i])
        # o =Optim.optimize(root_func, [0.0], Optim.LBFGS())
        # z[i] = ZCBYield(Optim.minimizer(o)[1],maturities[i])
    end

    # zero_vec = -log.(clamp.(discount_vec,0.00001,1)) ./ maturities
    return curve(bs,z)
end