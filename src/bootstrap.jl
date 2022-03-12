# bootstrapped class of curve methods

# Forward curves

"""
    ForwardStarting(curve,forwardstart)

Rebase a `curve` so that `discount`/`accumulation`/etc. are re-based so that time zero from the new curves perspective is the given `forwardstart` time.

# Examples

```julia-repl
julia> zero = [5.0, 5.8, 6.4, 6.8] ./ 100
julia> maturity = [0.5, 1.0, 1.5, 2.0]
julia> curve = Yields.Zero(zero, maturity)
julia> fwd = Yields.ForwardStarting(curve, 1.0)

julia> discount(curve,1,2)
0.9275624570410582

julia> discount(fwd,1) # `curve` has effectively been reindexed to `1.0`
0.9275624570410582
```

# Extended Help

While `ForwardStarting` could be nested so that, e.g. the third period's curve is the one-period forward of the second period's curve, it will be more efficient to reuse the initial curve from a runtime and compiler perspective.

`ForwardStarting` is not used to construct a curve based on forward rates. See  [`Forward`](@ref) instead.
"""
struct ForwardStarting{T,U} <: AbstractYield
    curve::U
    forwardstart::T
end

function discount(c::ForwardStarting, to)
    discount(c.curve, c.forwardstart, to + c.forwardstart)
end

function Base.zero(c::ForwardStarting, to,cf::C) where {C<:CompoundingFrequency}
    z = forward(c.curve,c.forwardstart,to+c.forwardstart)
    return convert(cf,z)
end

"""
    zero(curve,time)
    zero(curve,time,CompoundingFrequency)

Return the zero rate for the curve at the given time. If not specified, will use `Periodic(1)` compounding.
"""
Base.zero(c::YC, time) where {YC<:AbstractYield} = zero(c, time, Periodic(1))
function Base.zero(c::YC, time, cf::C) where {YC<:AbstractYield,C<:CompoundingFrequency}
    return convert(cf, Continuous(c.zero(time))) # c.zero is a curve of continuous rates represented as floats. explicitly wrap in continuous before converting
end

"""
    Constant(rate::Real, cf::CompoundingFrequency=Periodic(1))
    Constant(r::Rate)

Construct a yield object where the spot rate is constant for all maturities. If `rate` is not a `Rate` type, will assume `Periodic(1)` for the compounding frequency

# Examples

```julia-repl
julia> y = Yields.Constant(0.05)
julia> discount(y,2)
0.9070294784580498     # 1 / (1.05) ^ 2
```
"""
struct Constant{T} <: AbstractYield
    rate::T
end

function Constant(rate::T, cf::C = Periodic(1)) where {T<:Real,C<:CompoundingFrequency}
    return Constant(Rate(rate, cf))
end

Base.zero(c::Constant, time) = c.rate
Base.zero(c::Constant, time, cf::CompoundingFrequency) = convert(cf, c.rate)
rate(c::Constant) = c.rate
rate(c::Constant, time) = c.rate
discount(r::Constant, time) = discount(r.rate, time)
discount(r::Constant, from, to) = discount(r.rate, to - from)
accumulation(r::Constant, time) = accumulation(r.rate, time)
accumulation(r::Constant, from, to) = accumulation(r.rate, to - from)

"""
    Step(rates,times)

Create a yield curve object where the applicable rate is the effective rate of interest applicable until corresponding time. If `rates` is not a `Vector{Rate}`, will assume `Periodic(1)` type.

# Examples

```julia-repl
julia>y = Yields.Step([0.02,0.05], [1,2])

julia>rate(y,0.5)
0.02

julia>rate(y,1.5)
0.05

julia>rate(y,2.5)
0.05
```
"""
struct Step <: AbstractYield
    rates
    times
end

Step(rates) = Step(rates, collect(1:length(rates)))

function rate(y::Step, time)
    i = findfirst(t -> time <= t, y.times)
    if isnothing(i)
        return y.rates[end]
    else
        return y.rates[i]
    end
end

function discount(y::Step, time)
    v = 1 / (1 + y.rates[1])^min(y.times[1], time)

    if y.times[1] >= time
        return v
    end

    for i = 2:length(y.times)

        if y.times[i] >= time
            # take partial discount and break
            v /= (1 + y.rates[i])^(time - y.times[i-1])
            break
        else
            # take full discount and continue
            v /= (1 + y.rates[i])^(y.times[i] - y.times[i-1])
        end

    end

    return v
end


"""
    Zero(rates, maturities; interpolation=CubicSpline())

Construct a yield curve with given zero-coupon spot `rates` at the given `maturities`. If `rates` is not a `Vector{Rate}`, will assume `Periodic(1)` type.

See [`bootstrap`](@ref) for more on the `interpolation` parameter, which is set to `CubicSpline()` by default.
"""
function Zero(rates::Vector{<:Rate}, maturities; interpolation=CubicSpline())
    # bump to a constant yield if only given one rate
    length(rates) == 1 && return Constant(first(rates))
    return _zero_inner(rates,maturities,interpolation)
end

# zero is different than the other boostrapped curves in that it doesn't actually need to bootstrap 
# because the rate are already zero rates. Instead, we just cut straight to the 
# appropriate interpolation function based on the type dispatch.
function _zero_inner(rates::Vector{<:Rate}, maturities, interp::CubicSpline)
    continuous_zeros = rate.(convert.(Continuous(), rates))
    return YieldCurve(
        rates,
        maturities,
        cubic_interp([0.0; maturities],[first(continuous_zeros); continuous_zeros])
    )
end

function _zero_inner(rates::Vector{<:Rate}, maturities, interp::LinearSpline)
    continuous_zeros = rate.(convert.(Continuous(), rates))
    return YieldCurve(
        rates,
        maturities,
        linear_interp([0.0; maturities],[first(continuous_zeros); continuous_zeros])
    )
end

# fallback for user provied interpolation function
function _zero_inner(rates::Vector{<:Rate}, maturities, interp)
    continuous_zeros = rate.(convert.(Continuous(), rates))
    return YieldCurve(
        rates,
        maturities,
        interp([0.0; maturities],[first(continuous_zeros); continuous_zeros])
    )
end

#fallback if `rates` aren't `Rate`s. Assume `Periodic(1)` per Zero docstring
function Zero(rates, maturities; interpolation=CubicSpline())
    Zero(Periodic.(rates, 1), maturities; interpolation)
end

function Zero(rates; interpolation=CubicSpline())
    # bump to a constant yield if only given one rate
    maturities = collect(1:length(rates))
    return Zero(rates, maturities; interpolation)
end

"""
    Par(rates, maturities; interpolation=CubicSpline())

Construct a curve given a set of bond equivalent yields and the corresponding maturities. Assumes that maturities <= 1 year do not pay coupons and that after one year, pays coupons with frequency equal to the CompoundingFrequency of the corresponding rate.

See [`bootstrap`](@ref) for more on the `interpolation` parameter, which is set to `CubicSpline()` by default.

# Examples

```julia-repl

julia> par = [6.,8.,9.5,10.5,11.0,11.25,11.38,11.44,11.48,11.5] ./ 100
julia> maturities = [t for t in 1:10]
julia> curve = Par(par,maturities);
julia> zero(curve,1)
Rate(0.06000000000000005, Periodic(1))

```
"""
function Par(rates::Vector{<:Rate}, maturities; interpolation=CubicSpline())
    # bump to a constant yield if only given one rate
    if length(rates) == 1
        return Constant(rate[1])
    end
    return YieldCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise semi-annual
        # per Hull 4.7
        bootstrap(rates, maturities, [m <= 1 ? nothing : 1 / r.compounding.frequency for (r, m) in zip(rates, maturities)], interpolation)
    )
end

function Par(rates::Vector{T}, maturities; interpolation=CubicSpline()) where {T<:Real}
    return Par(Rate.(rates), maturities; interpolation)
end


"""
    Forward(rate_vector,maturities)

Takes a vector of 1-period forward rates and constructs a discount curve. If rate_vector is not a vector of `Rates` (ie is just a vector of `Float64` values), then
the assumption is that each value is `Periodic` rate compounded once per period.

# Examples

```julia-repl
julia> Yields.Forward( [0.01,0.02,0.03] );

julia> Yields.Forward( Yields.Continuous.([0.01,0.02,0.03]) );

```
"""
function Forward(rates::Vector{<:Rate}, maturities)
    # convert to zeros and pass to Zero
    disc_v = Vector{Float64}(undef, length(rates))

    v = 1.0

    for i = 1:length(rates)
        Δt = maturities[i] - (i == 1 ? 0 : maturities[i-1])
        v *= discount(rates[i], Δt)
        disc_v[i] = v
    end

    z = (1.0 ./ disc_v) .^ (1 ./ maturities) .- 1 # convert disc_v to zero
    return Zero(z, maturities)
end

# if rates isn't a vector of Rates, then we assume periodic rates compounded once per period.
function Forward(rates, maturities)
    return Forward(Yields.Periodic.(rates, 1), maturities)
end

Forward(rates) = Forward(rates, collect(1:length(rates)))

"""
    Yields.CMT(rates, maturities; interpolation=CubicSpline())

Takes constant maturity (treasury) yields (bond equivalent), and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual.

See [`bootstrap`](@ref) for more on the `interpolation` parameter, which is set to `CubicSpline()` by default.

# Examples

```
# 2021-03-31 rates from Treasury.gov
rates =[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
mats = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]
	
Yields.CMT(rates,mats)
```
"""
function CMT(rates::Vector{T}, maturities; interpolation=CubicSpline()) where {T<:Real}
    rs = map(zip(rates, maturities)) do (r, m)
        if m <= 1
            Rate(r, Periodic(1 / m))
        else
            Rate(r, Periodic(2))
        end
    end

    CMT(rs, maturities;interpolation)
end

function CMT(rates::Vector{<:Rate}, maturities; interpolation=CubicSpline())
    return YieldCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise semi-annual
        # per Hull 4.7
        bootstrap(rates, maturities, [m <= 1 ? nothing : 0.5 for m in maturities], interpolation)
    )
end


"""
    OIS(rates,maturities)
Takes Overnight Index Swap rates, and assumes that instruments <= one year maturity are settled once and other agreements are settled quarterly with a corresponding CompoundingFrequency.

See [`bootstrap`](@ref) for more on the `interpolation` parameter, which is set to `CubicSpline()` by default.

"""
function OIS(rates::Vector{T}, maturities; interpolation=CubicSpline()) where {T<:Real}
    rs = map(zip(rates, maturities)) do (r, m)
        if m <= 1
            Rate(r, Periodic(1 / m))
        else
            Rate(r, Periodic(4))
        end
    end

    return OIS(rs, maturities; interpolation)
end
function OIS(rates::Vector{<:Rate}, maturities ; interpolation=CubicSpline())
    return YieldCurve(
        rates,
        maturities,
        # assume that maturities less than or equal to 12 months are settled once, otherwise quarterly
        # per Hull 4.7
        bootstrap(rates, maturities, [m <= 1 ? nothing : 1 / 4 for m in maturities], interpolation)
    )
end


"""
    par(curve,time;frequency=2)

Calculate the par yield for maturity `time` for the given `curve` and `frequency`. Returns a `Rate` object with periodicity corresponding to the `frequency`. The exception to this is if `time` is less than what the payments allowed by frequency (e.g. a time `0.5` but with frequency `1`) will effectively assume frequency equal to 1 over `time`.
"""
function par(curve, time; frequency=2)
    mat_disc = discount(curve, 0, time)
    coup_times = coupon_times(time,frequency)
    coupon_pv = sum(discount(curve,0,t) for t in coup_times)
    Δt = step(coup_times)
    r = (1-mat_disc) / coupon_pv
    cfs = [t == last(coup_times) ? 1+r : r for t in coup_times]
    cfs = [-1;cfs]
    r = irr_newton(cfs,[0;coup_times])
    frequency_inner = min(1/Δt,max(1 / Δt, frequency))
    r = convert(Periodic(frequency_inner),r)
    return r
end

function coupon_times(time,frequency)
    Δt = min(1 / frequency,time)
    times = time:-Δt:0
    f = last(times)
    f += iszero(f) ? Δt : zero(f)
    l = first(times)
    return f:Δt:l
end

function irr_newton(cashflows, times)
    # use newton's method with hand-coded derivative
    f(r) =  sum(cf * exp(-r*t) for (cf,t) in zip(cashflows,times))
    f′(r) = sum(-t*cf * exp(-r*t) for (cf,t) in zip(cashflows,times) if t > 0)
    # r = Roots.solve(Roots.ZeroProblem((f,f′), 0.0), Roots.Newton())
    r = Roots.newton(x->(f(x),f(x)/f′(x)),0.0)
    return Yields.Periodic(exp(r)-1,1)

end