# make interest curve broadcastable so that you can broadcast over multiple`time`s in `interest_rate`
Base.Broadcast.broadcastable(ic::T) where {T<:AbstractYieldCurve} = Ref(ic)

function coupon_times(maturity,frequency)
    Δt = min(1 / frequency,maturity)
    times = maturity:-Δt:0
    f = last(times)
    f += iszero(f) ? Δt : zero(f)
    l = first(times)
    return f:Δt:l
end

# internal function (will be used in EconomicScenarioGenerators)
# defines the rate output given just the type of curve
__ratetype(::Type{Rate{T,U}}) where {T,U<:Periodic} = Yields.Rate{T, Periodic}
__ratetype(::Type{Rate{T,U}}) where {T,U<:Continuous} = Yields.Rate{T, Continuous}
__ratetype(curve::T) where {T<:AbstractYieldCurve} = __ratetype(typeof(curve))
__ratetype(::Type{T}) where {T<:AbstractYieldCurve} = Yields.Rate{Float64, typeof(DEFAULT_COMPOUNDING)}

# https://github.com/dpsanders/hands_on_julia/blob/master/during_sessions/Fractale%20de%20Newton.ipynb
newton(f, f′, x) = x - f(x) / f′(x)
function solve(g, g′, x0, max_iterations = 100)
    x = x0

    tolerance = 2 * eps(x0)
    iteration = 0

    while (abs(g(x) - 0) > tolerance && iteration < max_iterations)
        x = newton(g, g′, x)
        iteration += 1
    end

    return x
end


# convert to a given rate type if not already a rate
__coerce_rate(x::T,cf) where {T<:Rate} = return x
__coerce_rate(x,cf) = return Rate(x,cf)


abstract type InterpolationKind end

QuadraticSpline() = cubic_interp
LinearSpline() = linear_interp


"""
    _bootstrap(rates, maturities, settlement_frequency, interpolation_function)

Bootstrap the rates with the given maturities, treating the rates according to the periodic frequencies in settlement_frequency. 

`interpolator` is any function that will take two vectors of inputs and output points and return a function that will estimate an output given a scalar input. That is
`interpolator` should be: `interpolator(xs, ys) -> f(x)` where `f(x)` is the interpolated value of `y` at `x`. 

Built in `interpolator`s in Yields are: 
- `QuadraticSpline()`: Quadratic spline interpolation.
- `LinearSpline()`: Linear spline interpolation.

The default is `QuadraticSpline()`.
"""
function _bootstrap(rates, maturities, settlement_frequency, interpolation_function)
    discount_vec = zeros(length(rates)) # construct a placeholder discount vector matching maturities
    # we have to take the first rate as the starting point
    discount_vec[1] = discount(Constant(rates[1]), maturities[1])

    for t = 2:length(maturities)
        if isnothing(settlement_frequency[t])
            # no settlement before maturity
            discount_vec[t] = discount(Constant(rates[t]), maturities[t])
        else
            # need to account for the interim cashflows settled
            times = settlement_frequency[t]:settlement_frequency[t]:maturities[t]
            cfs = [rate(rates[t]) * settlement_frequency[t] for s in times]
            cfs[end] += 1

            function pv(v_guess)
                v = interpolation_function([[0.0]; maturities[1:t]], vcat(1.0, discount_vec[1:t-1], v_guess...))
                return sum(v.(times) .* cfs)
            end
            target_pv = sum(map(t2 -> discount(Constant(rates[t]), t2), times) .* cfs)
            root_func(v_guess) = pv(v_guess) - target_pv
            root_func′(v_guess) = ForwardDiff.derivative(root_func, v_guess)
            discount_vec[t] = solve(root_func, root_func′, rate(rates[t]))
        end

    end
    zero_vec = -log.(clamp.(discount_vec,0.00001,1)) ./ maturities
    return interpolation_function([0.0; maturities], [first(zero_vec); zero_vec])
end

# the ad-hoc approach to extrapoliatons is based on suggestion by author of 
# BSplineKit at https://github.com/jipolanco/BSplineKit.jl/issues/19
# this should not be exposed directly to user
struct _Extrap{I,L,R}
	int::I # the BSplineKit interpolation
	left::L # a tuple of (boundary, extrapolation function)
	right::R # a tuple of (boundary, extrapolation function)
end

function _wrap_spline(itp)

	S = BSplineKit.spline(itp)  # spline passing through data points
	B = BSplineKit.basis(S)     # B-spline basis
	
	a, b = BSplineKit.boundaries(B)  # left and right boundaries
	
	# For now, we construct the full spline S′(x).
	# There are faster ways of doing this that should be implemented...
	S′ = diff(S, BSplineKit.Derivative(1))
	
	return _Extrap(itp,
		(boundary = a, func = x->S(a) + S′(a)*(x-a)),
		(boundary = b, func = x->S(b) + S′(b)*(x-b)),
		
	)
end

function _interp(e::_Extrap,x)
	if x <= e.left.boundary
		return e.left.func(x)
	elseif x >= e.right.boundary
		return e.right.func(x)
	else
		return e.int(x)
	end
end

function linear_interp(xs, ys)
    int = BSplineKit.interpolate(xs, ys, BSplineKit.BSplineOrder(2))
    e = _wrap_spline(int)
    return x -> _interp(e, x)
end

function cubic_interp(xs, ys)
    order = min(length(xs),3) # in case the length of xs is less than the spline order
    int = BSplineKit.interpolate(xs, ys, BSplineKit.BSplineOrder(order))
    e = _wrap_spline(int)
    return x -> _interp(e, x)
end

# used to display simple type name in show method
# https://stackoverflow.com/questions/70043313/get-simple-name-of-type-in-julia?noredirect=1#comment123823820_70043313
name(::Type{T}) where {T} = (isempty(T.parameters) ? T : T.name.wrapper)

function Base.show(io::IO, curve::T) where {T<:AbstractYieldCurve}
    println() # blank line for padding
    r = zero(curve, 1)
    ylabel = isa(r.compounding, Continuous) ? "Continuous" : "Periodic($(r.compounding.frequency))"
    kind = name(typeof(curve))
    l = lineplot(
        t -> rate(zero(curve, t)),
        0.0, #from 
        30.0,  # to
        xlabel = "time",
        ylabel = ylabel,
        compact = true,
        name = "Zero rates",
        width = 60,
        title = "Yield Curve ($kind)"
    )
    show(io, l)
end