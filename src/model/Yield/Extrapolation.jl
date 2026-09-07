"""
    Yield.FlatForwardAt(forward)

Hold the instantaneous forward at the supplied, continuously-compounded scalar rate
beyond the last knot. Pass this object as `extrapolation` to `Yield.Spline`,
`Yield.MonotoneConvex`, `ZeroRateCurve`, `Yield.build_model`, or a spline `fit` call.

The last-knot discount factor and all interpolation through that knot are preserved.
The forward usually jumps at the boundary. Unlike `:flat_forward`, the supplied
forward is an independent assumption: changing or fitting knot rates keeps it fixed.
The value must be finite; negative rates and automatic differentiation are supported.

```julia
curve = ZeroRateCurve([0.02, 0.03, 0.04], [1.0, 10.0, 30.0];
    extrapolation=Yield.FlatForwardAt(0.035))
forward(curve, 40.0, 60.0)  # Continuous(0.035), up to rounding
```
"""
struct FlatForwardAt{F <: Real}
    forward::F
    function FlatForwardAt(forward::Real)
        isfinite(forward) || throw(ArgumentError("FlatForwardAt requires a finite forward rate."))
        f = float(forward)
        return new{typeof(f)}(f)
    end
end

Base.:(==)(a::FlatForwardAt, b::FlatForwardAt) = a.forward == b.forward
Base.isequal(a::FlatForwardAt, b::FlatForwardAt) = isequal(a.forward, b.forward)
Base.hash(p::FlatForwardAt, h::UInt) = hash(p.forward, hash(:FlatForwardAt, h))

# Policy is public configuration. The tail objects below are derived boundary data,
# rebuilt from the policy and knots on construction, fitting, and Accessors updates.
const __EXTRAPOLATION_METHODS = (:flat_forward, :flat_zero, :linear, :extension)
function __extrapolation_method(method)
    method isa FlatForwardAt && return method
    method isa Symbol && method in __EXTRAPOLATION_METHODS || throw(
        ArgumentError(
            "extrapolation must be one of $(join(__EXTRAPOLATION_METHODS, ", ")) or " *
                "Yield.FlatForwardAt(forward); got $(repr(method))."
        )
    )
    return method
end

function __monotone_extrapolation_method(method)
    method = __extrapolation_method(method)
    method === :extension && throw(
        ArgumentError(
            "extrapolation=:extension is only available for DataInterpolations-backed curves; " *
                "use :flat_forward, :flat_zero, :linear, or Yield.FlatForwardAt(forward) with MonotoneConvex."
        )
    )
    return method
end

function __interpolation_extrapolation(method)
    method = __extrapolation_method(method)
    extension = DataInterpolations.ExtrapolationType.Extension
    return method === :extension ? (; extrapolation = extension) : (; extrapolation_left = extension)
end

struct FlatForwardTail{T, Z, F}
    last_tenor::T
    last_zero::Z
    last_forward::F
end

function (e::FlatForwardTail)(t)
    # Avoid forming fₙ*t, which can overflow at very long finite horizons.
    return e.last_forward + (e.last_zero - e.last_forward) * (e.last_tenor / t)
end
__tail_forward(e::FlatForwardTail, t) = e.last_forward

struct FlatZeroTail{T, Z}
    last_tenor::T
    last_zero::Z
end
(e::FlatZeroTail)(t) = e.last_zero
__tail_forward(e::FlatZeroTail, t) = e.last_zero

struct LinearZeroTail{T, Z, D}
    last_tenor::T
    last_zero::Z
    last_derivative::D
end
(e::LinearZeroTail)(t) = e.last_zero + e.last_derivative * (t - e.last_tenor)
__tail_forward(e::LinearZeroTail, t) = e(t) + t * e.last_derivative

# `boundary()` supplies the interpolator's own endpoint forward and zero slope.
# Keeping it lazy avoids computing derivatives for flat-zero or supplied-forward
# policies, and lets MonotoneConvex use its native forward without a round-trip.
function __build_tail(method::Symbol, t, z, boundary)
    method === :flat_zero && return FlatZeroTail(t, z)
    endpoint = boundary()
    method === :linear && return LinearZeroTail(t, z, endpoint.slope)
    return FlatForwardTail(t, z, endpoint.forward)
end
__build_tail(method::FlatForwardAt, t, z, boundary) = FlatForwardTail(t, z, method.forward)

# A single callable adapter combines a DataInterpolations zero-rate interpolant
# with its tail. Native curves reuse the same tails in their own evaluation methods.
struct Extrapolated{I, E}
    interpolant::I
    tail::E
end
(e::Extrapolated)(t) = t <= e.tail.last_tenor ? e.interpolant(t) : e.tail(t)

function __extrapolate(interpolant, g::KnotGrid, method)
    method = __extrapolation_method(method)
    method === :extension && return Spline(interpolant)
    t, z = last(g.tenors), last(g.rates)
    boundary = () -> begin
        slope = DataInterpolations.derivative(interpolant, t)
        (forward = z + t * slope, slope = slope)
    end
    return Spline(Extrapolated(interpolant, __build_tail(method, t, z, boundary)))
end
