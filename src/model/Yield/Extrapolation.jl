"""
    Yield.FlatForwardAt(forward)

Hold the instantaneous forward at the supplied rate
beyond the last knot. Pass this object as `extrapolation` to `Yield.Spline`,
`Yield.MonotoneConvex`, `ZeroRateCurve`, `Yield.build_model`, or a spline `fit` call.

The last-knot discount factor and all interpolation through that knot are preserved.
The forward usually jumps at the boundary. Unlike `:flat_forward`, the supplied
forward is an independent assumption: changing or fitting knot rates keeps it fixed.
The value must be finite; negative rates and automatic differentiation are supported.
Bare numbers are continuously compounded, unlike `Yield.Constant`'s annual-effective
convention. A `FinanceCore.Rate` is converted from its stated compounding convention.

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

FlatForwardAt(forward::FinanceCore.Rate) =
    FlatForwardAt(FinanceCore.rate(convert(Continuous(), forward)))

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

# All policies share a coefficient layout so a runtime Symbol cannot change the
# constructed curve's type. The coefficients represent
# z(t) = α + β*tₙ/t + γ*(t-tₙ), f(t) = α + γ*(2t-tₙ).
# `linear` selects which nonconstant term can be present. Keeping that selection
# as a value avoids evaluating structurally absent terms (notably 0*Inf) and does
# not mistake a zero-valued Dual coefficient for an identically zero term.
struct CurveTail{T, C}
    last_tenor::T
    α::C
    β::C
    γ::C
    linear::Bool
end

# Include the tenor's numeric type in coefficient promotion even for flat tails;
# otherwise mixed-precision grids could still select different C types by policy.
function __curve_tail(t, α, β, γ, linear)
    a, b, c, _ = promote(α, β, γ, zero(t))
    return CurveTail{typeof(t), typeof(a)}(t, a, b, c, linear)
end

function (e::CurveTail)(t)
    e.linear && return e.α + e.γ * (t - e.last_tenor)
    # Avoid forming fₙ*t, which can overflow at very long finite horizons.
    return e.α + e.β * (e.last_tenor / t)
end

function __tail_forward(e::CurveTail, t)
    # Do not form 2t: it can overflow while γ*t is still representable.
    return e.linear ? e(t) + e.γ * t : e.α
end

# `boundary()` supplies the interpolator's own endpoint forward and zero slope.
# Keeping it lazy avoids computing derivatives for flat-zero or supplied-forward
# policies, and lets MonotoneConvex use its native forward without a round-trip.
function __build_tail(method::Symbol, t, z, boundary)
    method === :flat_zero && return __curve_tail(t, z, zero(z), zero(z), false)
    method in (:flat_forward, :linear) || throw(
        ArgumentError(
            "cannot build a financial tail for extrapolation=$(repr(method)); " *
                ":extension must delegate to the interpolant."
        )
    )
    endpoint = boundary()
    method === :linear && return __curve_tail(t, z, zero(z), endpoint.slope, true)
    return __curve_tail(t, endpoint.forward, z - endpoint.forward, zero(z), false)
end
__build_tail(method::FlatForwardAt, t, z, boundary) =
    __curve_tail(t, method.forward, z - method.forward, zero(z), false)

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
