"""
    Yield.FlatForwardAt(forward::FinanceCore.Rate)

Hold the instantaneous forward at the supplied rate beyond the last knot. Pass this object as
`extrapolation` to `ZeroRateCurve`, `Yield.Spline`, `Yield.MonotoneConvex`, `reconstruct`, or a
spline `fit` call.

`forward` must be a `FinanceCore.Rate`, such as `Continuous(0.035)` or `Periodic(0.035, 1)`,
so that its compounding convention is explicit; it is converted to and stored as a
continuously compounded rate in the `forward` field. There is no method for a bare number,
because packages read bare numbers differently (`Yield.Constant(0.035)` is annual effective).

The last-knot discount factor and all interpolation through that knot are preserved.
The forward usually jumps at the boundary. Unlike `:flat_forward`, the supplied
forward is an independent assumption: changing or fitting knot rates keeps it fixed.
The value must be finite; negative rates and automatic differentiation are supported.

```julia
curve = ZeroRateCurve([0.02, 0.03, 0.04], [1.0, 10.0, 30.0];
    extrapolation=Yield.FlatForwardAt(Continuous(0.035)))
forward(curve, 40.0, 60.0)  # Continuous(0.035), up to rounding
```
"""
struct FlatForwardAt{F <: Real}
    forward::F   # continuously compounded
    function FlatForwardAt(forward::FinanceCore.Rate)
        f = float(FinanceCore.rate(convert(Continuous(), forward)))
        isfinite(f) || throw(ArgumentError("FlatForwardAt requires a finite forward rate (got $forward)."))
        return new{typeof(f)}(f)
    end
end

Base.:(==)(a::FlatForwardAt, b::FlatForwardAt) = a.forward == b.forward
Base.isequal(a::FlatForwardAt, b::FlatForwardAt) = isequal(a.forward, b.forward)
Base.hash(p::FlatForwardAt, h::UInt) = hash(p.forward, hash(:FlatForwardAt, h))
Base.show(io::IO, p::FlatForwardAt) = print(io, "Yield.FlatForwardAt(Continuous(", p.forward, "))")

# Policy is public configuration. The tail objects below are derived boundary data,
# rebuilt from the policy and knots on construction, fitting, and Accessors updates. The policy
# is validated where it is used, in `__build_tail`.

# DataInterpolations handles the short end (a flat zero rate before the first knot) and, only
# for `:extension`, the long end; every other long-end policy is a `CurveTail` below.
function __interpolation_extrapolation(method)
    E = DataInterpolations.ExtrapolationType
    return method === :extension ? (; extrapolation_left = E.Constant, extrapolation_right = E.Extension) :
        (; extrapolation_left = E.Constant)
end

# The average (discrete) continuously compounded forward over the last knot interval,
# `(zₙtₙ - zₙ₋₁tₙ₋₁) / (tₙ - tₙ₋₁)`: the `:flat_forward` anchor for DataInterpolations-backed
# curves. With a single knot, the only forward the data imply is the zero rate itself.
function __last_discrete_forward(rates, tenors)
    n = length(tenors)
    n == 1 && return only(rates)
    zₙ, tₙ, zₘ, tₘ = rates[n], tenors[n], rates[n - 1], tenors[n - 1]
    return (zₙ * tₙ - zₘ * tₘ) / (tₙ - tₘ)
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
    if e.linear
        # A zero boundary slope makes the linear tail flat. Return that limit at
        # t = Inf instead of evaluating 0 * Inf. (`iszero` of a ForwardDiff dual also
        # requires zero partials, so a slope that carries derivatives is not dropped.)
        isinf(t) && iszero(e.γ) && return e.α
        return e.α + e.γ * (t - e.last_tenor)
    end
    # Avoid forming fₙ*t, which can overflow at very long finite horizons.
    return e.α + e.β * (e.last_tenor / t)
end

function __tail_forward(e::CurveTail, t)
    e.linear || return e.α
    isinf(t) && iszero(e.γ) && return e.α
    # Do not form 2t: it can overflow while γ*t is still representable.
    return e(t) + e.γ * t
end

# `discount` at t = Inf: the limit of exp(-z(t)·t) = exp(-(α·t + β·tₙ + γ·(t - tₙ)·t)). The sign of
# γ under a `:linear` tail, and then of α, decides it: 0 when positive, Inf when negative. When both
# are zero it is finite, exp(-β·tₙ), which z(Inf)·Inf = 0·Inf would lose as NaN. The signs are
# those of the primal values. A deciding coefficient that is zero there but carries partials has
# no derivative: a bump that moves it can move the limit to 0 or Inf.
function __discount_at_infinity(e::CurveTail)
    d = exp(-__finite_log_discount_at_infinity(e))
    limit = __tail_limit_sign(e)
    return limit > 0 ? zero(d) : limit < 0 ? oftype(d, Inf) : d
end

# The same limit as a cumulative log-discount: +Inf where the discount factor tends to 0, -Inf
# where it tends to Inf, and β·tₙ (0 under a `:linear` tail) where it is finite. It is not formed as
# `-log(__discount_at_infinity(e))`, so each result keeps its own numeric type and partials.
function __log_discount_at_infinity(e::CurveTail)
    L = __finite_log_discount_at_infinity(e)
    limit = __tail_limit_sign(e)
    return limit > 0 ? oftype(L, Inf) : limit < 0 ? oftype(L, -Inf) : L
end

__finite_log_discount_at_infinity(e::CurveTail) = e.linear ? zero(e.β) : e.β * e.last_tenor

# +1 if the discount factor at infinity is 0, -1 if it is Inf, 0 if it is finite.
function __tail_limit_sign(e::CurveTail)
    for (c, name) in (e.linear ? ((e.γ, "slope"), (e.α, "forward")) : ((e.α, "forward"),))
        v = __primal(c)
        v > 0 && return 1
        v < 0 && return -1
        iszero(c) || throw(
            ArgumentError(
                "discount(curve, Inf) has no derivative here: the tail's $name is zero, so a bump " *
                    "that moves it can send the discount factor at infinity to 0 or Inf. " *
                    "Differentiate at a finite time."
            )
        )
    end
    return 0
end

# `forward()` and `slope()` supply the curve's boundary anchor for `:flat_forward` and its
# left-hand zero-rate slope for `:linear`. Keeping them lazy avoids computing quantities
# a policy does not use (e.g. an interpolant derivative for `:flat_forward`).
function __build_tail(method::Symbol, t, z, forward, slope)
    method === :flat_zero && return __curve_tail(t, z, zero(z), zero(z), false)
    method === :linear && return __curve_tail(t, z, zero(z), slope(), true)
    method === :flat_forward || throw(
        ArgumentError(
            "extrapolation must be :flat_forward, :flat_zero, :linear, or Yield.FlatForwardAt(forward), " *
                "or :extension for a DataInterpolations-backed curve (not Spline.MonotoneConvex()); " *
                "got $(repr(method))."
        )
    )
    f = forward()
    return __curve_tail(t, f, z - f, zero(z), false)
end
__build_tail(method::FlatForwardAt, t, z, forward, slope) =
    __curve_tail(t, method.forward, z - method.forward, zero(z), false)

# A single callable adapter combines a DataInterpolations zero-rate interpolant with its
# tail. Every policy, including `:extension`, uses this one type, so a runtime policy
# Symbol does not change the constructed curve's type. With `extend = true`
# (`:extension`) the interpolant also evaluates beyond the last knot and the tail is unused.
struct Extrapolated{I, E}
    interpolant::I
    tail::E
    extend::Bool
end
(e::Extrapolated)(t) = (e.extend || t <= e.tail.last_tenor) ? e.interpolant(t) : e.tail(t)

# The zero-rate function of a `Yield.Spline`: the interpolant through the last knot, then the tail.
function __extrapolate(interpolant, g::KnotGrid, method)
    t, z = last(g.tenors), last(g.rates)
    if method === :extension
        # Same tail type as `:flat_zero`, but never evaluated.
        tail = __curve_tail(t, z, zero(z), zero(z), false)
        return Extrapolated(interpolant, tail, true)
    end
    forward = () -> __last_discrete_forward(g.rates, g.tenors)
    slope = () -> DataInterpolations.derivative(interpolant, t)
    return Extrapolated(interpolant, __build_tail(method, t, z, forward, slope), false)
end
