module Yield
import ..AbstractModel
import ..FinanceCore
import ..Spline as Sp
import ..DataInterpolations
import ..Bond: coupon_times

using ..FinanceCore: Continuous, Periodic, discount, accumulation, forward, pv, AbstractContract

export discount, zero, forward, par, pv, instantaneous_forward

abstract type AbstractYieldModel <: AbstractModel end

# Discount factor derived from the continuous zero rate — the single home for this logic,
# shared by every zero-native curve through the one-line `discount` stubs at each curve
# definition (CompositeYield, ScaledYield, the yield shifts, NelsonSiegel(Svensson),
# CairnsPritchard, MonotoneConvex). Curves with a cheaper direct formula (Constant,
# Yield.Spline) define their own `discount` instead. The affected curves all define a
# finite zero rate at t=0, so delegating to the rate-level `discount` also preserves the
# promoted numeric type carried by the curve parameters.
_discount_from_zero(c, t) = discount(Base.zero(c, t), t)

# Generic callable fallback: `curve(t) ≡ discount(curve, t)`. Covers every
# AbstractYieldModel subtype (Constant, Spline, CompositeYield, ScaledYield,
# TenorShift, ProjectedShift, NelsonSiegel, MonotoneConvex, ZeroRateCurve, …);
# each one routes through its own `discount`, so no per-type callable is needed.
(yc::AbstractYieldModel)(t) = FinanceCore.discount(yc, t)

"""
    Constant(rate)

A yield curve representing a flat term structure. `rate` can be a [`Rate`](@ref) object or a `Real` object.


If [`fit`](@ref FinanceModels.fit-Union{Tuple{F}, Tuple{Any, Any}, Tuple{Any, Any, F}} where F<:FinanceModels.Fit.Loss)ing with the default FinanceModels.jl settings, the solver will attempt to fit a discount rate with the range of: `-1.0 .. 1.0`
"""
struct Constant{R} <: AbstractYieldModel
    rate::R
end

function Constant(rate::R) where {R <: Real}
    return Constant(FinanceCore.Rate(rate))
end

Constant() = Constant(0.0)

FinanceCore.discount(c::Constant, t) = FinanceCore.discount(c.rate, t)

# The continuous zero rate of a flat curve is its (continuous) rate at every tenor,
# including t=0. Defining `zero` directly avoids the generic `-log(discount)/t`
# round-trip — which is `0/0 → NaN` at t=0 — and lets curves composed from a
# `Constant` stay in zero-rate space (see `CompositeYield`/`ScaledYield`).
Base.zero(c::Constant, t) = convert(Continuous(), c.rate)

struct Spline{U} <: AbstractYieldModel
    fn::U # here, fn is a map from time to instantaneous zero rate
end

# `discount`/`zero` below are the interface; the callable `(c::Spline)(t)` comes from
# the generic `AbstractYieldModel` fallback, which routes to this same `discount`.
function FinanceCore.discount(c::Spline, time)
    z = c.fn(time)
    return exp(-z * time)
end

# `c.fn(t)` is the continuous zero rate, so `zero` is a direct read of the
# interpolant. Avoids the generic `-log(discount)/t` round-trip (and its
# `0/0 → NaN` at t=0), which also makes composites/shifts built on a Spline cheap.
Base.zero(c::Spline, time) = Continuous(c.fn(time))

function Spline(b::Sp.BSpline, xs, ys)
    order = min(length(xs) - 1, b.order) # in case the length of xs is less than the spline order
    xs = float.(xs)
    knot_type = if length(xs) < 3
        :Uniform
    else
        :Average
    end

    return Spline(DataInterpolations.BSplineInterpolation(ys, xs, order, :Uniform, knot_type; extrapolation = DataInterpolations.ExtrapolationType.Extension))
end

function Spline(::Sp.PCHIP, xs, ys)
    return Spline(DataInterpolations.PCHIPInterpolation(ys, float.(xs); extrapolation = DataInterpolations.ExtrapolationType.Extension))
end

function Spline(::Sp.Akima, xs, ys)
    return Spline(DataInterpolations.AkimaInterpolation(ys, float.(xs); extrapolation = DataInterpolations.ExtrapolationType.Extension))
end

function Spline(b::Sp.PolynomialSpline, xs, ys)
    order = min(length(xs) - 1, b.order) # in case the length of xs is less than the spline order
    # `cache_parameters = true` precomputes per-segment parameters at construction so that evaluation is
    # read-only and therefore thread-safe — notably for `QuadraticSpline`, which is B-spline-based and
    # otherwise overwrites a shared internal coefficient buffer on every call. It also avoids recomputing
    # parameters on each evaluation. (Curves here are built immutably, so the "do not mutate u/t" caveat
    # of cached parameters does not apply.)
    if order == 1
        return Spline(DataInterpolations.LinearInterpolation(ys, xs; extrapolation = DataInterpolations.ExtrapolationType.Extension, cache_parameters = true))
    elseif order == 2
        return Spline(DataInterpolations.QuadraticSpline(ys, xs; extrapolation = DataInterpolations.ExtrapolationType.Extension, cache_parameters = true))
    else
        return Spline(DataInterpolations.CubicSpline(ys, xs; extrapolation = DataInterpolations.ExtrapolationType.Extension, cache_parameters = true))
    end
end

include("Yield/SmithWilson.jl")
include("Yield/NelsonSiegelSvensson.jl")
include("Yield/CairnsPritchard.jl")
include("Yield/MonotoneConvex.jl")

"""
    build_model(spline, tenors, rates)

Build a yield model from a `SplineCurve` type descriptor, tenor times, and rates.
Returns an `AbstractYieldModel` that supports `discount()` and callable `(t)` syntax.

Used by `ZeroRateCurve` and the AD pathway in ActuaryUtilities to construct the
interpolation model efficiently (once per gradient step rather than per discount call).
"""
build_model(spline::Sp.SplineCurve, tenors, rates) = Yield.Spline(spline, tenors, rates)
build_model(::Sp.MonotoneConvex, tenors, rates) = Yield.MonotoneConvex(collect(rates), collect(float.(tenors)))

include("Yield/ZeroRateCurve.jl")

## Generic and Fallbacks
"""
    discount(yc, to)
    discount(yc, from,to)

The discount factor for the yield curve `yc` for times `from` through `to`.
"""
FinanceCore.discount(yc::T, from, to) where {T <: AbstractYieldModel} = discount(yc, to) / discount(yc, from)

"""
    forward(yc, from, to)˚

The forward `Rate` implied by the yield curve `yc` between times `from` and `to`.
"""
function FinanceCore.forward(yc::T, from, to = from + 1) where {T <: AbstractYieldModel}
    # forward = log(DF(from)/DF(to)) / (to-from) = (z(to)·to − z(from)·from)/(to−from).
    # The `z·t` terms are `−log(DF)`, which is exactly 0 at t=0, so we guard t=0
    # rather than evaluate a (possibly singular) zero rate there. For zero-native
    # curves this is transcendental-free; discount-native curves (SmithWilson, the
    # short-rate models) recover the same value through `zero`'s generic fallback.
    zt(t) = iszero(t) ? zero(float(t)) : FinanceCore.rate(Base.zero(yc, t)) * t
    return Continuous((zt(to) - zt(from)) / (to - from))
end

"""
    par(curve,time;frequency=2)

Calculate the par yield for maturity `time` for the given `curve` and `frequency`. Returns a `Rate` object with periodicity corresponding to the `frequency`.

If `time` is shorter than one regular coupon period (e.g. `time=0.5` with `frequency=1`), the single stub payment implies a compounding frequency of `1/time`: the result is quoted as `Periodic(1/time)` when `1/time` is a (near-)integer, and otherwise an `ArgumentError` is thrown because the implied frequency cannot be represented as a `Periodic` rate.

# Examples

```julia-repl
julia> c = Yield.Constant(0.04);

julia> par(c,4)
Periodic(0.03960780543711406, 2)

julia> par(c,4;frequency=1)
Periodic(0.040000000000000036, 1)

julia> par(c,0.6;frequency=4)
Periodic(0.039413626195875295, 4)

julia> par(c,0.2;frequency=4)
Periodic(0.039374942589460726, 5)

julia> par(c,2.5)
Periodic(0.03960780543711406, 2)
```
"""
function par(curve, time; frequency = 2)
    coup_times = coupon_times(time, frequency)
    mat_disc = discount(curve, time)
    coupon_pv = sum(discount(curve, t) for t in coup_times)
    Δt = step(coup_times)
    r = (1 - mat_disc) / coupon_pv

    # Build cash flows: initial outflow of -1, then coupons r, final coupon+principal 1+r
    # Pre-allocate arrays for better performance
    n = length(coup_times)
    cfs = Vector{typeof(r)}(undef, n + 1)
    times = Vector{typeof(Δt)}(undef, n + 1)

    cfs[1] = -one(r)
    times[1] = zero(Δt)

    @inbounds for i in 1:n
        cfs[i + 1] = i == n ? 1 + r : r
        times[i + 1] = coup_times[i]
    end

    r = FinanceCore.internal_rate_of_return(cfs, times)
    frequency_inner = 1 / Δt  # Simplified from min(1 / Δt, max(1 / Δt, frequency))
    if !isinteger(round(frequency_inner, digits = 8))
        throw(
            ArgumentError(
                "par(curve, $time; frequency=$frequency) implies a coupon period of $Δt and a compounding frequency of 1/Δt = $frequency_inner, which is not an integer and cannot be represented as a `Periodic` rate. Choose a maturity commensurate with the coupon frequency."
            )
        )
    end
    r = convert(Periodic(frequency_inner), r)
    return r
end

"""
    zero(curve,time)

Return the zero rate for the curve at the given time.
"""
function Base.zero(c::YC, time) where {YC <: AbstractYieldModel}
    df = discount(c, time)
    r = -log(df) / time
    return Continuous(r)
end

"""
    accumulation(yc, from, to)

The accumulation factor for the yield curve `yc` for times `from` through `to`.
"""
function FinanceCore.accumulation(yc::AbstractYieldModel, time)
    return 1 ./ discount(yc, time)
end

function FinanceCore.accumulation(yc::AbstractYieldModel, from, to)
    return 1 ./ discount(yc, from, to)
end

## Curve Manipulations
"""
    CompositeYield(curve1,curve2,operation)

Combines two yield curves by applying `operation` to their continuous zero rates.

Given discount factors `DF₁(t)` and `DF₂(t)`, the continuous zero rates are
`z₁ = -log(DF₁)/t` and `z₂ = -log(DF₂)/t`, and the composite discount factor is
`exp(-op(z₁, z₂) * t)`.

For addition (`+`), this gives `DF(t) = DF₁(t) × DF₂(t)` (the no-arbitrage spread relationship).
For subtraction (`-`), this gives `DF(t) = DF₁(t) / DF₂(t)`.

Created via `+` and `-` on `AbstractYieldModel` objects. For scalar multiplication/division,
see [`ScaledYield`](@ref).

Composition is performed in continuous-zero-rate space: a `+`/`-` composite reads each
component's zero rate, combines them, and applies a single `exp` to form the discount
factor — it no longer pays the `log`/`exp` round-trip that earlier versions did. Composing
many curves in a hot loop is still marginally slower than pre-fitting a single combined
curve, but the gap is small.

Curves can be added or subtracted together, but note that this is not always the same thing
as adding or subtracting spreads with rates. If spreads and base rates are expressed as zero
rates, then the curve addition/subtraction has the same effect as re-fitting the yield model
with the rate+spread inputs added together first. Non-zero rates (e.g. par rates) do not have
this same property.

## Examples

```julia
rates = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
spreads = [0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
mats = [1 / 12, 2 / 12, 3 / 12, 6 / 12, 1, 2, 3, 5, 7, 10, 20, 30]


### Zero coupon rates/spreads

q_rf_z = ZCBYield.(rates,mats)
q_s_z = ZCBYield.(spreads,mats)
q_y_z = ZCBYield.(rates + spreads,mats)

c_rf_z = fit(Spline.Linear(),q_rf_z,Fit.Bootstrap())
c_s_z = fit(Spline.Linear(),q_s_z,Fit.Bootstrap())
c_y_z = fit(Spline.Linear(),q_y_z,Fit.Bootstrap())

# adding curves when the spreads were zero spreads works
@test discount(c_rf_z+c_s_z,20) ≈ discount(c_y_z,20)


### Par coupon rates/spreads

q_rf = CMTYield.(rates,mats)
q_s = CMTYield.(spreads,mats)
q_y = CMTYield.(rates + spreads,mats)

c_rf = fit(Spline.Linear(),q_rf,Fit.Bootstrap())
c_s = fit(Spline.Linear(),q_s,Fit.Bootstrap())
c_y = fit(Spline.Linear(),q_y,Fit.Bootstrap())

# adding curves when the spreads were par spreads does not work
@test !(discount(c_rf+c_s,20) ≈ discount(c_y,20))
```
"""
struct CompositeYield{T, U, V} <: AbstractYieldModel
    r1::T
    r2::U
    op::V
end


# Composition happens in continuous-zero-rate space: combine the components' zero
# rates with `op`, then form the discount factor with a single `exp`. This avoids
# the previous round-trip (discount → log → recompose → exp), collapsing the common
# Spline/Constant case from 3 `exp` + 2 `log` to a single `exp`.
function Base.zero(rc::CompositeYield, time)
    z1 = FinanceCore.rate(Base.zero(rc.r1, time))
    z2 = FinanceCore.rate(Base.zero(rc.r2, time))
    return Continuous(rc.op(z1, z2))
end
FinanceCore.discount(rc::CompositeYield, time) = _discount_from_zero(rc, time)

# ─── Yield shifts (TenorShift, ProjectedShift) ────────────────────────────

"""
    AbstractYieldShift <: AbstractYieldModel

Supertype for lazy zero-rate shift models: a curve produced by transforming
a base yield curve's zero rate via a user-supplied rule.

Two concrete subtypes:

- [`TenorShift`](@ref) — shift depends only on the tenor `t`. Use for parallel
  bumps, twists, butterflies — static curve transformations.
- [`ProjectedShift`](@ref) — shift depends on the tenor `t` *and* on a second
  time axis `τ` (projection / as-of / valuation-date time). Use for phase-in
  profiles (BMA SBA, IFRS17 macro scenarios) and any shift whose shape evolves
  across a projection horizon.

Both subtypes implement the standard `AbstractYieldModel` interface (`zero`,
`discount`, `forward`, `pv`).
"""
abstract type AbstractYieldShift <: AbstractYieldModel end

"""
    TenorShift(base, rule)

Lazy zero-rate transformation depending only on the tenor:
`z_new(t) = rule(z_base(t), t)`.

The `rule` function receives the base curve's `Continuous` zero rate and
the tenor, and returns a new rate. It is evaluated on demand — no discretization
or refitting. The base curve's analytic structure is fully preserved.

The rule function must have the signature `(z::Rate, t) -> Rate`. The return
value is type-asserted as `Rate` so that compounding convention is always
carried explicitly; rules returning a plain `Real` will raise a `TypeError`
at call time. Use `z + Continuous(0.01)`, `Periodic(0.04, 2)`, etc., and let
`Rate` arithmetic handle conversion to the curve's continuous representation.

Use this for static shifts — parallel bumps, twists, butterflies — that don't
depend on where you are in projection time. For shifts whose shape evolves
across a projection horizon, see [`ProjectedShift`](@ref).

# Constructing

The most ergonomic way to create a `TenorShift` is via the `+` operator
with an `AbstractYieldModel` and a two-argument function:

```julia
base = Yield.Constant(0.05)

# Parallel shift (+100 bp)
base + (z, t) -> z + Periodic(0.01, 1)

# Tenor-dependent twist (steepener that fades at 30y)
base + (z, t) -> z + Continuous(0.02 * max(0.0, 1.0 - t/30.0))
```

You can also construct directly:

```julia
TenorShift(base, (z, t) -> z + Continuous(0.01))
```

Note: The `+` operator dispatches on `Function`. For callable objects that are
not `Function` subtypes (e.g. custom structs with call syntax), use the direct
constructor: `TenorShift(base, my_callable)`.

`TenorShift` is a post-processing wrapper — it is not a fitting target.
ForwardDiff propagates correctly through the transform for sensitivity analysis,
but the rule function itself should be differentiable if used in an AD context.

`TransformedYield` is retained as a deprecated alias for `TenorShift`.

See also: [`ProjectedShift`](@ref), [`AbstractYieldShift`](@ref),
[`CompositeYield`](@ref), [`ScaledYield`](@ref).
"""
struct TenorShift{C<:AbstractYieldModel,F} <: AbstractYieldShift
    base::C
    rule::F
end

function Base.zero(s::TenorShift, t)
    z = Base.zero(s.base, t)
    return convert(Continuous(), s.rule(z, t)::FinanceCore.Rate)
end

"""
    ProjectedShift(base, rule, time)

Lazy zero-rate transformation that depends on the tenor `t` *and* a second
time axis `τ`:

    z_new(t) = rule(τ, z_base(t), t)

`τ` (stored as the `.time` field) is the **projection time** — the as-of or
valuation-date offset at which this curve is being evaluated. It is distinct
from the tenor `t`, which is time-to-maturity from `τ`.

The rule function must have the signature `(τ, z::Rate, t) -> Rate`. The
return value is type-asserted as `Rate` so that compounding convention is
always carried explicitly; rules returning a plain `Real` will raise a
`TypeError` at call time.

Use this for shifts whose shape evolves across a projection horizon — phase-in
profiles (BMA SBA, IFRS17 macro scenarios), runoff schedules, calendar-rolling
shocks. For static, tenor-only shifts, see [`TenorShift`](@ref).

# Constructing

There is no `+` operator sugar — `ProjectedShift` needs an explicit `τ`, which
fixing at composition time would defeat the purpose of storing the rule as a
year-independent first-class value. Always use the direct constructor:

```julia
base = Yield.Constant(0.05)

# −150 bp parallel shift, phased in linearly over 10 projection years.
phase_in = (τ, z, _) -> z + Continuous(-0.015 * min(τ, 10) / 10)

# Curve as seen at projection year 3 (30% phased in → -45 bp).
c3 = ProjectedShift(base, phase_in, 3.0)

# Curve as seen at projection year 10 (fully phased in → -150 bp).
c10 = ProjectedShift(base, phase_in, 10.0)
```

The intended pattern: store `rule` once as a first-class value, then call
`ProjectedShift(base, rule, τ)` at each `τ` in a projection loop.

See also: [`TenorShift`](@ref), [`AbstractYieldShift`](@ref).
"""
struct ProjectedShift{C<:AbstractYieldModel,F,T} <: AbstractYieldShift
    base::C
    rule::F
    time::T
end

function Base.zero(s::ProjectedShift, t)
    z = Base.zero(s.base, t)
    return convert(Continuous(), s.rule(s.time, z, t)::FinanceCore.Rate)
end
FinanceCore.discount(s::AbstractYieldShift, t) = _discount_from_zero(s, t)

# Deprecated alias for the previous name. Slated for removal one minor release after introduction.
Base.@deprecate_binding TransformedYield TenorShift

"""
    ScaledYield(curve, factor)

A yield model that scales the continuous zero rates of `curve` by a `Real` scalar `factor`.

Created via `curve * scalar` or `curve / scalar`. For example, `curve * 0.79` scales
all continuous zero rates by 0.79, which is useful for after-tax yield calculations.
"""
struct ScaledYield{T<:AbstractYieldModel, S<:Real} <: AbstractYieldModel
    curve::T
    factor::S
end

# Scaling is a multiply in continuous-zero-rate space, so derive `zero` directly; the
# discount factor (a single `exp`, no discount → log round-trip) follows from it.
function Base.zero(sy::ScaledYield, time)
    z = FinanceCore.rate(Base.zero(sy.curve, time))
    return Continuous(z * sy.factor)
end
FinanceCore.discount(sy::ScaledYield, time) = _discount_from_zero(sy, time)

"""
    ForwardStarting(curve,forwardstart)

Rebase a `curve` so that `discount`/`accumulation`/etc. are re-based so that time zero from the new curves perspective is the given `forwardstart` time.

# Examples

```julia-repl
julia> zero = [5.0, 5.8, 6.4, 6.8] ./ 100
julia> maturity = [0.5, 1.0, 1.5, 2.0]
julia> curve = ZeroRateCurve(zero, maturity)
julia> fwd = Yield.ForwardStarting(curve, 1.0)

julia> discount(curve,1,2)
0.9275624570410582

julia> discount(fwd,1) # `curve` has effectively been reindexed to `1.0`
0.9275624570410582
```

# Extended Help

While `ForwardStarting` could be nested so that, e.g. the third period's curve is the one-period forward of the second period's curve, it will be more efficient to reuse the initial curve from a runtime and compiler perspective.

`ForwardStarting` is not used to construct a curve based on forward rates. 
"""
struct ForwardStarting{T, U, V} <: AbstractYieldModel
    curve::U
    forwardstart::T
    discount_to_forwardstart::V
    function ForwardStarting(curve::U, forwardstart::T) where {T, U}
        df = FinanceCore.discount(curve, forwardstart)
        new{T, U, typeof(df)}(curve, forwardstart, df)
    end
end

function FinanceCore.discount(c::ForwardStarting, to)
    return FinanceCore.discount(c.curve, to + c.forwardstart) / c.discount_to_forwardstart
end

"""
    Yield.AbstractYieldModel + Yield.AbstractYieldModel

The addition of two yields will create a `CompositeYield`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the two curves will be added together.
"""
function Base.:+(a::AbstractYieldModel, b::AbstractYieldModel)
    return CompositeYield(a, b, +)
end

function Base.:+(a::Constant, b::Constant)
    z_a = FinanceCore.rate(convert(Continuous(), a.rate))
    z_b = FinanceCore.rate(convert(Continuous(), b.rate))
    return Constant(Continuous(z_a + z_b))
end

function Base.:+(a::T, b::Union{Real,Rate}) where {T <: AbstractYieldModel}
    return a + Constant(b)
end

function Base.:+(a::Union{Real,Rate}, b::T) where {T <: AbstractYieldModel}
    return Constant(a) + b
end

function Base.:+(a::AbstractYieldModel, f::Function)
    return TenorShift(a, f)
end

function Base.:+(f::Function, a::AbstractYieldModel)
    return TenorShift(a, f)
end

"""
    curve * scalar
    scalar * curve

Scale the continuous zero rates of `curve` by a `Real` scalar. Returns a [`ScaledYield`](@ref).

This is useful for after-tax yield calculations. For example, `curve * 0.79` produces a
curve whose continuous zero rate at every point is 79% of the original.

# Examples

```julia-repl
julia> m = Yield.Constant(Continuous(0.05)) * 0.79;

julia> discount(m, 1) ≈ exp(-0.05 * 0.79)
true
```
"""
function Base.:*(a::AbstractYieldModel, b::Real)
    return ScaledYield(a, b)
end

function Base.:*(a::Real, b::AbstractYieldModel)
    return ScaledYield(b, a)
end

"""
    Yield.AbstractYieldModel - Yield.AbstractYieldModel

The subtraction of two yields will create a `CompositeYield`. For `rate`, `discount`, and `accumulation` purposes the spot rates of the second curves will be subtracted from the first.
"""
function Base.:-(a::AbstractYieldModel, b::AbstractYieldModel)
    return CompositeYield(a, b, -)
end

function Base.:-(a::Constant, b::Constant)
    z_a = FinanceCore.rate(convert(Continuous(), a.rate))
    z_b = FinanceCore.rate(convert(Continuous(), b.rate))
    return Constant(Continuous(z_a - z_b))
end

function Base.:-(a::T, b::Union{Real,Rate}) where {T <: AbstractYieldModel}
    return a - Constant(b)
end

function Base.:-(a::Union{Real,Rate}, b::T) where {T <: AbstractYieldModel}
    return Constant(a) - b
end

"""
    curve / scalar

Scale the continuous zero rates of `curve` by `1/scalar`. Returns a [`ScaledYield`](@ref).

This is useful for grossing-up a yield to a pre-tax equivalent.

# Examples

```julia-repl
julia> m = Yield.Constant(Continuous(0.05)) / 0.79;

julia> discount(m, 1) ≈ exp(-0.05 / 0.79)
true
```
"""
function Base.:/(a::AbstractYieldModel, b::Real)
    return ScaledYield(a, inv(b))
end


end
