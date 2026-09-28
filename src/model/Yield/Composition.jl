## Curve Manipulations
"""
    CompositeYield(curve1, curve2, op)

Combines two yield curves by adding (`op = +`) or subtracting (`op = -`) their continuous zero
rates. Created via `+` and `-` on `AbstractYieldModel` objects; for scalar multiplication or
division, see [`ScaledYield`](@ref).

Given discount factors `DF₁(t)` and `DF₂(t)` with continuous zero rates `z₁` and `z₂`, the
composite discount factor is `exp(-(z₁ ± z₂) t)`:

- `+` gives `DF(t) = DF₁(t) × DF₂(t)`, the product of the two discount factors (for example a
  base curve and a spread);
- `-` gives `DF(t) = DF₁(t) / DF₂(t)`, their quotient.

These are the operations that compose the curves' factors: every interval factor of the result
is the product (or quotient) of the components' interval factors, so, for example,
`ForwardStarting(a + b, τ)` prices like `ForwardStarting(a, τ) + ForwardStarting(b, τ)`. Other
operations on the zero rates (`max`, `*`, …) are not accepted: they would build a new curve from
zero rates measured from time 0, not a composition of the two curves. For a pointwise
transformation of a curve's zero rates, use a [`TenorShift`](@ref), `curve + ((z, t) -> ...)`,
or define a curve type with its own `zero`.

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
struct CompositeYield{T, U, V <: Union{typeof(+), typeof(-)}} <: AbstractYieldModel
    r1::T
    r2::U
    op::V
end


# Composition happens in log-discount space: `+` adds the components' cumulative log-discounts
# (multiplying their discount factors) and `-` subtracts them (dividing), then a single `exp`
# forms the discount factor. The zero rate is the same combination of the components' zero rates.
function Base.zero(rc::CompositeYield, time)
    z1 = FinanceCore.rate(Base.zero(rc.r1, time))
    z2 = FinanceCore.rate(Base.zero(rc.r2, time))
    return Continuous(rc.op(z1, z2))
end
# At t = Inf the components' tails are combined before the limit is taken. The components' L is
# inlined at these calls: a composite inlined into a loop (`pv`'s `map`) otherwise copies each
# component onto the stack for every call. Inlining a Spline's L into every caller instead makes
# the generic interval too large to inline into bond pricing.
function __log_discount(rc::CompositeYield, time)
    __at_infinity(time) && return __log_discount_at_infinity(rc)
    return rc.op((@inline __log_discount(rc.r1, time)), (@inline __log_discount(rc.r2, time)))
end
function FinanceCore.discount(rc::CompositeYield, time)
    __at_infinity(time) && return __discount_at_infinity(rc)
    return exp(-__log_discount(rc, time))
end
__log_tail(rc::CompositeYield) = __combine_tails(rc.op, __log_tail(rc.r1), __log_tail(rc.r2))

# A wrapper's interval combines its components' intervals, so each keeps its own form: a log-native
# curve's difference of log-discounts (L(to) alone from 0), Smith–Wilson's signed ratio, a rebased
# curve's base interval, the ratio D(to)/D(from) of a curve that defines only `discount` (whose D(0)
# need not be 1). Those stay finite where both of a wrapper's discount factors underflow. At an
# infinite endpoint the components' limits alone lose information (flat forwards of 4% and −2% have
# L = Inf and −Inf), so the wrapper's own L there comes from its combined tail.
__at_infinity_either(from, to) = isinf(from) || isinf(to)
function __log_interval(rc::CompositeYield, from, to)
    __at_infinity_either(from, to) && return __log_discount(rc, to) - __log_discount(rc, from)
    return rc.op(__log_interval(rc.r1, from, to), __log_interval(rc.r2, from, to))
end

"""
    ScaledYield(curve, factor)

A yield model that scales the continuous zero rates of `curve` by a `Real` scalar `factor`.

Created via `curve * scalar` or `curve / scalar`. For example, `curve * 0.79` scales
all continuous zero rates by 0.79, which is useful for after-tax yield calculations.
"""
struct ScaledYield{T <: AbstractYieldModel, S <: Real} <: AbstractYieldModel
    curve::T
    factor::S
end

# Scaling multiplies the zero rate, and so the cumulative log-discount, by `factor` (it raises
# the discount factor to that power); the discount factor is a single `exp` of it.
function Base.zero(sy::ScaledYield, time)
    z = FinanceCore.rate(Base.zero(sy.curve, time))
    return Continuous(z * sy.factor)
end
# The curve's L is inlined at this call, as for `CompositeYield`.
function __log_discount(sy::ScaledYield, time)
    __at_infinity(time) && return __log_discount_at_infinity(sy)
    return sy.factor * @inline(__log_discount(sy.curve, time))
end
function FinanceCore.discount(sy::ScaledYield, time)
    __at_infinity(time) && return __discount_at_infinity(sy)
    return exp(-__log_discount(sy, time))
end
__log_tail(sy::ScaledYield) = __scale_tail(sy.factor, __log_tail(sy.curve))
function __log_interval(sy::ScaledYield, from, to)
    __at_infinity_either(from, to) && return __log_discount(sy, to) - __log_discount(sy, from)
    return sy.factor * __log_interval(sy.curve, from, to)
end
function FinanceCore.discount(w::Union{CompositeYield, ScaledYield}, from, to)
    d = exp(-__log_interval(w, from, to))
    return from == to ? one(d) : d
end

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
struct ForwardStarting{T, U} <: AbstractYieldModel
    curve::U
    forwardstart::T
end

# The rebased curve's factors are the base curve's interval factors from `forwardstart`, so the
# base curve's own interval method applies (finite in the far tail for the built-in curves).
FinanceCore.discount(c::ForwardStarting, to) = FinanceCore.discount(c.curve, c.forwardstart, to + c.forwardstart)
FinanceCore.discount(c::ForwardStarting, from, to) =
    FinanceCore.discount(c.curve, from + c.forwardstart, to + c.forwardstart)
# L of the rebased curve is the base curve's log-discount over [forwardstart, forwardstart + t], so it
# is real wherever the rebased discount factor is positive, even where the base curve's own discount
# factors are negative.
function __log_discount(c::ForwardStarting, t)
    __at_infinity(t) && return __log_discount_at_infinity(c)
    return __log_interval(c.curve, c.forwardstart, t + c.forwardstart)
end
__log_interval(c::ForwardStarting, from, to) = __log_interval(c.curve, from + c.forwardstart, to + c.forwardstart)
__log_tail(c::ForwardStarting) = __shift_tail(__log_tail(c.curve), c.forwardstart, __log_discount(c.curve, c.forwardstart))

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

function Base.:+(a::T, b::Union{Real, Rate}) where {T <: AbstractYieldModel}
    return a + Constant(b)
end

function Base.:+(a::Union{Real, Rate}, b::T) where {T <: AbstractYieldModel}
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

function Base.:-(a::T, b::Union{Real, Rate}) where {T <: AbstractYieldModel}
    return a - Constant(b)
end

function Base.:-(a::Union{Real, Rate}, b::T) where {T <: AbstractYieldModel}
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
