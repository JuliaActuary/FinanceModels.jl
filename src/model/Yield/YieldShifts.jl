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

# Parallel shift of +100 bp annual effective
base + (z, t) -> Continuous(Periodic(1)(z) + 0.01)

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
struct TenorShift{C <: AbstractYieldModel, F} <: AbstractYieldShift
    base::C
    rule::F
end

function Base.zero(s::TenorShift, t)
    z = convert(Continuous(), Base.zero(s.base, t))   # the rule receives a Continuous rate
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
struct ProjectedShift{C <: AbstractYieldModel, F, T} <: AbstractYieldShift
    base::C
    rule::F
    time::T
end

function Base.zero(s::ProjectedShift, t)
    z = convert(Continuous(), Base.zero(s.base, t))   # the rule receives a Continuous rate
    return convert(Continuous(), s.rule(s.time, z, t)::FinanceCore.Rate)
end
FinanceCore.discount(s::AbstractYieldShift, t) = _discount_from_zero(s, t)
# At t = 0 the forward is the zero rate's limit: the rule applied to the base curve's `zero(base, 0)`,
# with every parameter's derivative. Differentiating L there would need the base's zero rate at a dual
# time of 0, which a curve without its own zero rate refuses. Elsewhere, the derivative of L.
__instantaneous_forward(s::AbstractYieldShift, t) =
    iszero(t) ? FinanceCore.rate(Base.zero(s, t)) : __log_discount_derivative(s, t)
__log_discount(s::AbstractYieldShift, t) = __zero_log_discount(s, t)
__log_native(::AbstractYieldShift) = true

# Deprecated alias for the previous name. Slated for removal one minor release after introduction.
Base.@deprecate_binding TransformedYield TenorShift
