# ALM projections: the balance movements of assets with principal, under projection assumptions.

"""
    AssetProjection(; pricing = nothing, default_times = nothing, rng = nothing)

A [`ProjectionKind`](@ref) for asset-liability projections of assets with principal: fixed bonds,
amortizing loans, pools, their [`Credit.Defaultable`](@ref FinanceModels.Credit.Defaultable)
versions and their [`Bond.Callable`](@ref FinanceModels.Bond.Callable) and
[`Bond.Putable`](@ref FinanceModels.Bond.Putable) versions. `collect(Projection(asset, ctx,
AssetProjection(...)))` gives one row per payment date, a named tuple

    (; time, opening, interest, principal, prepayment, defaulted, recovery, redeemed, redemption, closing)

per unit of the asset's initial principal:

- `opening`, `closing`: the principal outstanding after the previous payment and after this one;
- `interest`, `principal`, `prepayment`: paid by the surviving asset at this date;
- `defaulted`, `recovery`: the principal outstanding at defaults since the previous payment, and
  the recovery paid on it (the default loss is `defaulted - recovery`);
- `redeemed`, `redemption`: the principal called or put at this date, and the cash paid for it.

Each row reconciles, `opening - principal - prepayment - defaulted - redeemed == closing` to
rounding, and its `closing` is the next row's `opening`. Contractual payments stop after a default,
an exercise or the last payment.

`ctx` holds the projection's assumptions, which need not be the pricing ones: its store holds the
credit and prepayment assumptions the asset reads by key, and its valuation model is the rate path
when the asset has an exercise right (a `RatePath` from [`simulate`](@ref), whose grid must include
the exercise dates).

**Defaults.** Without `default_times` and `rng`, defaults are expected: survival weights each
payment, and `defaulted` and `recovery` are the expected amounts over the interval, reported at its
end. With them, an issuer defaults at a sampled time: `default_times[key]` if given (`Inf` for no
default), otherwise drawn from `rng` by inverting the cumulative hazard and stored in
`default_times`, so that every asset of the same issuer in the projection, and in later projections
with the same `AssetProjection`, defaults at the same time. A default row at that time pays the
recovery and closes the asset. Pools always use expected defaults: a pool's loans don't default
together.

**Exercise.** An exercisable asset's decisions come from [`exercise_policy`](@ref)`(pricing,
asset)`, computed once per asset and kept in the `AssetProjection`: at each exercise date the asset
is exercised by the policy at the short rate of the path then, `short_rate(path, t)`. A decision
reads only the current rate. It is exact when the path comes from the pricing Hull–White model, and
otherwise applies the pricing model's decisions to the observed rate.

An `AssetProjection` holds sampled default times and policies as they are first needed, so give
each thread its own.

```julia
pool = Credit.Defaultable(Bond.Pool(Bond.LevelPayment(0.06, 12, 10), "speed"), "ACME")
store = Dict("speed" => Prepayment.CPR(0.08), "ACME" => Credit.Assumptions(Credit.ConstantHazard(0.01); recovery = 0.5))
rows = collect(Projection(pool, Models(curve, store), AssetProjection()))

# a callable bond along simulated Hull–White paths, with the pricing lattice's decisions
paths = simulate(hw; n_scenarios = 100, timestep = 1 / 12, horizon = 10.0)
kind = AssetProjection(; pricing = Models(ShortRate.HullWhiteTree(hw)))
[collect(Projection(callable, Models(path), kind)) for path in paths]
```
"""
struct AssetProjection{P, D, R} <: ProjectionKind
    pricing::P
    default_times::D
    rng::R
    policies::IdDict{Any, ExercisePolicy}
end
function AssetProjection(; pricing = nothing, default_times = nothing, rng = nothing)
    times = (default_times === nothing && rng !== nothing) ? Dict{Any, Float64}() : default_times
    return AssetProjection(pricing, times, rng, IdDict{Any, ExercisePolicy}())
end

# The exercise policy of `asset`, computed from the pricing context the first time it is needed.
function __projection_policy(kind::AssetProjection, asset)
    kind.pricing === nothing && throw(
        ArgumentError(
            "projecting an exercisable asset needs its exercise decisions: pass `pricing`, a valuation context " *
                "whose model is a ShortRate.HullWhiteTree, to AssetProjection"
        )
    )
    return get!(() -> exercise_policy(kind.pricing, asset), kind.policies, asset)
end

# The sampled default time of the issuer `key`, drawn once and shared.
function __sampled_default(kind::AssetProjection, key, hazard)
    haskey(kind.default_times, key) && return kind.default_times[key]
    kind.rng === nothing && throw(KeyError(key))
    τ = Credit.__default_time(hazard, randexp(kind.rng))
    kind.default_times[key] = τ
    return τ
end

@inline function __asset_projection(rf, val, p)
    spec = __asset_spec(p.contract)
    kind = p.kind
    rows = __asset_rows(spec.base, p.model)
    credit = spec.key === nothing ? nothing : p.model[spec.key]
    sampled = credit !== nothing && kind.default_times !== nothing && !(spec.base isa Bond.Pool)
    τ = sampled ? __sampled_default(kind, spec.key, credit.hazard) : nothing
    policy = spec.exercise === nothing ? nothing : __projection_policy(kind, p.contract)
    exercise_at = spec.exercise === nothing ? Dict{Int, Int}() : Dict(row => e for (e, row) in enumerate(__exercise_rows(spec.exercise, rows)))
    R = credit === nothing ? false : credit.recovery
    T = promote_type(
        typeof(first(rows).payment), typeof(first(rows).continuing), typeof(first(rows).opening), typeof(R),
        credit === nothing ? Bool : typeof(Credit.survival(credit.hazard, 0.0, 1.0))
    )
    alive = one(T)   # the share of the asset alive just after the previous payment
    t_prev = zero(first(rows).time)
    z = zero(T)
    for (i, r) in enumerate(rows)
        opening = alive * r.opening
        if sampled && τ <= r.time
            val = @next(rf, val, (; time = τ, opening, interest = z, principal = z, prepayment = z, defaulted = opening, recovery = R * opening, redeemed = z, redemption = z, closing = z))
            return complete(rf, val)
        end
        q = (credit === nothing || sampled) ? z : Credit.default_probability(credit.hazard, t_prev, r.time)
        defaulted = opening * q
        alive *= 1 - q
        interest, principal, prepayment = alive * r.interest, alive * r.principal, alive * r.prepayment
        alive *= r.continuing
        closing = alive * r.closing
        e = get(exercise_at, i, nothing)
        if e !== nothing && policy(r.time, short_rate(valuation_model(p.model), r.time))
            redemption = spec.exercise.prices[e] * closing
            val = @next(rf, val, (; time = r.time, opening, interest, principal, prepayment, defaulted, recovery = R * defaulted, redeemed = closing, redemption, closing = z))
            return complete(rf, val)
        end
        val = @next(rf, val, (; time = r.time, opening, interest, principal, prepayment, defaulted, recovery = R * defaulted, redeemed = z, redemption = z, closing))
        t_prev = r.time
    end
    return complete(rf, val)
end

for C in (Bond.Fixed, Bond.Amortizing, Bond.Pool, Credit.Defaultable, Bond.Callable, Bond.Putable)
    @eval @inline Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: $C, M, K <: AssetProjection} = __asset_projection(rf, val, p)
end
