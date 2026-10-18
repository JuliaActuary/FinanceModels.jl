# Exercisable assets on a Hull–White trinomial lattice, and the exercise policy for projections.

# ── The asset behind its wrappers ──────────────────────────────────────────────────────────

# An exercisable or defaultable asset unwrapped: the asset with principal (`base`), the credit key
# (`nothing` without default risk) and the exercise right (`nothing`, or the `Bond.Callable` or
# `Bond.Putable`). The wrappers compose in any order; a feature given twice throws.
function __asset_spec(c)
    base, key, exercise = c, nothing, nothing
    while base isa Union{Credit.Defaultable, Bond.Callable, Bond.Putable}
        if base isa Credit.Defaultable
            key === nothing || throw(ArgumentError("the asset is Credit.Defaultable twice (keys $(repr(key)) and $(repr(base.key)))"))
            key, base = base.key, base.asset
        else
            exercise === nothing || throw(ArgumentError("the asset has two exercise schedules; one call or put schedule per asset is supported"))
            exercise, base = base, base.asset
        end
    end
    return (; base, key, exercise)
end

# The row of each exercise time: a payment date of the asset before its maturity.
function __exercise_rows(exercise, rows)
    return map(exercise.times) do t
        i = findfirst(r -> abs(__primal(r.time) - __primal(t)) <= 1.0e-9 * max(1, abs(__primal(t))), rows)
        i === nothing && throw(ArgumentError("exercise time $t is not a payment date of the asset; exercise is on payment dates"))
        i == length(rows) && throw(ArgumentError("exercise time $t is the asset's maturity; exercise at maturity is not supported"))
        i
    end
end

# ── The lattice ───────────────────────────────────────────────────────────────────────────

# A trinomial lattice of the Hull–White state x = r - α(t), on the levels `times` (t₀ = 0, …, t_N).
# Level ℓ (1-based) holds the nodes `lo[ℓ]:hi[ℓ]`, at x = j·dx[ℓ]. From node j of level ℓ the state
# branches to nodes k-1, k, k+1 of level ℓ+1 (`k[ℓ][j]`) with probabilities `p[ℓ][j]` (down, middle,
# up). Over the step from level ℓ, the rate at node j is α[ℓ] + j·dx[ℓ].
struct __Lattice{A}
    times::Vector{Float64}
    lo::Vector{Int}
    hi::Vector{Int}
    dx::Vector{A}
    k::Vector{Vector{Int}}
    p::Vector{Vector{NTuple{3, A}}}
    α::Vector{A}
end

# The node grid of the levels `times` for the model `hw`, with every level spanning at least
# ±`width` in x (0 for the reachable nodes only). Hull–White (1996) variable steps: the spacing of
# level ℓ+1 is √3 standard deviations of the step into it, and each node branches to the three nodes
# nearest its conditional mean, with probabilities that match its mean and variance (positive for any
# step and any mean reversion). The topology depends on primal values only; the probabilities and
# spacings carry the model's partials. α is calibrated by forward induction on the state prices,
# normalized by each level's discount factor, so that the lattice reprices the curve at every level.
function __build_lattice(hw::ShortRate.HullWhite, times::Vector{Float64}, width)
    a, σ, curve = hw.a, hw.σ, hw.curve
    n = length(times)
    A = typeof(float(one(a) * one(σ) * Yield.__log_interval(curve, times[1], times[end])))
    dual = __ad_depth(A) > 0
    lo, hi = zeros(Int, n), zeros(Int, n)
    dx = zeros(A, n)
    ks = Vector{Vector{Int}}(undef, n - 1)
    ps = Vector{Vector{NTuple{3, A}}}(undef, n - 1)
    deterministic = iszero(__primal(σ))
    for ℓ in 1:(n - 1)
        if deterministic
            # one state, x = 0: every node of a level coincides, and the middle branch keeps x = 0
            m = hi[ℓ] - lo[ℓ] + 1
            ks[ℓ], ps[ℓ] = zeros(Int, m), fill((zero(A), one(A), zero(A)), m)
            lo[ℓ + 1], hi[ℓ + 1] = -1, 1
            continue
        end
        Δt = times[ℓ + 1] - times[ℓ]
        dx[ℓ + 1] = σ * sqrt(3 * __decay_integral(2a, Δt))
        decay = exp(-a * Δt)
        nodes = lo[ℓ]:hi[ℓ]
        kℓ = Vector{Int}(undef, length(nodes))
        pℓ = Vector{NTuple{3, A}}(undef, length(nodes))
        for (i, j) in enumerate(nodes)
            μ = j * dx[ℓ] * decay / dx[ℓ + 1]   # the conditional mean, in units of the next spacing
            k = round(Int, __primal(μ))
            ξ = μ - k
            (dual && abs(__primal(ξ)) == 0.5) && throw(
                ArgumentError("a lattice node's mean lies exactly between two nodes, where the lattice's branching switches; its derivatives are undefined there")
            )
            kℓ[i] = k
            pℓ[i] = (1 / 6 + (ξ^2 - ξ) / 2, 2 / 3 - ξ^2, 1 / 6 + (ξ^2 + ξ) / 2)
        end
        J = width > 0 ? ceil(Int, width / __primal(dx[ℓ + 1])) : 0
        lo[ℓ + 1], hi[ℓ + 1] = min(minimum(kℓ) - 1, -J), max(maximum(kℓ) + 1, J)
        ks[ℓ], ps[ℓ] = kℓ, pℓ
    end
    α = Vector{A}(undef, n - 1)
    q = A[one(A)]   # state prices of level ℓ divided by the discount factor to tℓ
    for ℓ in 1:(n - 1)
        Δt = times[ℓ + 1] - times[ℓ]
        ΔL = Yield.__log_interval(curve, times[ℓ], times[ℓ + 1])
        x(i) = (lo[ℓ] + i - 1) * dx[ℓ]
        α[ℓ] = (log(sum(q[i] * exp(-x(i) * Δt) for i in eachindex(q))) + ΔL) / Δt
        next = zeros(A, hi[ℓ + 1] - lo[ℓ + 1] + 1)
        g = exp(ΔL)
        for i in eachindex(q)
            w = q[i] * exp(-(α[ℓ] + x(i)) * Δt) * g
            b = ks[ℓ][i] - lo[ℓ + 1] + 1
            pd, pm, pu = ps[ℓ][i]
            next[b - 1] += w * pd
            next[b] += w * pm
            next[b + 1] += w * pu
        end
        q = next
    end
    return __Lattice(times, lo, hi, dx, ks, ps, α)
end

# The levels of a lattice for an asset: time 0 and every payment date, each gap split into equal
# steps no longer than `max_step`.
function __lattice_times(rows, max_step)
    times = [0.0]
    for r in rows
        a, b = last(times), float(__primal(r.time))
        m = max(1, ceil(Int, (b - a) / max_step - 1.0e-9))
        for j in 1:m
            push!(times, j == m ? b : a + j * (b - a) / m)
        end
    end
    return times
end

# ── Rolling an asset back ─────────────────────────────────────────────────────────────────

# (1 - e⁻ˣ)/x, the average of e⁻ᵘ over [0, x], with its limit 1 at 0.
__mean_decay(x) = inv(Bond.__y_over_expm1(x))

# The pieces of the step [t0, t1] between the hazard curve's breakpoints:
# `(start, length, intensity, survival from t0 to start)`.
function __hazard_pieces(h, t0, t1)
    cuts = __recovery_cuts(t0, t1, Credit.__breakpoints(h))
    return map(1:(length(cuts) - 1)) do i
        s, e = cuts[i], cuts[i + 1]
        (s, e - s, Credit.hazard(h, s), Credit.survival(h, t0, s))
    end
end

# The value per unit of the asset alive at time 0, and the continuation values (per unit continuing,
# before exercise) at the levels in `keep`. `events[ℓ]` is the row paid at level ℓ (or nothing);
# `exercise[ℓ]` the redemption value of each continuing unit and whether the issuer (call) or holder
# (put) decides; `opening[ℓ]` the principal outstanding over the step from level ℓ.
function __rollback(lat, rows, events, exercise, opening, credit; keep = ())
    n = length(lat.times)
    # every input a node value can carry partials from: the lattice, the payments and balances, the
    # redemption prices and the credit assumptions
    T = promote_type(
        eltype(lat.α), typeof(first(rows).payment), typeof(first(rows).continuing), eltype(opening),
        mapreduce(e -> typeof(e.redemption), promote_type, values(exercise); init = Bool),
        credit === nothing ? Bool : typeof(credit.recovery * Credit.survival(credit.hazard, 0.0, lat.times[end]))
    )
    dual = __ad_depth(T) > 0
    V = zeros(T, lat.hi[n] - lat.lo[n] + 1)   # nothing remains after the last payment
    kept = Dict{Int, Vector{T}}()
    for ℓ in n:-1:2
        C = V
        p = events[ℓ]
        if p !== nothing
            r = rows[p]
            ex = get(exercise, ℓ, nothing)
            ex === nothing || ℓ in keep && (kept[ℓ] = copy(V))
            C = map(V) do v
                if ex !== nothing
                    (dual && __primal(v) == __primal(ex.redemption)) && throw(
                        ArgumentError("an exercise decision at time $(lat.times[ℓ]) is exactly at indifference; its derivative is undefined")
                    )
                    v = ex.call ? min(v, ex.redemption) : max(v, ex.redemption)
                end
                r.payment + r.continuing * v
            end
        end
        t0, t1 = lat.times[ℓ - 1], lat.times[ℓ]
        Δt = t1 - t0
        pieces = credit === nothing ? nothing : __hazard_pieces(credit.hazard, t0, t1)
        survival = credit === nothing ? one(T) : Credit.survival(credit.hazard, t0, t1)
        base = credit === nothing ? zero(T) : credit.recovery * opening[ℓ - 1]
        nodes = lat.lo[ℓ - 1]:lat.hi[ℓ - 1]
        Vp = Vector{T}(undef, length(nodes))
        for (i, j) in enumerate(nodes)
            rate = lat.α[ℓ - 1] + j * lat.dx[ℓ - 1]
            b = lat.k[ℓ - 1][i] - lat.lo[ℓ] + 1
            pd, pm, pu = lat.p[ℓ - 1][i]
            v = exp(-rate * Δt) * survival * (pd * C[b - 1] + pm * C[b] + pu * C[b + 1])
            if credit !== nothing
                # recovery of the outstanding principal at default within the step, at this node's rate:
                # ∫ λ S(t0, u) e^{-r(u - t0)} du, piece by piece, each an exact exponential
                for (s, len, λ, S) in pieces
                    v += base * λ * S * exp(-rate * (s - t0)) * len * __mean_decay((λ + rate) * len)
                end
            end
            Vp[i] = v
        end
        V = Vp
    end
    return only(V), kept
end

# The inputs of a rollback for the exercisable asset `c` valued in `ctx`.
function __exercise_inputs(ctx, c, max_step)
    spec = __asset_spec(c)
    rows = __asset_rows(spec.base, ctx)
    credit = spec.key === nothing ? nothing : ctx[spec.key]
    ex_rows = __exercise_rows(spec.exercise, rows)
    times = __lattice_times(rows, max_step)
    level = Dict(float(__primal(r.time)) => findfirst(==(float(__primal(r.time))), times) for r in rows)
    n = length(times)
    events = Vector{Union{Nothing, Int}}(nothing, n)
    for (p, r) in enumerate(rows)
        events[level[float(__primal(r.time))]] = p
    end
    call = spec.exercise isa Bond.Callable
    exercise = Dict(
        level[float(__primal(rows[p].time))] => (; call, redemption = spec.exercise.prices[e] * rows[p].closing)
            for (e, p) in enumerate(ex_rows)
    )
    # the principal outstanding over each step: the opening balance of the next payment
    opening = Vector{typeof(first(rows).opening)}(undef, n - 1)
    p = 1
    for ℓ in 1:(n - 1)
        while times[ℓ + 1] > float(__primal(rows[p].time))
            p += 1
        end
        opening[ℓ] = rows[p].opening
    end
    return (; spec, rows, credit, times, events, exercise, opening)
end

# Second derivatives through an exercise decision would leave out the moving exercise boundary, and
# the lattice has no derivative in σ at σ = 0 (it degenerates to one node per level).
function __check_lattice_derivatives(hw, depth)
    depth >= 2 && throw(
        ArgumentError(
            "second derivatives through an exercise decision are not supported: differentiating the lattice along its " *
                "decisions leaves out the moving exercise boundary (a callable bond's negative convexity). " *
                "Measure convexity by revaluing at bumped inputs."
        )
    )
    (iszero(__primal(hw.σ)) && !iszero(hw.σ)) && throw(
        ArgumentError("an exercisable asset's value is not differentiable in σ at σ = 0, where the lattice has one node per level")
    )
    return nothing
end

function __exercise_value(tree::ShortRate.HullWhiteTree, ctx, c)
    inputs = __exercise_inputs(ctx, c, tree.max_step)
    lat = __build_lattice(tree.model, inputs.times, 0.0)
    value, _ = __rollback(lat, inputs.rows, inputs.events, inputs.exercise, inputs.opening, inputs.credit)
    __check_lattice_derivatives(tree.model, __ad_depth(value))
    return value
end
__exercise_value(m, ctx, c) = throw(
    ArgumentError(
        "an exercisable asset is valued on a lattice: use ShortRate.HullWhiteTree(model) as the valuation model, not $(nameof(typeof(m)))"
    )
)

FinanceCore.present_value(ctx, c::Union{Bond.Callable, Bond.Putable}) = __exercise_value(valuation_model(ctx), ctx, c)
FinanceCore.present_value(ctx, c::Credit.Defaultable{<:Union{Bond.Callable, Bond.Putable}}) = __exercise_value(valuation_model(ctx), ctx, c)

# A lattice valuation model prices closed-form contracts by its Hull–White model's formulas.
closed_form(m::ShortRate.HullWhiteTree, c::__ClosedFormContract) = closed_form(m.model, c)

# ── The exercise policy ───────────────────────────────────────────────────────────────────

"""
    ExercisePolicy

The exercise decisions of an asset's lattice, from [`exercise_policy`](@ref). Calling it with an
exercise time and the short rate observed then, `policy(t, r)`, returns `true` when the asset is
exercised: called when its redemption value is below the value of continuing, put when it is above.
The rate is a number, read as continuously compounded, or a `Rate`.

The policy reads only the time and the current rate, never later rates. It maps the rate to the
lattice's state with the Hull–White drift of the pricing model at `t`, and interpolates the
continuation value linearly between the lattice's nodes. Its nodes span ±8 standard deviations of
the state at maturity at every exercise time; a rate beyond them throws an `ArgumentError`. The
decisions are exact for rate paths of the same Hull–White model; for paths of another model they
are the pricing model's decisions at the observed rate.
"""
struct ExercisePolicy{M}
    call::Bool
    model::M
    times::Vector{Float64}
    redemption::Vector{Float64}
    lo::Vector{Int}
    dx::Vector{Float64}
    continuation::Vector{Vector{Float64}}
end

"""
    exercise_policy(ctx, asset) -> ExercisePolicy

The exercise decisions that value the exercisable `asset` in `ctx`, whose valuation model is a
[`ShortRate.HullWhiteTree`](@ref FinanceModels.ShortRate.HullWhiteTree): the continuation value at
each exercise time as a function of the short rate, on a lattice that spans ±8 standard deviations of
the state at maturity at every exercise time, so that projections along simulated rate paths can
apply them. The asset's value is [`present_value`](@ref FinanceCore.present_value)`(ctx, asset)`.

```julia
policy = exercise_policy(Models(ShortRate.HullWhiteTree(hw)), bond)
policy(5.0, 0.07)   # whether the bond is exercised at time 5 when the short rate is 7%
```
"""
function exercise_policy(ctx, c)
    tree = valuation_model(ctx)
    tree isa ShortRate.HullWhiteTree || __exercise_value(tree, ctx, c)
    hw = tree.model
    inputs = __exercise_inputs(ctx, c, tree.max_step)
    width = 8 * __primal(hw.σ) * sqrt(__primal(__decay_integral(2hw.a, last(inputs.times))))
    lat = __build_lattice(hw, inputs.times, width)
    levels = sort!(collect(keys(inputs.exercise)))
    _, kept = __rollback(lat, inputs.rows, inputs.events, inputs.exercise, inputs.opening, inputs.credit; keep = levels)
    return ExercisePolicy(
        inputs.spec.exercise isa Bond.Callable, hw, [inputs.times[ℓ] for ℓ in levels],
        [float(__primal(inputs.exercise[ℓ].redemption)) for ℓ in levels], [lat.lo[ℓ] for ℓ in levels],
        [float(__primal(lat.dx[ℓ])) for ℓ in levels], [map(v -> float(__primal(v)), kept[ℓ]) for ℓ in levels],
    )
end

function (p::ExercisePolicy)(t, r)
    i = findfirst(s -> abs(s - __primal(t)) <= 1.0e-9 * max(1, abs(s)), p.times)
    i === nothing && throw(ArgumentError("time $t is not an exercise time of the policy; its exercise times are $(p.times)"))
    x = __primal(__continuous(r)) - __primal(__hw_alpha(p.model, p.times[i]))
    v = __continuation(p, i, x)
    return p.call ? p.redemption[i] < v : p.redemption[i] > v
end

# The continuation value at state x of exercise time i, linear between the lattice's nodes.
function __continuation(p::ExercisePolicy, i, x)
    values, dx = p.continuation[i], p.dx[i]
    if iszero(dx)   # zero volatility: the one deterministic state
        abs(x) <= 1.0e-10 || throw(ArgumentError("the short rate is off the zero-volatility lattice's path by $x"))
        return values[1 - p.lo[i]]
    end
    u = x / dx - p.lo[i] + 1   # the fractional node index
    1 <= u <= length(values) || throw(
        ArgumentError(
            "the short rate at time $(p.times[i]) is beyond the exercise policy's lattice (state $x outside " *
                "[$(p.lo[i] * dx), $((p.lo[i] + length(values) - 1) * dx)])"
        )
    )
    j = min(floor(Int, u), length(values) - 1)
    w = u - j
    return (1 - w) * values[j] + w * values[j + 1]
end
