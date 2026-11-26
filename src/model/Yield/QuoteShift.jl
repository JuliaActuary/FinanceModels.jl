# ─── Quote shifts ─────────────────────────────────────────────────────────

# The bumps of a quote shift: continuous zero rates interpolated linearly between the tenors and
# held flat beyond both ends, so the bump at one tenor is a hat from the previous tenor to the next.
# The solve's trial curves skip validation, like a fit's; the returned curve is built through
# `ZeroRateCurve`.
__quote_bumps(b, ts) = ZeroRateCurve(b, ts, Sp.Linear(); extrapolation = :flat_zero)
__trial_quote_bumps(b, ts) = __build(Sp.Linear(), KnotGrid(Unchecked(), b, ts); extrapolation = :flat_zero)

# A copy of `b` with entry `k` replaced by `x`, in an element type holding both.
function __with_bump(b, k, x)
    v = Vector{promote_type(eltype(b), typeof(x))}(b)
    v[k] = x
    return v
end

# A quote's payments as sorted `(time, amount)` pairs of primal values, one per time, with zero
# amounts dropped. Payments at one time share their discount factor, so only their sum counts. An
# instrument whose cashflows depend on a model, such as a floating leg, fails in `collect`.
function __quote_payments(instrument)
    ps = sort!([(__primal(cf.time), __primal(cf.amount)) for cf in collect(instrument)]; by = first)
    merged = empty(ps)
    for (s, a) in ps
        if !isempty(merged) && first(last(merged)) == s
            merged[end] = (s, last(merged[end]) + a)
        else
            push!(merged, (s, a))
        end
    end
    return filter!(p -> !iszero(last(p)), merged)
end

# The structural conditions (S1)–(S3) of `quote_shift` for the quote at tenor `t`, whose bump moves
# only the payments in `(t_prev, t]`.
function __check_quote_payments(payments, t_prev, t, target)
    i = findfirst(p -> first(p) > t, payments)
    isnothing(i) || throw(
        ArgumentError(
            "quote_shift: the quote at tenor $t (target $target) pays at time $(first(payments[i])), after its tenor. " *
                "Each quote may pay only up to its own tenor, so that a later tenor's bump cannot move it."
        )
    )
    window = [a for (s, a) in payments if s > t_prev]
    isempty(window) && throw(
        ArgumentError(
            "quote_shift: the quote at tenor $t (target $target) pays nothing after the previous tenor $t_prev, " *
                "so its bump cannot move it."
        )
    )
    (issorted(window; by = sign) && last(window) > 0) || throw(
        ArgumentError(
            "quote_shift: the quote at tenor $t (target $target) pays $window after the previous tenor $t_prev. " *
                "Those amounts must change sign at most once, from negative to positive, and end positive, " *
                "for the target to determine its bump uniquely."
        )
    )
    return nothing
end

# A bracket around the unique root of the primal residual `r` of the bump at tenor `t`. The residual
# is positive below the root and negative above it. An exact root at zero gives the degenerate
# bracket `(0, 0)`, which the solve never uses: its secant stops at the starting point. Otherwise the
# probes step toward the root by `2ⁱ/t`; the bump moves each payment's log discount by at most `x·t`,
# so the first probe changes each payment's value by at most a factor of e. Upward, every term
# decays and stays finite. Downward, terms grow and can overflow before the root; a non-finite
# probe is never a bracket end, even `+Inf`, and the search bisects back toward the last finite one.
function __quote_bracket(r, t)
    x0 = zero(float(t))
    g0 = r(x0)
    iszero(g0) && return (x0, x0)
    dir = g0 > 0 ? one(x0) : -one(x0)
    lo, step = x0, inv(float(t))
    while true
        x = dir * step
        v = r(x)
        isfinite(v) || break
        sign(v) == sign(g0) || return minmax(lo, x)
        lo, step = x, 2 * step
    end
    hi = dir * step
    m = (lo + hi) / 2
    while m != lo && m != hi
        v = r(m)
        if !isfinite(v)
            hi = m
        elseif sign(v) == sign(g0)
            lo = m
        else
            return minmax(lo, m)
        end
        m = (lo + hi) / 2
    end
    throw(
        ArgumentError(
            "quote_shift could not bracket the bump for the quote at tenor $t within floating-point limits: " *
                "the residual keeps its sign up to $lo and is not finite at $hi."
        )
    )
end

"""
    quote_shift(rule, curve, family, tenors)

Shift the quotes of `family` at `tenors` by `rule`, and return `curve` plus continuous zero-rate
bumps at `tenors` that reprice each shifted quote.

`rule(q, t)` receives the quote `q = implied_quote(curve, family, t)` at each tenor `t` and returns
the new quote, in the family's own units: a rate number for `ParYield`, `CMTYield` or `ZCBYield`,
a price for `ZCBPrice`. `family` is a quote constructor taking `(quote, maturity)`, as in
[`implied_quote`](@ref FinanceModels.Yield.implied_quote). `tenors` are positive, strictly
increasing, and not dual numbers.

The result is `curve + ZeroRateCurve(bumps, tenors, Spline.Linear(); extrapolation = :flat_zero)`,
a [`CompositeYield`](@ref). Each bump is a continuous zero-rate hat, linear between neighboring
tenors and flat beyond the first and last, so one tenor gives a parallel shift. The bumps are
solved one tenor at a time, shortest first: the quote at a tenor depends only on the bumps up to it.

# Contract

- **Targets.** Each shifted quote reprices on the result within solver tolerance, at the tenors
  only. Between tenors, quotes move approximately, and forwards jump at the tenors.
- **Identity.** A rule returning `q` unchanged gives back `curve` within solver tolerance, not
  bit for bit.
- **Shape.** The result does not keep `curve`'s shape guarantees, such as positive forwards.
- **Derivatives.** First-order ForwardDiff derivatives are exact, by the implicit function theorem
  at each solve. They flow through `curve`'s parameters and the rule's new quotes (including a
  quote's price), also when a shock's value is zero. The value of a dual result equals the
  non-dual result exactly. Nested dual numbers throw an `ArgumentError`.
- **Accepted quotes.** Each quote is checked before its bump is solved. With `t₀ = 0` and the
  quote at `tₖ`, zero amounts ignored, it must:
  1. pay only at or before `tₖ`;
  2. pay something in `(tₖ₋₁, tₖ]`;
  3. there, change sign at most once, from negative to positive, and end positive;
  4. have a price above the value of its payments at or before `tₖ₋₁`, which the bump cannot
     move (payments at time 0 count).

  Then, with positive discount factors, the target determines the bump uniquely. Every standard
  family meets 1–3 at any rates where its final payment is positive, including negative rates
  with negative coupons. A violation throws an `ArgumentError` naming the tenor and the
  condition, as does a bump that cannot be bracketed within floating-point limits.

The derivative at a zero shock is local quote risk on `curve`. To get risk to quotes that
calibrate a curve, differentiate through [`fit`](@ref FinanceModels.fit) instead.

# Examples

```julia
tenors = [1, 2, 3, 5, 7, 10, 20, 30]
par_family(y, t) = ParYield(y, t; frequency = 2)
curve = fit(Spline.MonotoneConvex(), par_family.(0.04, tenors))

# every quoted par yield +100bp
up = quote_shift(curve, par_family, tenors) do q, t
    q + 0.01
end
implied_quote(up, par_family, 10.0)  # ≈ 0.05

# a 10-year CMT key-rate bump
quote_shift(curve, CMTYield, tenors) do q, t
    t == 10 ? q + 0.0001 : q
end
```

See also [`implied_quote`](@ref FinanceModels.Yield.implied_quote),
[`TenorShift`](@ref FinanceModels.Yield.TenorShift), [`fit`](@ref FinanceModels.fit).
"""
function quote_shift(rule::R, curve::AbstractYieldModel, family::F, tenors) where {R, F}
    ts = __owned_float_vector(tenors)
    __ad_depth(eltype(ts)) == 0 || throw(
        ArgumentError(
            "quote_shift cannot differentiate with respect to the tenors (got $ts); " *
                "only the curve and the rule's new quotes may carry dual numbers."
        )
    )
    grid = KnotGrid(zeros(eltype(ts), length(ts)), ts, Sp.Linear(); who = "quote_shift")   # validates the tenors
    zero_bumps = __build(Sp.Linear(), grid; extrapolation = :flat_zero)
    iszero(first(ts)) && throw(
        ArgumentError("quote_shift: a quote maturing at t = 0 does not depend on its bump; remove the zero tenor.")
    )
    news = map(t -> rule(implied_quote(curve, family, t), t), ts)
    targets = map(family, news, ts)
    primal_targets = map((x, t) -> family(__primal(x), t), news, ts)
    # The bumps hold every partial of every residual: the curve's, and those of each new quote,
    # which some families carry only in the price.
    B = mapreduce(q -> typeof(FinanceCore.present_value(curve + zero_bumps, q.instrument) - q.price), promote_type, targets)
    bumps = zeros(B, length(ts))
    primal = __PrimalCurve(curve)
    who = "quote_shift"
    hint = "quote_shift differentiates through the curve and the rule's new quotes only, not through `family`."
    for k in eachindex(ts)
        t, t_prev = ts[k], k == 1 ? zero(ts[k]) : ts[k - 1]
        q, qp = targets[k], primal_targets[k]
        pb = map(__primal, bumps)
        value(c, b, instrument) = FinanceCore.present_value(c + __trial_quote_bumps(b, ts), instrument)
        payments = __quote_payments(qp.instrument)
        __check_quote_payments(payments, t_prev, t, news[k])
        solved = primal + __trial_quote_bumps(pb, ts)
        fixed = sum(a * discount(solved, s) for (s, a) in payments if s <= t_prev; init = zero(qp.price))
        qp.price > fixed || throw(
            ArgumentError(
                "quote_shift: the target $(news[k]) at tenor $t is unreachable: its price $(qp.price) does not exceed " *
                    "$fixed, the value of its payments at or before $t_prev, which its bump cannot move."
            )
        )
        # Solve relative to the quote's size, as `implied_quote` does.
        magnitude(x) = max(abs(value(primal, __with_bump(pb, k, x), qp.instrument)), abs(qp.price))
        s0 = magnitude(zero(t))
        σ = iszero(s0) ? one(s0) : s0
        g(x) = (value(curve, __with_bump(bumps, k, x), q.instrument) - q.price) / σ
        g_primal(x) = (value(primal, __with_bump(pb, k, x), qp.instrument) - qp.price) / σ
        scale(x) = magnitude(x) / σ
        bracket = __quote_bracket(__PrimalResidual(g_primal, who, hint), t)
        bumps[k] = __implicit_root(g, g_primal, zero(t); bracket, scale, who, hint)
    end
    return curve + __quote_bumps(bumps, ts)
end
