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

# The fixed cashflows of a quote's instrument: a contract's own, or those of every member of a
# portfolio, which is an array or a tuple of contracts (or of portfolios). An empty portfolio has
# none. An instrument whose cashflows depend on a model, such as a floating leg, fails in `collect`.
__quote_cashflows(contract) = collect(contract)
__quote_cashflows(portfolio::Union{AbstractArray, Tuple}) = mapreduce(__quote_cashflows, vcat, portfolio; init = Union{}[])

# `(time, amount)` pairs, sorted by the value of the time, with the amounts of equal times summed.
# Payments at one time share their discount factor, so only their sum counts. `==` between dual
# numbers also compares their partials, so two dates that coincide in value but move differently
# stay apart: their payments offset in value only, not in derivative.
function __sum_by_time(dated)
    summed = empty(dated)
    for (s, a) in dated
        i = lastindex(summed)
        while i >= firstindex(summed) && __primal(first(summed[i])) == __primal(s) && first(summed[i]) != s
            i -= 1
        end
        if i >= firstindex(summed) && first(summed[i]) == s
            summed[i] = (s, last(summed[i]) + a)
        else
            push!(summed, (s, a))
        end
    end
    return summed
end

# A quote's payments with their dual parts, one per time. The amounts are floating-point numbers,
# so that integer amounts cannot overflow in a sum.
__quote_payments(instrument) =
    __sum_by_time(sort!([(cf.time, float(cf.amount)) for cf in __quote_cashflows(instrument)]; by = __primal ∘ first))

# The payments the solve works on: values only, one per time.
__primal_payments(payments) = __sum_by_time([(__primal(s), __primal(a)) for (s, a) in payments])

# Structural condition (S1) of `quote_shift` for the quote at tenor `t`, on its payments with their
# dual parts. A dual number is zero only if its partials are too, and equal only to a number with
# the same partials. So a payment after the tenor that is zero in value only, or that another
# offsets only where their dates coincide, is rejected like any other: a later tenor's bump would
# move its derivative after this quote's solve.
#
# A payment date may move with a differentiated quantity (a family whose dates depend on the quote),
# but not when it is exactly the tenor. The bumps have a kink in time at every tenor, and a dual
# time there reads the piece that starts at it: the next tenor's, in which this quote may not pay
# at all. The derivative would describe a move the quote cannot make.
function __check_late_payments(payments, t, target)
    for (s, a) in payments
        iszero(a) && continue
        (__primal(s) == t && s != t) && throw(
            ArgumentError(
                "quote_shift: the quote at tenor $t pays exactly at its tenor on a date that moves with a differentiated " *
                    "quantity (the quote, the curve or the shock). Its derivative is not defined there: any later, the payment " *
                    "would belong to the next tenor's bump. A moving payment date must stay before its tenor."
            )
        )
        __primal(s) > t && throw(
            ArgumentError(
                "quote_shift: the quote at tenor $t (target $target) pays at time $(__primal(s)), after its tenor. " *
                    "Each quote may pay only up to its own tenor, so that a later tenor's bump cannot move it."
            )
        )
    end
    return nothing
end

# The structural conditions (S2) and (S3) of `quote_shift` for the quote at tenor `t`, whose bump
# moves only the payments in `(t_prev, t]`. They are about the root's value, so they read the
# primal payments, zero amounts dropped.
function __check_quote_payments(payments, t_prev, t, target)
    window = [a for (s, a) in payments if t_prev < s <= t && !iszero(a)]
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
# bracket `(0, 0)`, which is the root. Otherwise the probes step toward the root by `2ⁱ/t`; the bump
# moves each payment's log discount by at most `x·t`, so the first probe changes each payment's
# value by at most a factor of e. Upward, every term decays and stays finite. Downward, terms grow
# and can overflow before the root; a non-finite probe is never a bracket end, even `+Inf`, and the
# search bisects back toward the last finite one. The probes themselves are capped at the largest
# finite number, which is also the first probe of a tenor so short that `1/t` overflows its type: a
# residual that keeps its sign and stays finite all the way there has no bracket.
function __quote_bracket(r, t)
    x0 = zero(float(t))
    g0 = r(x0)
    iszero(g0) && return (x0, x0)
    dir = g0 > 0 ? one(x0) : -one(x0)
    top = floatmax(x0)
    lo, step = x0, min(inv(float(t)), top)
    hi = dir * step
    while true
        v = r(hi)
        isfinite(v) || break
        sign(v) == sign(g0) || return minmax(lo, hi)
        lo = hi
        step == top && throw(
            ArgumentError(
                "quote_shift could not bracket the bump for the quote at tenor $t within floating-point limits: " *
                    "the residual keeps its sign up to $lo, the largest finite probe."
            )
        )
        step = min(2 * step, top)
        hi = dir * step
    end
    # halves, so that the midpoint of two large probes does not overflow
    m = lo / 2 + hi / 2
    while m != lo && m != hi
        v = r(m)
        if !isfinite(v)
            hi = m
        elseif sign(v) == sign(g0)
            lo = m
        else
            return minmax(lo, m)
        end
        m = lo / 2 + hi / 2
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
a price for `ZCBPrice`. `curve` is a yield model. `family` is a quote constructor taking
`(quote, maturity)`, as in [`implied_quote`](@ref FinanceModels.Yield.implied_quote). Its
instrument has fixed cashflows, and may be a portfolio: an array or a tuple of contracts, whose
payments are taken together. `tenors` are positive, strictly increasing, and not dual numbers.

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
  condition, as does a bump that cannot be bracketed within floating-point limits. Under
  differentiation, condition 1 also rejects a payment after `tₖ` that is zero in value but not
  in derivative (an amount, or two offsetting payments whose dates coincide but move
  differently): a later bump would move that derivative.

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
        target = __primal(news[k])
        dual_payments = __quote_payments(q.instrument)
        __check_late_payments(dual_payments, t, target)
        payments = __primal_payments(dual_payments)
        __check_quote_payments(payments, t_prev, t, target)
        solved = primal + __trial_quote_bumps(pb, ts)
        fixed = sum(a * discount(solved, s) for (s, a) in payments if s <= t_prev; init = zero(qp.price))
        qp.price > fixed || throw(
            ArgumentError(
                "quote_shift: the target $target at tenor $t is unreachable: its price $(qp.price) does not exceed " *
                    "$fixed, the value of its payments at or before $t_prev, which its bump cannot move."
            )
        )
        # Solve relative to the quote's size: the gross value of its payments, each in absolute
        # value. The net value is no measure of size where payments offset one another, and a
        # target price of zero would leave no tolerance at the root at all. The residual's own unit
        # `σ` is the larger of that size before the bump and the target price: a price many orders
        # above today's value, divided by today's value alone, overflows before the search starts.
        sizes = [(s, abs(a)) for (s, a) in payments]
        function magnitude(x)
            c = primal + __trial_quote_bumps(__with_bump(pb, k, x), ts)
            return sum(a * discount(c, s) for (s, a) in sizes; init = zero(qp.price))
        end
        s0 = max(magnitude(zero(t)), abs(qp.price))
        σ = iszero(s0) ? one(s0) : s0
        g(x) = (value(curve, __with_bump(bumps, k, x), q.instrument) - q.price) / σ
        g_primal(x) = (value(primal, __with_bump(pb, k, x), qp.instrument) - qp.price) / σ
        scale(x) = magnitude(x) / σ
        # The slope is the bump's own effect, so it is measured against the payments the bump
        # moves. Against all of them, large payments that offset one another before the previous
        # tenor would make a sound slope look like none.
        moved = [(s, a) for (s, a) in sizes if s > t_prev]
        function slope_scale(x)
            c = primal + __trial_quote_bumps(__with_bump(pb, k, x), ts)
            return sum(a * discount(c, s) for (s, a) in moved; init = zero(qp.price)) / σ
        end
        # The root is solved inside its bracket. A secant search from zero stops where this
        # residual is small in absolute terms, which a long quote under a large shock reaches far
        # from its root: `exp(-x·t)` is below every absolute tolerance once `x·t` passes about 30.
        residual = __PrimalResidual(g_primal, who, hint)
        lo, hi = __quote_bracket(residual, t)
        x = lo == hi ? lo : Roots.find_zero(residual, (lo, hi), Roots.A42())
        bumps[k] = __implicit_step(g, residual, x; scale, slope_scale, who)
    end
    return curve + __quote_bumps(bumps, ts)
end
