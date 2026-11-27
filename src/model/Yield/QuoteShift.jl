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

# Structural condition (S1) of `quote_shift` for the quote at tenor `t`, on the instrument's
# cashflows one at a time, before any are summed. A later tenor's bump moves every discount factor
# after `t`, so a payment there moves after this quote's solve. That holds for two payments on one
# later date that offset each other (the quote's value adds them in the instrument's own order,
# where their difference need not survive), and for a payment of nothing: its amount can carry a
# derivative, and a discount factor that a later bump raises to `Inf` makes its value `NaN`.
#
# A payment date may not carry a derivative. The bumps' zero rates have a kink in time at every
# tenor, so a date that moves onto one has a derivative on each side, and ForwardDiff takes the
# side from the direction of the first seed: one gradient would mix the two. A payment of nothing
# at or before `t` has no value on any date. (A dual number is zero only if its partials are too.)
function __check_quote_cashflows(cashflows, t, target)
    for cf in cashflows
        s = cf.time
        __primal(s) > t && throw(
            ArgumentError(
                "quote_shift: the quote at tenor $t (target $target) pays at time $(__primal(s)), after its tenor. " *
                    "Each quote may pay only up to its own tenor, so that a later tenor's bump cannot move it."
            )
        )
        iszero(cf.amount) && continue
        s == __primal(s) || throw(
            ArgumentError(
                "quote_shift: the quote at tenor $t (target $target) pays at time $(__primal(s)), a date that depends on " *
                    "a differentiated quantity (the quote, the curve or the shock). Payment dates must be fixed: the " *
                    "bumps have a kink in time at every tenor."
            )
        )
    end
    return nothing
end

# `(time, amount, gross)` of a quote's cashflows, values only, sorted by time, the amounts of one
# time summed: payments at one time share their discount factor. `gross` is the sum of their
# absolute values, which is the amount's size where two of them offset each other. The amounts are
# floating-point numbers, so that integer amounts cannot overflow in a sum. The conditions on the
# root (S2 to S5) and its sizes read these.
function __quote_payments(cashflows)
    dated = sort!([(__primal(cf.time), float(__primal(cf.amount))) for cf in cashflows]; by = first)
    summed = [(s, a, abs(a)) for (s, a) in dated[1:0]]
    for (s, a) in dated
        if !isempty(summed) && summed[end][1] == s
            summed[end] = (s, summed[end][2] + a, summed[end][3] + abs(a))
        else
            push!(summed, (s, a, abs(a)))
        end
    end
    return summed
end

# The last payment of a quote that is not zero: positive by (S3), and always moved by the bump.
__last_payment(payments) = payments[findlast(p -> !iszero(p[2]), payments)]

# Two sizes of the quote at tenor `t` on curve `c`, for (S5). `gross` is the value of all its
# payments, each on its own and in absolute value: the valuation's rounding scales with it.
# `exposure` is the derivative of the quote's value to the log discount factor that the bump sets
# at `t`, so the value's slope in the bump is `-exposure·t`. The bump is a hat in the zero rate, so
# it moves the log discount factor at `s` by `s·h(s)`, where `h` rises from 0 at `t_prev` to 1 at
# `t` (and is 1 up to the first tenor): a payment just after `t_prev` barely moves. At a root that
# (S3) and (S4) make unique, the exposure is positive.
function __quote_sizes(c, payments, t_prev, t, first_tenor)
    gross = exposure = zero(payments[1][3] * discount(c, payments[1][1]))
    for (s, a, g) in payments
        d = discount(c, s)
        gross += g * d
        h = first_tenor ? one(s) : s > t_prev ? (s - t_prev) / (t - t_prev) : zero(s)
        exposure += a * d * (s / t) * h
    end
    return gross, exposure
end

# The unit a quote's residual is measured in on curve `c`, where its value is `net`: the mean of
# the sizes of that value, of the price, and of the value of the quote's last payment. The first
# two bound the residual; the third keeps the unit positive where a price of zero meets a value of
# zero. A mean, taken term by term, is smooth in the curve and cannot overflow while its terms are
# finite. The unit carries no dual part.
__quote_unit(c, net, price, (s_last, a_last, _)) = __primal(abs(net) / 3 + abs(price) / 3 + a_last * discount(c, s_last) / 3)

# The structural conditions (S2) and (S3) of `quote_shift` for the quote at tenor `t`, whose bump
# moves only the payments in `(t_prev, t]`. They are about the root's value, so they read the
# payments' values, zero amounts dropped.
function __check_quote_payments(payments, t_prev, t, target)
    window = [a for (s, a, _) in payments if t_prev < s <= t && !iszero(a)]
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

- **Targets.** Each shifted quote reprices on the result, by `present_value` of its instrument,
  within `sqrt(eps)` of the quote's size (the mean of the magnitudes of its price, of its value
  and of the value of its last payment), or the call throws an `ArgumentError`. That holds at
  the tenors only: between them quotes move approximately, and forwards jump at the tenors.
- **Identity.** A rule returning `q` unchanged gives back `curve` within solver tolerance, not
  bit for bit.
- **Shape.** The result does not keep `curve`'s shape guarantees, such as positive forwards.
- **Derivatives.** First-order ForwardDiff derivatives are exact where `curve` and the rule are
  differentiable, by the implicit function theorem at each solve. They flow through `curve`'s
  parameters and the rule's new quotes (including a quote's price), also when a shock's value
  is zero. Each bump is solved on the values of its inputs, and the derivative step does not
  move it. Nested dual numbers throw an `ArgumentError`, as does a payment date that depends on
  a differentiated quantity: the bumps' zero rates have a kink in time at every tenor, where
  the derivative would depend on the direction of the move.
- **Accepted quotes.** Each quote is checked as its bump is solved. With `t₀ = 0` and the quote
  at `tₖ`, it must:
  1. pay only at or before `tₖ`, every payment taken on its own, a payment of nothing
     included: two payments after `tₖ` that offset each other are two payments after `tₖ`;
  2. pay something in `(tₖ₋₁, tₖ]`, its payments summed by date;
  3. there, change sign at most once, from negative to positive, and end positive;
  4. have a price above the value of its payments at or before `tₖ₋₁`, which the bump cannot
     move (payments at time 0 count);
  5. be resolved by its valuation: at the solved bump, its value moves with the log discount
     factor at `tₖ` by at least `sqrt(eps)` of the value of all its payments, each taken on
     its own.

  Then, with positive discount factors, the target determines the bump uniquely, to about half
  the digits of the valuation or better. Every standard family meets 1–3 at any rates where its
  final payment is positive, including negative rates with negative coupons, and 5 unless its
  final payment is discounted to nearly nothing (a 50-year par yield on a 60% curve). Condition
  5 is what rejects payments that offset one another beyond the precision of the valuation, and
  a quote that pays only just after `tₖ₋₁`, where the bump hardly reaches: rounding would
  decide either bump. A violation throws an `ArgumentError` naming the tenor and the condition,
  as does a bump that cannot be bracketed within floating-point limits.

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
        target = __primal(news[k])
        cashflows = __quote_cashflows(q.instrument)
        __check_quote_cashflows(cashflows, t, target)
        payments = __quote_payments(cashflows)
        __check_quote_payments(payments, t_prev, t, target)
        solved = primal + __trial_quote_bumps(pb, ts)
        fixed = sum(a * discount(solved, s) for (s, a, _) in payments if s <= t_prev; init = zero(qp.price))
        qp.price > fixed || throw(
            ArgumentError(
                "quote_shift: the target $target at tenor $t is unreachable: its price $(qp.price) does not exceed " *
                    "$fixed, the value of its payments at or before $t_prev, which its bump cannot move."
            )
        )
        # The residual is the quote's value, by the same valuation that gave the starting quote,
        # less its price, in the quote's unit at each trial bump. A unit fixed before the solve
        # would not do: a large shock moves the value by hundreds of orders of magnitude, and
        # against today's value a far smaller target price rounds to zero (a bump that discounts
        # the quote to nothing would pass for its root) and a far larger one overflows.
        trial(c, b, y) = c + __trial_quote_bumps(__with_bump(b, k, y), ts)
        tail = __last_payment(payments)
        function g_primal(y)
            c = trial(primal, pb, y)
            v = FinanceCore.present_value(c, qp.instrument)
            return (v - qp.price) / __quote_unit(c, v, qp.price, tail)
        end
        # The root is solved inside its bracket. Far from the root this residual is flat, at about
        # ±1 (the value or the price is nothing against the other), so a secant search from zero
        # has no slope to follow under a large shock.
        residual = __PrimalResidual(g_primal, who, hint)
        lo, hi = __quote_bracket(residual, t)
        x = lo == hi ? lo : Roots.find_zero(residual, (lo, hi), Roots.A42())
        # The contract's target: at the root the quote reprices within `sqrt(eps)` of its unit. This
        # is the valuation the caller can repeat, on the curve and the target as given, with their
        # dual parts. Later bumps are zero up to this tenor, where every payment of this quote
        # lies (S1), so this trial curve prices the quote as the result will.
        at_root = trial(curve, bumps, x)
        net = FinanceCore.present_value(at_root, q.instrument)
        unit = __quote_unit(at_root, net, q.price, tail)
        miss = net - q.price
        abs(__primal(miss)) <= sqrt(eps(typeof(unit))) * unit || throw(
            ArgumentError(
                "quote_shift: the quote at tenor $t (target $target) does not reprice: at the bump $x, the nearest the " *
                    "solve can find, its value is $(__primal(net)) against a price of $(qp.price). Its payments may " *
                    "offset one another beyond the precision of the valuation."
            )
        )
        # (S5) A value can also meet its price while its rounding hides what the bump moves: the
        # bump is then whatever the rounding allows. Below the smallest normal number a value's
        # rounding no longer shrinks with it, and sizes that overflowed measure nothing.
        gross, exposure = __quote_sizes(trial(primal, pb, x), payments, t_prev, t, k == 1)
        E = typeof(exposure)
        (isfinite(exposure) && exposure >= sqrt(eps(E)) * max(gross, floatmin(E))) || throw(
            ArgumentError(
                "quote_shift: the quote at tenor $t (target $target) does not determine its bump: its value moves by " *
                    "$exposure per unit of the log discount factor at its tenor, less than `sqrt(eps)` of the $gross " *
                    "that all its payments are worth, each on its own. It pays too little near its tenor, or its " *
                    "payments offset one another, and the bump would be left to the rounding of the valuation."
            )
        )
        # The implicit function theorem at the root: the dual parts of the miss over the value's
        # slope in the bump, `-exposure·t`. Nothing is divided by the unit here, so a partial far
        # larger or smaller than the quote neither overflows nor vanishes.
        depth = __ad_depth(miss)
        depth <= 1 || throw(ArgumentError("$who supports first-order ForwardDiff derivatives only; nested dual numbers are not supported"))
        bumps[k] = depth == 0 ? x : x + (miss - __primal(miss)) / (exposure * t)
    end
    return curve + __quote_bumps(bumps, ts)
end
