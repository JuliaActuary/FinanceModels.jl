# Expected cashflows and values of defaultable assets, and the calibration of hazard curves.

# ── Recovery quadrature ───────────────────────────────────────────────────────────────────
#
# Recovery accrues continuously between payments, at the rate R·B·λ(u)·S(u)·D(u). Each payment
# interval is split at the hazard curve's breakpoints, each piece into equal parts no longer than a
# month, and each part integrated by 4-point Gauss–Legendre. The nodes depend on primal times only, so
# a dual run takes the same nodes as the primal run, and each node is an expected recovery cashflow:
# discounting them is the quadrature.
const __GL4_NODES = (-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
const __GL4_WEIGHTS = (0.34785484513745385, 0.6521451548625461, 0.6521451548625461, 0.34785484513745385)
const __RECOVERY_STEP = 1 / 12

# The cut points of the interval [a, b]: its ends and the breakpoints strictly inside, in primal time.
function __recovery_cuts(a, b, breaks)
    cuts = [a]
    for t in breaks
        __primal(a) < __primal(t) < __primal(b) && push!(cuts, t)
    end
    push!(cuts, b)
    return cuts
end

# The quadrature nodes `u` and weights `w` of the interval [a, b], as `(u, w)` pairs.
function __recovery_nodes(a, b, breaks)
    cuts = __recovery_cuts(a, b, breaks)
    nodes = Tuple{typeof(float(a + b)), typeof(float(b - a))}[]
    for k in 1:(length(cuts) - 1)
        lo, hi = cuts[k], cuts[k + 1]
        m = max(1, ceil(Int, __primal(hi - lo) / __RECOVERY_STEP - 1.0e-9))
        h = (hi - lo) / m
        for j in 1:m
            left = lo + (j - 1) * h
            for (x, w) in zip(__GL4_NODES, __GL4_WEIGHTS)
                push!(nodes, (left + h * (1 + x) / 2, h * w / 2))
            end
        end
    end
    return nodes
end

# ── Expected cashflows of a defaultable asset ─────────────────────────────────────────────

# The expected cashflows of `Credit.Defaultable(asset, key)` under the context's assumptions: each
# payment weighted by survival (and by the pool's continuing share), and, between payments, the
# recovery of the principal outstanding at default, at the quadrature nodes. No payment follows a
# default.
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Credit.Defaultable, M, K <: CashflowProjection}
    d = p.contract
    a = p.model[d.key]
    rows = __asset_rows(d.asset, p.model)
    h, R = a.hazard, a.recovery
    breaks = Credit.__breakpoints(h)
    # the mass alive just before the current payment: survival times the pool's continuing share
    alive = one(R * Credit.survival(h, zero(first(rows).time)))
    t_prev = zero(first(rows).time)
    for r in rows
        base = alive * R * r.opening
        for (u, w) in __recovery_nodes(t_prev, r.time, breaks)
            val = @next(rf, val, Cashflow(base * w * Credit.hazard(h, u) * Credit.survival(h, t_prev, u), u))
        end
        s = Credit.survival(h, t_prev, r.time)
        val = @next(rf, val, Cashflow(alive * s * r.payment, r.time))
        alive *= s * r.continuing
        t_prev = r.time
    end
    return complete(rf, val)
end

# ── Calibrating a hazard curve ────────────────────────────────────────────────────────────

# Primal copies of the contracts a hazard calibration reads (see `__primal_contract`).
function __primal_contract(d::Credit.Defaultable)
    a = __primal_contract(d.asset)
    return a === d.asset ? d : Credit.Defaultable(a, d.key)
end
function __primal_contract(b::Bond.Amortizing)
    foreach(__no_dual_time, b.times)
    all(x -> __ad_depth(x) == 0, (b.rate, first(b.principal), first(b.accruals))) && return b
    return Bond.Amortizing(__primal(b.rate), b.times, map(__primal, b.principal), map(__primal, b.accruals))
end
function __primal_contract(p::Bond.Pool)
    l = __primal_contract(p.loan)
    return l === p.loan ? p : Bond.Pool(l, p.key)
end

const __HAZARD_FIT_HINT = "Hazard calibration differentiates quote prices, the coupons and amounts of the bonds, the " *
    "recovery fraction and the discount curve's parameters; a dual number elsewhere is not supported."

"""
    fit(Credit.Assumptions(Credit.PiecewiseHazard(); recovery), quotes, Fit.Bootstrap(); discount)

Bootstrap an issuer's hazard curve from the prices of its bonds. Each quote is a
`Quote(price, Credit.Defaultable(bond, key))` of a nonoptional bond (a `Bond.Fixed`,
`Bond.Amortizing` or `Bond.Pool`) of one issuer `key`, with its full (dirty) price in the bond's
principal units. The returned [`Credit.Assumptions`](@ref FinanceModels.Credit.Assumptions) has a
[`Credit.PiecewiseHazard`](@ref FinanceModels.Credit.PiecewiseHazard) with one segment ending at each
quote maturity and the template's `recovery`; `discount`, a yield curve, discounts and is held fixed,
as is the recovery.

The segments are solved shortest maturity first, each to reprice its quote: a bond depends only on
the intensities up to its maturity. A quote whose price needs a negative intensity on its segment
(a price above what default-free survival over that segment would give) throws an
`ArgumentError`, as do two quotes with one maturity and quotes of different issuers. Prices that
fall as the intensity rises are the usual case; for a deep-discount bond with a high recovery
assumption the price can rise with the intensity, and the solved segment is then the root the
solver reaches from the previous segment's intensity.

Derivatives of the result, first and second order, are exact with respect to the quote prices,
the bonds' coupons and amounts, the recovery and the parameters of `discount`, by the implicit
function theorem at each segment.

```julia
bonds = [Credit.Defaultable(Bond.Fixed(0.05, Periodic(2), t), "ACME") for t in (2, 5, 10)]
quotes = Quote.([0.99, 0.97, 0.93], bonds)
acme = fit(Credit.Assumptions(Credit.PiecewiseHazard(); recovery = 0.4), quotes, Fit.Bootstrap(); discount = curve)
present_value(Models(curve, Dict("ACME" => acme)), bonds[2])   # ≈ 0.97
```
"""
function fit(template::Credit.Assumptions{<:Credit.PiecewiseHazard}, quotes, ::Fit.Bootstrap; discount::Yield.AbstractYieldModel)
    qs, pqs = __calibration_quotes(quotes)
    keys = unique(q.instrument.key for q in pqs)
    length(keys) == 1 ||
        throw(ArgumentError("a hazard curve is calibrated to one issuer's bonds; the quotes are of $(length(keys)) keys: $keys"))
    key = only(keys)
    order = sortperm(pqs; by = maturity)
    qs, pqs = qs[order], pqs[order]
    tenors = [float(maturity(q)) for q in pqs]
    for i in 2:length(tenors)
        tenors[i] > tenors[i - 1] || throw(
            ArgumentError("two quotes mature at $(tenors[i]); a hazard bootstrap solves one segment per maturity")
        )
    end
    R, pR = template.recovery, __primal(template.recovery)
    primal_discount = Yield.__PrimalCurve(discount)
    trial(λs) = Credit.PiecewiseHazard(Yield.Unchecked(), λs, tenors[1:length(λs)])
    context(curve, λs, recovery) = Models(curve, Dict(key => Credit.Assumptions(trial(λs), recovery)))
    # The intensities hold every partial of every residual: the curve's, the recovery's and each quote's.
    B = mapreduce(q -> typeof(present_value(context(discount, zeros(length(tenors)), R), q.instrument) - q.price), promote_type, qs)
    λs = zeros(B, length(tenors))
    who = "fit(Credit.PiecewiseHazard(), Bootstrap)"
    for i in eachindex(tenors)
        q, qp = qs[i], pqs[i]
        pλ = map(__primal, λs[1:(i - 1)])
        seed = i == 1 ? 0.0 : last(pλ)
        value(x) = present_value(context(primal_discount, vcat(pλ, x), pR), qp.instrument)
        # the residual relative to the quote's size, as a yield bootstrap does
        s = max(abs(qp.price), abs(value(seed)))
        σ = iszero(s) ? one(s) : s
        g(x) = (present_value(context(discount, vcat(λs[1:(i - 1)], x), R), q.instrument) - q.price) / σ
        g_primal(x) = (value(x) - qp.price) / σ
        λ = __implicit_root(g, g_primal, seed; bracket = (-1.0, 5.0), who, hint = __HAZARD_FIT_HINT, scale = x -> one(σ), max_order = 2)
        __primal(λ) >= 0 || throw(
            ArgumentError(
                "the quote of $(key) maturing at $(tenors[i]) needs a negative default intensity ($(__primal(λ))) on its " *
                    "segment: its price $(qp.price) is above what survival without default risk there would give"
            )
        )
        λs[i] = λ
    end
    return Credit.Assumptions(Credit.PiecewiseHazard(λs, tenors), R)
end
