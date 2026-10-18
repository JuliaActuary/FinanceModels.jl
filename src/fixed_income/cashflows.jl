# Cashflow projections of amortizing loans and loan pools.

# An amortizing loan pays each scheduled payment's interest and principal together.
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Bond.Amortizing, M, K}
    for row in Bond.__schedule(p.contract)
        val = @next(rf, val, Cashflow(row.interest + row.principal, row.time))
    end
    return complete(rf, val)
end

# The fraction of the remaining balance that prepays at each payment of `rows`, from the assumption
# `pre`: one per payment, each after the interval since the previous payment (from 0).
function __prepayment_fractions(pre, rows)
    Prepayment.__check_length(pre, length(rows))
    return [Prepayment.__smm(pre, i, rows[i].time - (i == 1 ? zero(rows[i].time) : rows[i - 1].time)) for i in eachindex(rows)]
end

# ── An asset's payments per unit outstanding ───────────────────────────────────────────────
#
# Every principal-bearing asset reduces to one row per payment, per unit of the asset alive just
# before it: `opening`, the principal outstanding since the previous payment (recovered on default);
# `payment`, what the payment pays; `continuing`, the fraction that remains after it (below one when a
# pool prepays); and `closing`, the principal of each remaining unit. Defaults, exercise and ALM
# reporting act on these rows.

# `ctx` holds the models the asset reads by key (a pool's prepayment).
function __asset_rows(b::Union{Bond.Fixed, Bond.Amortizing}, ctx)
    return map(Bond.__schedule(b)) do r
        (; r.time, r.opening, payment = r.interest + r.principal, continuing = one(r.opening), r.closing, r.interest, r.principal, prepayment = zero(r.opening))
    end
end
function __asset_rows(p::Bond.Pool, ctx)
    rows = Bond.__schedule(p.loan)
    smms = __prepayment_fractions(ctx[p.key], rows)
    return map(rows, smms) do r, smm
        prepayment = smm * r.closing
        (; r.time, r.opening, payment = r.interest + r.principal + prepayment, continuing = 1 - smm, r.closing, r.interest, r.principal, prepayment)
    end
end

# A pool's cashflows per unit of its initial balance: at each payment, the share that has not
# prepaid pays its scheduled payment and prepays its fraction of the remaining balance (see `Bond.Pool`).
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Bond.Pool, M, K <: CashflowProjection}
    rows = __asset_rows(p.contract, p.model)
    F = one(first(rows).continuing)
    for r in rows
        val = @next(rf, val, Cashflow(F * r.payment, r.time))
        F *= r.continuing
    end
    return complete(rf, val)
end

# A prescribed prepayment assumption counts payments and their intervals, so it is the same seen from
# any start (`Forward`).
__rebase(m::Prepayment.AbstractPrepayment, t) = m
