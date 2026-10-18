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

# A pool's cashflows per unit of its initial balance: at each payment, the share `F` that has not
# prepaid pays its scheduled payment and prepays `smm` of its remaining balance (see `Bond.Pool`).
@inline function Transducers.__foldl__(rf, val, p::Projection{C, M, K}) where {C <: Bond.Pool, M, K <: CashflowProjection}
    rows = Bond.__schedule(p.contract.loan)
    smms = __prepayment_fractions(p.model[p.contract.key], rows)
    F = one(eltype(smms))
    for (row, smm) in zip(rows, smms)
        val = @next(rf, val, Cashflow(F * (row.interest + row.principal + smm * row.closing), row.time))
        F *= 1 - smm
    end
    return complete(rf, val)
end

# A prescribed prepayment assumption counts payments and their intervals, so it is the same seen from
# any start (`Forward`).
__rebase(m::Prepayment.AbstractPrepayment, t) = m
