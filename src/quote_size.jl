# Whether a contract's cashflows are fixed: the same whatever the curve, so that its price is linear in the
# curve's discount factors. Every quote type FinanceModels defines is a `Cashflow` or a `Bond.Fixed`.
__fixed_cashflows(c) = false
__fixed_cashflows(::FinanceCore.Cashflow) = true
__fixed_cashflows(::Bond.Fixed) = true
__fixed_cashflows(c::FinanceCore.Composite) = __fixed_cashflows(c.a) && __fixed_cashflows(c.b)
__fixed_cashflows(v::Union{AbstractVector, Tuple}) = all(__fixed_cashflows, v)

# Quote residuals are measured against gross discounted projected cashflows, including
# curve-dependent payments. A known value avoids a second projection when every amount
# is nonnegative. Pricing-only custom contracts retain their legacy price/value scale:
# their callback does not disclose gross payments. This fallback is kept here so the
# released bootstrap and implied-quote APIs do not acquire a projection requirement.
__nonnegative_cashflows(c) = false
__nonnegative_cashflows(c::FinanceCore.Cashflow) = c.amount >= 0
__nonnegative_cashflows(b::Bond.Fixed) = b.coupon_rate >= 0
__nonnegative_cashflows(c::FinanceCore.Composite) = __nonnegative_cashflows(c.a) && __nonnegative_cashflows(c.b)
__nonnegative_cashflows(v::Union{AbstractVector, Tuple}) = all(__nonnegative_cashflows, v)

# Cashflow families whose projection is known to FinanceModels. Unknown custom contracts
# keep the pricing-only contract; fixed-cashflow extensions can declare __fixed_cashflows.
__projects_cashflows(c) = __fixed_cashflows(c)
__projects_cashflows(::Bond.Floating) = true
__projects_cashflows(c::FinanceCore.Composite) = __projects_cashflows(c.a) && __projects_cashflows(c.b)
__projects_cashflows(c::Forward) = __projects_cashflows(c.instrument)
__projects_cashflows(v::Union{AbstractVector, Tuple}) = all(__projects_cashflows, v)
__projects_cashflows(e::Transducers.Eduction{<:Any, <:FinanceCore.AbstractContract}) = __projects_cashflows(e.coll)

function __quote_size(q, curve, value)
    __projects_cashflows(q.instrument) || return float(max(abs(value), abs(q.price)))
    __nonnegative_cashflows(q.instrument) && return float(abs(value))
    payments = Projection(q.instrument, curve, CashflowProjection())
    return Transducers.foldxl(
        +, payments |> Transducers.Map(cf -> abs(cf.amount) * discount(curve, cf.time));
        init = zero(float(q.price))
    )
end
