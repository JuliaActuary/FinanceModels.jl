# Quote residuals are measured against gross discounted projected cashflows, including
# curve-dependent payments. A known value avoids a second projection when every amount
# is nonnegative. Pricing-only custom contracts retain their legacy price/value scale:
# their callback does not disclose gross payments. This fallback is kept here so the
# released bootstrap and implied-quote APIs do not acquire a projection requirement.
__nonnegative_cashflows(c) = false
__nonnegative_cashflows(c::FinanceCore.Cashflow) = c.amount >= 0
__nonnegative_cashflows(b::FinanceModels.Bond.Fixed) = b.coupon_rate >= 0
__nonnegative_cashflows(c::FinanceCore.Composite) = __nonnegative_cashflows(c.a) && __nonnegative_cashflows(c.b)
__nonnegative_cashflows(v::Union{AbstractVector, Tuple}) = all(__nonnegative_cashflows, v)

# Cashflow families whose projection is known to FinanceModels. Unknown custom contracts
# keep the pricing-only contract; fixed-cashflow extensions can declare __fixed_cashflows.
__projects_cashflows(c) = FinanceModels.__fixed_cashflows(c)
__projects_cashflows(::FinanceModels.Bond.Floating) = true
__projects_cashflows(c::FinanceCore.Composite) = __projects_cashflows(c.a) && __projects_cashflows(c.b)
__projects_cashflows(c::FinanceModels.Forward) = __projects_cashflows(c.instrument)
__projects_cashflows(v::Union{AbstractVector, Tuple}) = all(__projects_cashflows, v)

function __quote_size(q, curve; value = nothing)
    if !__projects_cashflows(q.instrument)
        v = value === nothing ? FinanceCore.present_value(curve, q.instrument) : value
        return float(max(abs(v), abs(q.price)))
    end
    value !== nothing && __nonnegative_cashflows(q.instrument) && return float(abs(value))
    payments = FinanceModels.Projection(q.instrument, curve, FinanceModels.CashflowProjection())
    return Transducers.foldxl(
        +, payments |> Transducers.Map(cf -> abs(cf.amount) * discount(curve, cf.time));
        init = zero(float(q.price))
    )
end
