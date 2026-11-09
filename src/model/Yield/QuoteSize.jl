# Quote residuals are measured against gross discounted projected cashflows, including
# curve-dependent payments. A known value avoids a second projection when every amount
# is nonnegative; the definition is the same in either evaluation path.
__nonnegative_cashflows(c) = false
__nonnegative_cashflows(c::FinanceCore.Cashflow) = c.amount >= 0
__nonnegative_cashflows(b::FinanceModels.Bond.Fixed) = b.coupon_rate >= 0
__nonnegative_cashflows(c::FinanceCore.Composite) = __nonnegative_cashflows(c.a) && __nonnegative_cashflows(c.b)
__nonnegative_cashflows(v::Union{AbstractVector, Tuple}) = all(__nonnegative_cashflows, v)

function __quote_size(q, curve; value = nothing)
    value !== nothing && __nonnegative_cashflows(q.instrument) && return float(abs(value))
    payments = FinanceModels.Projection(q.instrument, curve, FinanceModels.CashflowProjection())
    return Transducers.foldxl(
        +, payments |> Transducers.Map(cf -> abs(cf.amount) * discount(curve, cf.time));
        init = zero(float(q.price))
    )
end
