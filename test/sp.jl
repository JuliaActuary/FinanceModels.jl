using FinanceModels
using FinanceCore
using BenchmarkTools
c = Cashflow(1.0, 1.0)


p = Projection(NullModel(), c, CashflowProjection())

collect(p)

p = Projection(NullModel(), Bond.Fixed(0.05, Periodic(1), 3.0), CashflowProjection())

collect(p)

p = Projection(NullModel(), Forward(5.0, Bond.Fixed(0.05, Periodic(1), 3.0)), CashflowProjection())

collect(p)


FinanceModels.pv(Yield.Constant(0.05), Bond.Fixed(0.05, Periodic(1), 3.0))

y = Yield.Constant(0.05)
Yield.discount(y, 5)
p = Projection(
    Dict("SOFR" => Yield.Constant(0.05)),
    Bond.Floating(0.02, Periodic(1), 3.0, "SOFR"),
    CashflowProjection(),
)

collect(p)



qs = [
    Quote(1.0, Bond.Fixed(0.05, Periodic(1), 3.0)),
    Quote(1.0, Bond.Fixed(0.07, Periodic(1), 3.0)),
]
fit(Yield.Constant, Fit.Loss(x -> x^2), qs)

c = FinanceModels.Composite(Bond.Fixed(0.05, Periodic(1), 3.0), Bond.Fixed(0.1, Periodic(4), 3.0))

p = Projection(NullModel(), CashflowProjection(), c)