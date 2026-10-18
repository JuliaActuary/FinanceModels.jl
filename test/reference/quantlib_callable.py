"""Reference values for test/lattice.jl: callable and putable fixed-rate bonds under Hull–White.

QuantLib's TreeCallableFixedRateBondEngine on a flat continuously compounded curve. Every date
(curve, coupons, calls) uses 30/360 bond basis, so coupon and call dates fall at exact half-years,
matching FinanceModels' `Bond.Fixed(c, Periodic(2), 10)` on `Yield.Constant(Continuous(r))`. QuantLib
applies a call to the continuation value and then adds the coupon, as FinanceModels does.

Run with QuantLib 1.44: `python quantlib_callable.py`. Prices are per unit of face value.
"""

import QuantLib as ql

today = ql.Date(15, ql.January, 2026)
ql.Settings.instance().evaluationDate = today
basis = ql.Thirty360(ql.Thirty360.BondBasis)


def price(rate, a, sigma, coupon, first_exercise, kind, steps):
    curve = ql.YieldTermStructureHandle(ql.FlatForward(today, rate, basis, ql.Continuous))
    maturity = today + ql.Period(10, ql.Years)
    schedule = ql.Schedule(today, maturity, ql.Period(ql.Semiannual), ql.NullCalendar(),
                           ql.Unadjusted, ql.Unadjusted, ql.DateGeneration.Backward, False)
    exercises = ql.CallabilitySchedule()
    for d in schedule:
        if today + ql.Period(int(first_exercise * 12), ql.Months) <= d < maturity:
            exercises.append(ql.Callability(ql.BondPrice(100.0, ql.BondPrice.Clean), kind, d))
    bond = ql.CallableFixedRateBond(0, 100.0, schedule, [coupon], basis, ql.Unadjusted, 100.0, today, exercises)
    bond.setPricingEngine(ql.TreeCallableFixedRateBondEngine(ql.HullWhite(curve, a, sigma), steps))
    return bond.NPV() / 100.0


cases = [
    (0.04, 0.10, 0.010, 0.05, 3.0, ql.Callability.Call),
    (0.04, 0.10, 0.010, 0.05, 3.0, ql.Callability.Put),
    (0.02, 0.03, 0.015, 0.03, 2.0, ql.Callability.Call),
]
for rate, a, sigma, coupon, first, kind in cases:
    name = "call" if kind == ql.Callability.Call else "put"
    values = [price(rate, a, sigma, coupon, first, kind, n) for n in (1000, 2000, 4000)]
    print(f"r={rate} a={a} sigma={sigma} coupon={coupon} first={first} {name}: " +
          " ".join(f"{v:.12f}" for v in values))
