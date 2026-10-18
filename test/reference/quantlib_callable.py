"""Reference values for test/lattice.jl and test/published_examples.jl: QuantLib's
TreeCallableFixedRateBondEngine under Hull–White.

QuantLib applies a call or put to the continuation value and then adds that date's coupon, as
FinanceModels does, and both lattices branch to the nearest nodes without edge truncation, so equal
time grids give equal values.

1. Callable and putable fixed-rate bonds on a flat continuously compounded curve (test/lattice.jl).
   Every date (curve, coupons, calls) uses 30/360 bond basis, so coupon and call dates fall at exact
   half-years, matching `Bond.Fixed(c, Periodic(2), 10)` on `Yield.Constant(Continuous(r))`.
   QuantLib's N steps over ten years are FinanceModels' `max_step = 10 / N`.
2. Hull and White (1996), "Using Hull–White interest rate trees", Journal of Derivatives 3(3): a
   3-year European put, strike 63, on a zero paying 100 at 9 years, on the paper's DM zero curve
   (Table 1: continuously compounded, Actual/365, linear between maturities), a = 0.1, σ = 0.01
   (test/published_examples.jl). The put is a zero with one put date less the zero. QuantLib's 3N
   uniform steps to 9 years are FinanceModels' `max_step = 3 / N`.
3. QuantLib's own test case, test-suite/callablebonds.cpp `testCached` (test/lattice.jl), on
   FinanceModels' grid for `max_step = 1/24` (each interval between payment dates split into
   ceil(length / max_step) equal steps), passed as an explicit TimeGrid. On its own 240-step grid
   QuantLib caches clean prices 110.60975477 (callable) and 115.16559362 (putable), which this script
   also prints.

Run with QuantLib 1.44: `python quantlib_callable.py`. Prices are per 100 of face value.
"""

import math

import QuantLib as ql

# 1. Callable and putable fixed-rate bonds

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
    return bond.NPV()


cases = [
    (0.04, 0.10, 0.010, 0.05, 3.0, ql.Callability.Call),
    (0.04, 0.10, 0.010, 0.05, 3.0, ql.Callability.Put),
    (0.02, 0.03, 0.015, 0.03, 2.0, ql.Callability.Call),
]
for rate, a, sigma, coupon, first, kind in cases:
    name = "call" if kind == ql.Callability.Call else "put"
    values = [price(rate, a, sigma, coupon, first, kind, n) / 100 for n in (1000, 2000, 4000)]
    print(f"r={rate} a={a} sigma={sigma} coupon={coupon} first={first} {name}: " +
          " ".join(f"{v:.12f}" for v in values))

# 2. Hull and White (1996): a 3-year put on a 9-year zero

today = ql.Date(8, ql.July, 1994)
ql.Settings.instance().evaluationDate = today
dc = ql.Actual365Fixed()
days = [0, 3, 31, 62, 94, 185, 367, 731, 1096, 1461, 1826, 2194, 2558, 2922, 3287, 3653]
rates = [5.01772, 5.01772, 4.98284, 4.97234, 4.96157, 4.99058, 5.09389, 5.79733, 6.30595, 6.73464,
         6.94816, 7.08807, 7.27527, 7.30852, 7.39790, 7.49015]  # flat before 3 days
curve = ql.ZeroCurve([today + d for d in days], [r / 100 for r in rates], dc, ql.NullCalendar(),
                     ql.Linear(), ql.Continuous, ql.Annual)
ts = ql.YieldTermStructureHandle(curve)
hw = ql.HullWhite(ts, 0.1, 0.01)
print("Hull-White 1996 analytic put: %.10f" % (100 * hw.discountBondOption(ql.Option.Put, 0.63, 3.0, 9.0)))
schedule = ql.Schedule(today, today + 3285, ql.Period(ql.Once), ql.NullCalendar(), ql.Unadjusted,
                       ql.Unadjusted, ql.DateGeneration.Backward, False)
puts = ql.CallabilitySchedule()
puts.append(ql.Callability(ql.BondPrice(63.0, ql.BondPrice.Clean), ql.Callability.Put, today + 1095))
bond = ql.CallableFixedRateBond(0, 100.0, schedule, [0.0], dc, ql.Unadjusted, 100.0, today, puts)
for n in (10, 30, 50, 100, 200):
    bond.setPricingEngine(ql.TreeCallableFixedRateBondEngine(hw, 3 * n, ts))
    print("Hull-White 1996 tree put, %d steps to expiry: %.12f" % (n, bond.NPV() - 100 * ts.discount(9.0)))

# 3. QuantLib's test-suite callable bond (testCached)

today = ql.Date(3, ql.June, 2004)
ql.Settings.instance().evaluationDate = today
calendar = ql.TARGET()
settlement = calendar.advance(today, 3, ql.Days)
ts = ql.YieldTermStructureHandle(ql.FlatForward(settlement, 0.032, dc))
model = ql.HullWhite(ts)  # a = 0.1, σ = 0.01
issue = calendar.adjust(today - 100)
maturity = calendar.advance(issue, 10, ql.Years)
schedule = ql.MakeSchedule(issue, maturity, ql.Period(ql.Semiannual), calendar=calendar,
                           convention=ql.ModifiedFollowing, rule=ql.DateGeneration.Backward)
call_dates = [calendar.advance(issue, i, ql.Years) for i in range(2, 10, 2)]
put_dates = [calendar.advance(issue, i, ql.Years) for i in range(1, 10, 2)]


def exercises(dates, price, kind):
    s = ql.CallabilitySchedule()
    for d in dates:
        s.append(ql.Callability(ql.BondPrice(price, ql.BondPrice.Clean), kind, d))
    return s


def bond_with(exercise):
    return ql.CallableFixedRateBond(3, 10000.0, schedule, [0.05], basis, ql.ModifiedFollowing, 100.0, issue, exercise)


straight = ql.FixedRateBond(3, 100.0, schedule, [0.05], basis, ql.ModifiedFollowing, 100.0, issue)
flows = [cf for cf in straight.cashflows() if cf.date() > settlement]
print("testCached days from settlement:", [cf.date() - settlement for cf in flows[:-1]])
print("testCached 30/360 coupon days:", [round(cf.amount() / 5 * 360) for cf in flows[:-1]])
mandatory = sorted(set(dc.yearFraction(settlement, cf.date()) for cf in flows))
max_step = 1 / 24
grid, previous = [0.0], 0.0
for t in mandatory:
    n = max(1, math.ceil((t - previous) / max_step - 1e-9))
    grid += [previous + k * (t - previous) / n for k in range(1, n)] + [t]
    previous = t
for name, exercise in [("callable", exercises(call_dates, 110.0, ql.Callability.Call)),
                       ("putable", exercises(put_dates, 100.0, ql.Callability.Put))]:
    b = bond_with(exercise)
    b.setPricingEngine(ql.TreeCallableFixedRateBondEngine(model, 240, ts))
    cached = b.cleanPrice()
    b.setPricingEngine(ql.TreeCallableFixedRateBondEngine(model, ql.TimeGrid(grid), ts))
    print("testCached %s: cached clean %.8f; dirty on the max_step = 1/24 grid (%d steps) %.12f" %
          (name, cached, len(grid) - 1, b.dirtyPrice()))
