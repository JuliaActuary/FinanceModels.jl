# A curve that defines only `discount`, with D(0) ≠ 1: intervals keep D(to)/D(from) semantics.
struct __HalfFlatCurve <: Yield.AbstractYieldModel end
FinanceCore.discount(::__HalfFlatCurve, t) = 0.5 * exp(-0.03 * t)

# A curve with negative discount factors (as a Smith-Wilson fit to arbitrary prices can have): its
# interval factors are still the exact ratio, but it has no log-discount at a single time, only over
# an interval whose factor is positive.
struct __NegativeDiscountCurve <: Yield.AbstractYieldModel end
FinanceCore.discount(::__NegativeDiscountCurve, t) = -exp(-0.03 * t)
