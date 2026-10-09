# Yield Curve Arithmetic

## Curves are deflators, not rates

A yield curve is fundamentally a **price deflator** — it tells you what a future dollar is worth today. When you combine two curves (e.g. a base curve and a spread curve), you're compounding two deflators: a dollar deflated by the base curve and then by the spread curve. Mathematically, discount factors multiply:

```math
D_{\text{combined}}(t) = D_{\text{base}}(t) \cdot D_{\text{spread}}(t)
```

Multiplying discount factors is equivalent to adding rates in **continuous (log) space** — `exp(a) × exp(b) = exp(a + b)` exactly. `+` on curves means "compound the deflators."

## Demonstration

!!! note
    `Yield.Constant` is a subtype of `AbstractYieldModel`, which represents a discount curve associated with a constant force of interest through time.

```julia
using FinanceModels, FinanceCore

base   = Yield.Constant(0.05)   # 5% annual effective
spread = Yield.Constant(0.02)   # 2% annual effective (e.g. credit spread)
combined = base + spread

# The combined curve multiplies discount factors:
t = 10.0
df_base   = discount(base, t)      # 1.05^-10 = 0.61391
df_spread = discount(spread, t)    # 1.02^-10 = 0.82035
df_combined = discount(combined, t) # 0.50362

df_base * df_spread  # 0.50362 — matches combined

# The combined rate is NOT simply 7%:
discount(Yield.Constant(0.07), t)  # 0.50835 — slightly different

# The difference: compounding 5% and 2% as deflators gives
# (1.05 × 1.02) - 1 = 7.1%, not 7%:
1.05 * 1.02  # 1.071
```

### Why the difference?

Adding nominal rates ignores the **cross-term**. If you invest \$1 at the base rate and the spread applies on top:

```math
(1 + r_{\text{base}})(1 + r_{\text{spread}}) = 1 + r_{\text{base}} + r_{\text{spread}} + r_{\text{base}} \, r_{\text{spread}}
```

The cross-term (`0.05 × 0.02 = 0.001`) is what gets dropped when you simply add `0.05 + 0.02 = 0.07`. Over a single year this is 10 bps; over 10 years it compounds.

In continuous space, the cross-term doesn't exist — `exp(a) × exp(b) = exp(a + b)` exactly — so adding continuous zero rates is the correct way to combine deflators. `curve_a + curve_b` creates a [`CompositeYield`](@ref FinanceModels.Yield.CompositeYield) that, at each time `t`:

1. Extracts continuous zero rates from each curve: ``z_i = -\log D_i(t) / t``
2. Adds the zero rates: ``z = z_1 + z_2``
3. Returns the combined DF: ``D(t) = e^{-z \, t}``

Subtraction works analogously — `curve_a - curve_b` divides discount factors, giving the implied spread curve.

## Working with spread curves

Use curve subtraction to find the spread between two curves:

```julia
base   = Yield.Constant(0.05)
target = Yield.Constant(0.07)

# The spread that, when added to base, reproduces target's discount factors:
spread = target - base
zero(spread,1)  # Continuous(0.01887...) — not exactly 2% in nominal terms

# Verify round-trip:
discount(base + spread, 10.0) ≈ discount(target, 10.0)  # true

# You can also add the spread Rate directly:
discount(base + spread.rate, 10.0) ≈ discount(target, 10.0)  # true
```

This pattern — using curve subtraction and adding the result back — ensures you're always working in the correct space, regardless of which rate convention (`Periodic`, `Continuous`, etc.) the inputs use.

## Curves versus rates

The key is that combining `Rate`s is not the same thing as combining curves (`AbstractYieldModel`s):

```julia
a = Yield.Constant(Periodic(0.05,1) + Periodic(0.02,1))
b = Yield.Constant(Periodic(0.05,1)) + Yield.Constant(Periodic(0.02,1))
a != b # true

discount(a, 10)  # 0.50835 — rate addition (drops cross-term)
discount(b, 10)  # 0.50362 — deflator compounding (correct)
```

This is the point made above - combining discount curves results in a different discount path than a curve of combined rates. See the [Migration Guide](@ref) for details on how this changed in v5.

## Scaling curves: `*` and `/`

`curve * α` scales every continuous zero rate by a scalar factor `α`. This is the "power-of-DF" operation:

```math
D_\alpha(t) = e^{-\alpha \, z(t) \, t} = \bigl[D(t)\bigr]^\alpha
```

This creates a [`ScaledYield`](@ref FinanceModels.Yield.ScaledYield).

A common use case is after-tax yields. If the tax rate is 21%, the after-tax curve is `curve * 0.79`:

```julia
pretax = Yield.Constant(Continuous(0.05))
aftertax = pretax * 0.79

discount(aftertax, 10) ≈ exp(-0.05 * 0.79 * 10)  # true
```

Division is the inverse — `curve / α` scales by `1/α`, useful for grossing up to a pre-tax equivalent:

```julia
aftertax = Yield.Constant(Continuous(0.0395))
pretax = aftertax / 0.79

discount(pretax, 10) ≈ exp(-0.0395 / 0.79 * 10)  # true
```

Multiplication is only defined between a curve and a scalar. Multiplying two curves together is not a meaningful operation and will raise a `MethodError`.

## Convenience: curves with scalars and rates

You can add or subtract a scalar or a [`Rate`](@ref) directly to/from a curve. The scalar is wrapped in a [`Yield.Constant`](@ref FinanceModels.Yield.Constant):

```julia
m = fit(Spline.Linear(), ZCBYield.([0.04, 0.05, 0.06], [1, 5, 10]), Fit.Bootstrap())

# These are equivalent:
m + 0.01
m + Yield.Constant(0.01)

# With an explicit rate type:
m + Continuous(0.01)
m + Yield.Constant(Continuous(0.01))
```

!!! note "Default convention for bare scalars"
    A bare number like `0.01` is interpreted as an annual effective rate (`Periodic(1)`) when wrapped in `Constant`. To be explicit about compounding, pass a `Rate` object: `Continuous(0.01)` or `Periodic(0.01, 2)`.

## Blending curves

[`Yield.Blend`](@ref FinanceModels.Yield.Blend) combines curves with weights that sum to 1, in one of three spaces. `Yield.Blend(a, b, w, space)` weights `a` by `w` and `b` by `1 - w`; `Yield.Blend(curves, weights, space)` takes any number of curves. A weight is a number, or a function of tenor `t -> w`.

| Space | The blend | Use it for |
|:------|:----------|:-----------|
| `Yield.DiscountFactors()` | ``D(t) = \sum_i w_i D_i(t)`` | a mixture: present value is linear in the curves, as for an expected discount factor over scenarios |
| `Yield.ZeroRates()` | ``z(t) = \sum_i w_i z_i(t)`` | the same curve as `w * a + (1 - w) * b`, also with a weight that varies with tenor |
| `Yield.ForwardRates(period)` | ``f(t) = \sum_i w_i(kp)\, f_i(t)`` for ``kp \le t < (k+1)p`` | grading one curve's forwards into another's, period by period |

```julia
market = ZeroRateCurve([0.03, 0.035, 0.04], [1.0, 10.0, 20.0])
ultimate = Yield.Constant(Continuous(0.045))

# the market curve's annual forwards, graded into the ultimate rate between 20 and 60 years
grade(t) = clamp((60 - t) / 40, 0, 1)
graded = Yield.Blend(market, ultimate, grade, Yield.ForwardRates(1.0))

# an equally weighted mixture of simulated paths: their expected discount factor
paths = simulate(ShortRate.Vasicek(0.1, 0.03, 0.01, Continuous(0.03)); n_scenarios = 1000, timestep = 0.25, horizon = 30.0)
expected = Yield.Blend(collect(paths), fill(1 / 1000, 1000), Yield.DiscountFactors())
```

The spaces behave differently:

- **Discount factors.** The weights are shares of their sum and must lie in [0, 1]. With constant weights, `pv(blend, cfs) ≈ Σ wᵢ·pv(cᵢ, cfs)` for cashflows valued from time 0, and the long end follows the curve with the lowest rate.
  - Rebased to a later time `s` (`ForwardStarting`), the mixture's weights become each curve's share of the discount factor at `s`, `wᵢ·Dᵢ(s)/D(s)`.
  - An expected discount factor gives the expected present value only of cashflows that don't vary by scenario.
- **Zero rates.** The weights are used as given, so a negative weight extrapolates. With a weight that varies with tenor, the forward rate picks up the term ``w'(t)\,(L_a(t) - L_b(t))``, so it can bump where the weight changes.
- **Forward rates.** Each period's forwards are blended with the weights at the period's start, so the forward is smooth within a period and jumps at its boundaries. The blend is built from the curves' interval factors, so it starts at `discount(blend, 0) == 1`. A weight that varies with tenor costs about `t / period` interval factors of each curve per evaluation at `t`.

Weights that don't sum to 1, to rounding, throw an `ArgumentError`: number weights when the blend is built, and weights that vary with tenor where they are read. With weights that vary with tenor, `discount(blend, Inf)` is `NaN`, since the weights at infinity are not known.

## Operation summary

| Expression | Result type | Semantics |
|:-----------|:------------|:----------|
| `a + b` | `CompositeYield` | Compound deflators: ``D_a \cdot D_b`` |
| `a - b` | `CompositeYield` | Spread between curves: ``D_a / D_b`` |
| `curve * α` | `ScaledYield` | ``D^\alpha`` (CZR scaling) |
| `curve / α` | `ScaledYield` | ``D^{1/\alpha}`` (CZR scaling) |
| `Yield.Blend(a, b, w, Yield.DiscountFactors())` | `Blend` | ``w D_a + (1-w) D_b`` |
| `Yield.Blend(a, b, w, Yield.ZeroRates())` | `Blend` | ``z = w z_a + (1-w) z_b`` |
| `Yield.Blend(a, b, w, Yield.ForwardRates(p))` | `Blend` | ``f = w f_a + (1-w) f_b``, `w` read once per period |

!!! warning "Par-rate spreads cannot be composed additively"
    Curve arithmetic operates on zero rates. If your base rates and spreads are quoted as **par rates**, you cannot simply add the spread curve to the base curve and get the same result as fitting a single curve to the combined par rates. Par rates depend on the path of rates at earlier tenors, so they must be converted to zero rates (e.g. via bootstrap) before composition. See the examples in [`CompositeYield`](@ref FinanceModels.Yield.CompositeYield) for a demonstration of this difference.

## Performance

[`CompositeYield`](@ref FinanceModels.Yield.CompositeYield) and [`ScaledYield`](@ref FinanceModels.Yield.ScaledYield) evaluate both underlying curves on every call to `discount`. If you are using a composite curve in a hot loop, consider pre-computing the combined zero rates at the tenors you need and fitting a single spline curve to the result.
