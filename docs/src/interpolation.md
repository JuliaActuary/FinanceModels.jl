# Interpolation Methods for `ZeroRateCurve`

`ZeroRateCurve` accepts an optional third argument specifying the interpolation method. The choice of interpolation affects **forward curve smoothness**, **key rate duration locality**, and **performance** when used with automatic differentiation (e.g. via `sensitivities()` in ActuaryUtilities.jl).

## Available Methods

| Method | Interior smoothness | Locality | Description |
|--------|-----------|----------|-------------|
| `Spline.MonotoneConvex()` | C1 (smooth) | Best among smooth | **Default.** Finance-aware. Positive forwards, best KRD locality, fastest AD. |
| `Spline.PCHIP()` | C1 (smooth) | Local | Monotonicity-preserving, local. Good general-purpose alternative. |
| `Spline.Akima()` | C1 (smooth) | Local | Local, resistant to outlier oscillation. |
| `Spline.Linear()` | C0 (kinked) | Perfectly local | Simplest. Kinks in forward curve at tenor points. |
| `Spline.Cubic()` | C2 (smoothest) | Global | Natural cubic spline. Smoothest, but bumping one rate affects the whole curve. Thread-safe. |
| `Spline.BSpline(n)` | Varies | Global | nth-order B-spline. Opt-in basis for least-squares *fitting*; **not thread-safe** for concurrent evaluation on a shared curve (shared internal buffer). |

Since v6.0, `Spline.Linear()`, `Spline.Quadratic()`, and `Spline.Cubic()` use piecewise polynomial interpolants (`PolynomialSpline`) that are safe to evaluate concurrently across threads. Evaluating a local polynomial piece does not imply local knot-rate sensitivities: cubic coefficients depend on the whole grid. For a *global* B-spline (e.g. as a least-squares fitting basis) request `Spline.BSpline(d)` explicitly; it is **not** thread-safe for concurrent evaluation on a shared curve.

```julia
using FinanceModels

rates = [0.02, 0.03, 0.035, 0.04, 0.045]
tenors = [1.0, 2.0, 5.0, 10.0, 20.0]

zrc = ZeroRateCurve(rates, tenors)                              # default: MonotoneConvex
zrc_pchip = ZeroRateCurve(rates, tenors, Spline.PCHIP())        # PCHIP
zrc_lin = ZeroRateCurve(rates, tenors, Spline.Linear())          # linear
zrc_cub = ZeroRateCurve(rates, tenors, Spline.Cubic())           # natural cubic spline
zrc_aki = ZeroRateCurve(rates, tenors, Spline.Akima())           # Akima
zrc_flat_zero = ZeroRateCurve(rates, tenors, Spline.Cubic();
    extrapolation=:flat_zero)
```

## Extrapolation Beyond the Last Knot

Every knot-based curve defaults to the same `extrapolation=:flat_forward` rule.
The resulting terminal forward level depends on the interpolation method. For a final knot `(tₙ, zₙ)`, let
`fₙ = zₙ + tₙz′(tₙ⁻)` be the left-hand instantaneous forward there. Then for `t > tₙ`,

```math
z(t) = f_n + (z_n - f_n)\frac{t_n}{t}.
```

This preserves the discount factor at `tₙ`, keeps the instantaneous forward continuous (though
not necessarily differentiable) at the boundary, and prevents a quadratic or cubic's final
polynomial piece from dominating at long horizons. The extrapolated zero rate is not flat
unless `fₙ = zₙ`; it converges to `fₙ` as the horizon increases. This policy applies to
`Yield.Spline`, `ZeroRateCurve`, and `Yield.MonotoneConvex`. A callable-only `Yield.Spline(fn)`
has no knots and therefore does not apply an extrapolation policy.

The terminal forward is a consequence of the interpolant's endpoint derivative, rather
than an independently selected long-term rate. For example, zero rates
`[0.02, 0.025, 0.03, 0.035, 0.04, 0.042]` at `[1, 2, 5, 10, 20, 30]` produce terminal
forwards of about 4.35% with PCHIP, 4.76% with Cubic, and 6.57% with `BSpline(3)`.
The default prevents polynomial continuation from driving the tail; it does not constrain
the terminal forward's level or sign for arbitrary interpolants and inputs.

These are general curve-extension choices, not implementations of an accounting basis
or prescribed regulatory extrapolation. For convergence to an independently chosen
ultimate forward rate, see [`Yield.SmithWilson`](@ref), which accepts `ufr` and a
convergence-speed parameter `α`; model selection and calibration remain explicit choices.

Pass the `extrapolation` keyword to `ZeroRateCurve`, `Yield.Spline`, or a spline `fit` call
to select the long-end behavior:

| Value | Long-end zero rate | Notes |
|-------|--------------------|-------|
| `:flat_forward` | `fₙ + (zₙ-fₙ)tₙ/t` | **Default.** Preserves forward continuity. |
| `Yield.FlatForwardAt(f)` | `f + (zₙ-f)tₙ/t` | Independent, continuously-compounded forward assumption; usually introduces a boundary jump. |
| `:flat_zero` | `zₙ` | Holds the zero rate constant; usually introduces a forward jump at `tₙ`. |
| `:linear` | `zₙ + z′(tₙ⁻)(t-tₙ)` | Extends the zero rate at its boundary slope. |
| `:extension` | Final interpolation piece | Restores the former DataInterpolations behavior and can become extreme far beyond the grid. Not available for `Spline.MonotoneConvex()`. |

```julia
curve = Yield.Spline(Spline.Cubic(), tenors, rates; extrapolation=:flat_zero)
zrc = ZeroRateCurve(rates, tenors, Spline.Cubic(); extrapolation=:linear)
quotes = ZCBYield.(Continuous.(rates), tenors)
fitted = fit(Spline.Linear(), quotes, Fit.Bootstrap(); extrapolation=:extension)
```

The selected value is available as `curve.extrapolation` on `ZeroRateCurve` and native
`Yield.MonotoneConvex` curves. It is preserved by fitting and `Accessors.@set`; derived
boundary quantities are recomputed when knots change. Fitting `Spline.MonotoneConvex()`
returns a native `Yield.MonotoneConvex` with `.rates`, `.times`, and
`Yield.instantaneous_forward` for every supported policy.

An explicit forward can be chosen directly or computed from the last discrete forward:

```julia
assumed = ZeroRateCurve(rates, tenors; extrapolation=Yield.FlatForwardAt(0.035))
n = length(tenors)
last_discrete = (rates[n]*tenors[n] - rates[n-1]*tenors[n-1]) / (tenors[n] - tenors[n-1])
discrete_tail = ZeroRateCurve(rates, tenors;
    extrapolation=Yield.FlatForwardAt(last_discrete))
```

`FlatForwardAt` treats a bare number as continuously compounded (unlike
`Yield.Constant`, which treats a bare number as annual effective). It also accepts rate
objects: `Yield.FlatForwardAt(Periodic(0.035, 1))` converts an annual-effective 3.5% rate,
while `Yield.FlatForwardAt(Continuous(0.035))` is equivalent to the bare-number example.
The normalized continuous value is stored in the policy and stays fixed when knot rates
are bumped or fitted. Recompute `last_discrete` and reconstruct the policy if that anchor
should instead track changed knots. Appending a synthetic far knot does not generally
impose a chosen terminal instantaneous forward and can also change interior interpolation.

## Key Tradeoffs

### Forward Curve Smoothness

The instantaneous forward rate `f(t) = r(t) + t · r'(t)` should be smooth for stochastic models that differentiate the forward curve (e.g. Hull-White θ(t) calibration). Linear interpolation creates discontinuous jumps in `f(t)` at tenor points, while PCHIP, MonotoneConvex, Akima, and CubicSpline produce smooth forward curves.

The following code evaluates forward rates near the 2yr and 5yr tenor points to illustrate the difference:

```julia
using FinanceModels
using FinanceCore: discount

rates = [0.02, 0.025, 0.03, 0.035, 0.04]
tenors = [1.0, 2.0, 5.0, 10.0, 20.0]

eval_points = [1.9, 1.99, 2.0, 2.01, 2.1, 4.9, 4.99, 5.0, 5.01, 5.1]

# Numerical instantaneous forward from the curve's continuous zero rate.
function fwd_from_zero(model, t)
    h = 1e-6
    z = rate(zero(model, t))
    dz = (rate(zero(model, t+h)) - rate(zero(model, t-h))) / (2h)
    z + t * dz
end

for (name, descriptor) in [
    ("Linear", Spline.Linear()),
    ("PCHIP", Spline.PCHIP()),
    ("MonotoneConvex", Spline.MonotoneConvex()),
    ("Akima", Spline.Akima()),
    ("CubicSpline", Spline.Cubic()),
]
    model = Yield.build_model(descriptor, tenors, rates)
    fwd(t) = fwd_from_zero(model, t)
    println("\n--- $name: forward rate f(t) ---")
    for t in eval_points
        println("  t=$(lpad(round(t, digits=2), 5)):  f=$(round(fwd(t)*100, digits=4))%")
    end
end
```

Results:

| Method | f(1.99) | f(2.0) | f(2.01) | Jump? | f(4.99) | f(5.0) | f(5.01) | Jump? |
|--------|---------|--------|---------|-------|---------|--------|---------|-------|
| **Linear** | 3.49% | 3.17% | 2.84% | **Yes** | 3.83% | 3.67% | 3.50% | **Yes** |
| **PCHIP** | 3.05% | 3.05% | 3.05% | No | 3.63% | 3.64% | 3.64% | No |
| **MonotoneConvex** | 3.08% | 3.08% | 3.09% | No | 3.58% | 3.58% | 3.59% | No |
| **Akima** | 2.96% | 2.94% | 2.95% | No | 3.54% | 3.54% | 3.55% | No |
| **CubicSpline** | 3.32% | 3.33% | 3.33% | No | 3.32% | 3.32% | 3.33% | No |

MonotoneConvex additionally guarantees positive continuous forward rates when input rates imply positive forwards — a property unique to this method among those listed ([Hagan & West, 2006](https://doi.org/10.1080/13504860600829233)).

### Key Rate Duration Locality

When computing key rate durations (KRDs), bumping one zero rate should ideally affect only nearby discount factors. The table below shows `∂rate(t)/∂r₃` — the sensitivity of the interpolated rate at various times to a bump in the 5yr rate (rate index 3, with tenors at 1, 2, 5, 10, 20):

```julia
using FinanceModels
using FinanceCore: discount
using ForwardDiff

rates = [0.02, 0.03, 0.035, 0.04, 0.045]
tenors = [1.0, 2.0, 5.0, 10.0, 20.0]

eval_points = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]

for (name, descriptor) in [
    ("Linear", Spline.Linear()),
    ("PCHIP", Spline.PCHIP()),
    ("MonotoneConvex", Spline.MonotoneConvex()),
    ("Akima", Spline.Akima()),
    ("CubicSpline", Spline.Cubic()),
]
    rate_at(r, t, pt) = rate(zero(Yield.build_model(descriptor, t, r), pt))
    println("\n--- $name: ∂rate(t)/∂r₃  (bump at 5yr) ---")
    for pt in eval_points
        g = ForwardDiff.gradient(r -> rate_at(r, tenors, pt), rates)
        println("  t=$(lpad(pt,4)):  $(round(g[3], digits=4))")
    end
end
```

Results (sensitivity of interpolated rate to 5yr rate bump):

| t | **Linear** | **PCHIP** | **MonotoneConvex** | **Akima** | **CubicSpline** |
|---|-----------|-----------|-------------------|-----------|-----------------|
| 0.5 | 0.0 | -0.10 | 0.0 | -0.11 | 0.02 |
| 1.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 1.5 | 0.0 | -0.08 | -0.02 | -0.12 | -0.02 |
| 2.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 3.0 | 0.33 | 0.50 | 0.45 | 0.65 | 0.26 |
| **5.0** | **1.0** | **1.0** | **1.0** | **1.0** | **1.0** |
| 7.0 | 0.6 | 0.64 | 0.53 | 0.63 | 0.95 |
| 10.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 15.0 | 0.0 | **-0.23** | **-0.08** | **-0.42** | **-0.56** |
| 20.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

**Linear** is perfectly local within the knot grid — zero sensitivity outside adjacent intervals. **MonotoneConvex** has the best locality among smooth methods (only -0.08 at t=15 vs -0.23 for PCHIP, -0.42 for Akima, and -0.56 for CubicSpline). All smooth methods have zero sensitivity at the exact tenor points (t=1, 2, 10, 20) because the interpolation passes through those data points exactly.

Beyond the final knot, locality depends on the extrapolation policy as well. With
`:flat_forward`, for fixed tenors,

```math
\frac{\partial z(t)}{\partial r_i}
= \frac{t_n}{t}\,\mathbf{1}_{i=n}
+ \left(1-\frac{t_n}{t}\right)\frac{\partial f_n}{\partial r_i}.
```

The endpoint forward's sensitivities therefore persist throughout the tail. Cubic can
spread alternating positive and negative sensitivities across the whole grid; even linear
interpolation propagates sensitivity to its last two knots beyond the last interval.
Narrow final intervals can amplify those sensitivities. With `FlatForwardAt(f)` held
fixed, only the final knot rate contributes to the extrapolated zero rate, with weight
`tₙ/t`. These statements concern zero-rate sensitivities; cashflow PV sensitivities also
depend on maturity and discounting.

### Performance

All methods are fast enough for interactive use. The table below shows end-to-end `sensitivities()` timing from [ActuaryUtilities.jl](https://github.com/JuliaActuary/ActuaryUtilities.jl), which includes gradient + Hessian + result packaging in a single call:

```julia
using ActuaryUtilities, FinanceModels, Printf

rates5 = [0.02, 0.025, 0.03, 0.035, 0.04]
tenors5 = [1.0, 2.0, 5.0, 10.0, 20.0]
cfs5 = [5.0, 5.0, 5.0, 5.0, 105.0]

for (name, spline) in [
    ("PCHIP", Spline.PCHIP()),
    ("MonotoneConvex", Spline.MonotoneConvex()),
    ("Linear", Spline.Linear()),
    ("Akima", Spline.Akima()),
    ("Cubic", Spline.Cubic()),
]
    zrc = ZeroRateCurve(rates5, tenors5, spline)
    sensitivities(zrc, cfs5, tenors5)  # warmup

    N = 5_000
    t0 = time_ns()
    for _ in 1:N; sensitivities(zrc, cfs5, tenors5); end
    elapsed = (time_ns() - t0) / 1e3 / N
    @printf("  %-20s  %7.1f μs\n", name, elapsed)
end
```

`sensitivities()` (5 tenors):

| Method | Time |
|--------|------|
| **MonotoneConvex** | **5.9 μs** |
| **Linear** | 5.3 μs |
| **PCHIP** | 10.1 μs |
| **Cubic** | 10.1 μs |
| **Akima** | 15.2 μs |

`sensitivities()` (12 tenors):

| Method | Time |
|--------|------|
| **MonotoneConvex** | **40.3 μs** |
| **PCHIP** | 69.4 μs |
| **Akima** | 102.1 μs |
| **Cubic** | 112.8 μs |
| **Linear** | 131.2 μs |

MonotoneConvex is fastest at both sizes. At 12 tenors the advantage is substantial — roughly 2x faster than PCHIP and 3x faster than Linear.

## Recommendations

- **`Spline.MonotoneConvex()`** (default): Best for finance applications. Guarantees positive continuous forward rates, best KRD locality among smooth methods (-0.08 vs -0.23 for PCHIP), and fastest AD performance. Based on [Hagan & West (2006)](https://doi.org/10.1080/13504860600829233).
- **`Spline.PCHIP()`**: Good general-purpose alternative. Smooth forward curves, local sensitivity, monotonicity-preserving.
- **`Spline.Linear()`**: Use when you need localized sensitivities within the knot grid and don't need smooth forwards.
- **`Spline.Akima()`**: Alternative to PCHIP with different behavior near inflection points. Slightly more non-local leakage than PCHIP.
- **`Spline.Cubic()`**: Use when curve smoothness matters most and you accept non-local KRD effects (e.g. negative duration at distant tenors from a local rate bump).
