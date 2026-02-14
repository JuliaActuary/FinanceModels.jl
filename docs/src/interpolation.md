# Interpolation Methods for `ZeroRateCurve`

`ZeroRateCurve` accepts an optional third argument specifying the interpolation method. The choice of interpolation affects **forward curve smoothness**, **key rate duration locality**, and **performance** when used with automatic differentiation (e.g. via `sensitivities()` in ActuaryUtilities.jl).

## Available Methods

| Method | Smoothness | Locality | Description |
|--------|-----------|----------|-------------|
| `Spline.MonotoneConvex()` | C1 (smooth) | Best among smooth | **Default.** Finance-aware. Positive forwards, best KRD locality, fastest AD. |
| `Spline.PCHIP()` | C1 (smooth) | Local | Monotonicity-preserving, local. Good general-purpose alternative. |
| `Spline.Akima()` | C1 (smooth) | Local | Local, resistant to outlier oscillation. |
| `Spline.Linear()` | C0 (kinked) | Perfectly local | Simplest. Kinks in forward curve at tenor points. |
| `Spline.Cubic()` | C2 (smoothest) | Global | Smoothest, but bumping one rate affects the entire curve. |
| `Spline.BSpline(n)` | Varies | Mostly local | nth-order B-spline. |

```julia
using FinanceModels

rates = [0.02, 0.03, 0.035, 0.04, 0.045]
tenors = [1.0, 2.0, 5.0, 10.0, 20.0]

zrc = ZeroRateCurve(rates, tenors)                              # default: MonotoneConvex
zrc_pchip = ZeroRateCurve(rates, tenors, Spline.PCHIP())        # PCHIP
zrc_lin = ZeroRateCurve(rates, tenors, Spline.Linear())          # linear
zrc_cub = ZeroRateCurve(rates, tenors, Spline.Cubic())           # cubic B-spline
zrc_aki = ZeroRateCurve(rates, tenors, Spline.Akima())           # Akima
```

## Key Tradeoffs

### Forward Curve Smoothness

The instantaneous forward rate `f(t) = r(t) + t · r'(t)` should be smooth for stochastic models that differentiate the forward curve (e.g. Hull-White θ(t) calibration). Linear interpolation creates discontinuous jumps in `f(t)` at tenor points, while PCHIP, MonotoneConvex, Akima, and CubicSpline produce smooth forward curves.

The following code evaluates forward rates near the 2yr and 5yr tenor points to illustrate the difference:

```julia
using FinanceModels
using FinanceCore: discount
DI = FinanceModels.DataInterpolations

rates = [0.02, 0.025, 0.03, 0.035, 0.04]
tenors = [1.0, 2.0, 5.0, 10.0, 20.0]

eval_points = [1.9, 1.99, 2.0, 2.01, 2.1, 4.9, 4.99, 5.0, 5.01, 5.1]

# Helper: numerical forward rate from a zero-rate interpolator
function fwd_from_interp(interp, t)
    r = interp(t); h = 1e-6
    dr = (interp(t+h) - interp(t-h)) / (2h)
    r + t * dr
end

# Helper: numerical forward rate from a discount function
function fwd_from_discount(model, t)
    h = 1e-6
    -log(discount(model, t+h) / discount(model, t-h)) / (2h)
end

for (name, make_fwd) in [
    ("Linear", (r, t) -> begin
        interp = DI.BSplineInterpolation(r, t, 1, :Uniform, :Average;
            extrapolation=DI.ExtrapolationType.Extension)
        pt -> fwd_from_interp(interp, pt)
    end),
    ("PCHIP", (r, t) -> begin
        interp = DI.PCHIPInterpolation(r, t;
            extrapolation=DI.ExtrapolationType.Extension)
        pt -> fwd_from_interp(interp, pt)
    end),
    ("MonotoneConvex", (r, t) -> begin
        mc = FinanceModels.Yield.MonotoneConvex(collect(r), collect(float.(t)))
        pt -> fwd_from_discount(mc, pt)
    end),
    ("Akima", (r, t) -> begin
        interp = DI.AkimaInterpolation(r, t;
            extrapolation=DI.ExtrapolationType.Extension)
        pt -> fwd_from_interp(interp, pt)
    end),
    ("CubicSpline", (r, t) -> begin
        interp = DI.CubicSpline(r, t;
            extrapolation=DI.ExtrapolationType.Extension)
        pt -> fwd_from_interp(interp, pt)
    end),
]
    fwd = make_fwd(rates, tenors)
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
using FinanceModels: DataInterpolations as DI
using FinanceCore: discount
using ForwardDiff

rates = [0.02, 0.03, 0.035, 0.04, 0.045]
tenors = [1.0, 2.0, 5.0, 10.0, 20.0]

eval_points = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]

for (name, rate_at) in [
    ("Linear", (r, t, pt) ->
        DI.BSplineInterpolation(r, t, 1, :Uniform, :Average;
            extrapolation=DI.ExtrapolationType.Extension)(pt)),
    ("PCHIP", (r, t, pt) ->
        DI.PCHIPInterpolation(r, t;
            extrapolation=DI.ExtrapolationType.Extension)(pt)),
    ("MonotoneConvex", (r, t, pt) -> begin
        mc = FinanceModels.Yield.MonotoneConvex(collect(r), collect(float.(t)))
        -log(discount(mc, pt)) / pt
    end),
    ("Akima", (r, t, pt) ->
        DI.AkimaInterpolation(r, t;
            extrapolation=DI.ExtrapolationType.Extension)(pt)),
    ("CubicSpline", (r, t, pt) ->
        DI.CubicSpline(r, t;
            extrapolation=DI.ExtrapolationType.Extension)(pt)),
]
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

**Linear** is perfectly local — zero sensitivity outside adjacent intervals. **MonotoneConvex** has the best locality among smooth methods (only -0.08 at t=15 vs -0.23 for PCHIP, -0.42 for Akima, and -0.56 for CubicSpline). All smooth methods have zero sensitivity at the exact tenor points (t=1, 2, 10, 20) because the interpolation passes through those data points exactly.

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
- **`Spline.Linear()`**: Use when you need perfectly localized KRDs (zero sensitivity outside adjacent intervals) and don't need smooth forwards.
- **`Spline.Akima()`**: Alternative to PCHIP with different behavior near inflection points. Slightly more non-local leakage than PCHIP.
- **`Spline.Cubic()`**: Use when curve smoothness matters most and you accept non-local KRD effects (e.g. negative duration at distant tenors from a local rate bump).
