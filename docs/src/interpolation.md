# Interpolation Methods for `ZeroRateCurve`

`ZeroRateCurve` accepts an optional third argument specifying the interpolation method. The choice of interpolation affects **forward curve smoothness**, **knot-rate sensitivity locality**, and the cost of constructing and evaluating curves. Knot-rate sensitivity changes the interpolation inputs. ActuaryUtilities 6 `KeyRates` instead adds triangular continuous-zero shifts without refitting the curve; its shock shapes do not depend on the interpolant.

`ZeroRateCurve` returns a `Yield.MonotoneConvex` for `Spline.MonotoneConvex()` and a `Yield.Spline` for the other methods. Both are [`Yield.AbstractInterpolatedZeroCurve`](@ref FinanceModels.Yield.AbstractInterpolatedZeroCurve)s, as are the curves that spline [`fit`](@ref FinanceModels.fit)s return: read their knots with [`knot_rates`](@ref FinanceModels.Yield.knot_rates) and [`knot_tenors`](@ref FinanceModels.Yield.knot_tenors), and build a changed curve (other rates, tenors, method, or extrapolation) with [`reconstruct`](@ref FinanceModels.Yield.reconstruct).

## Available Methods

| Method | Zero-rate smoothness inside the knot grid | Knot-rate locality | Description |
|--------|-----------|----------|-------------|
| `Spline.MonotoneConvex()` | C1 | Local | **Default.** Hagan–West interpolation with continuous, nonnegative forwards when the discrete forwards are nonnegative. |
| `Spline.PCHIP()` | C1 | Local | Preserves monotonicity of the zero-rate inputs. |
| `Spline.Akima()` | C1 | Local | Uses nearby slopes to limit oscillation. |
| `Spline.Linear()` | C0 | Adjacent intervals | Simple; instantaneous forwards jump at knots. |
| `Spline.Quadratic()` | C1 | Global | Piecewise quadratic zero rates. |
| `Spline.Cubic()` | C2 | Global | Natural cubic spline; a knot change affects the whole grid. |
| `Spline.BSpline(d)` | Depends on degree and knots | Global | Degree-`d` B-spline; **not thread-safe** for concurrent evaluation on a shared curve (shared internal buffer). |

These smoothness statements concern tenor, not derivatives with respect to knot rates. C1 zero rates give continuous instantaneous forwards; they need not give differentiable forwards at the knots. Extrapolation and the short-end boundary can introduce additional forward jumps.

Since v6.0, `Spline.Linear()`, `Spline.Quadratic()`, and `Spline.Cubic()` return `PolynomialSpline` interpolants, which are fast and safe to evaluate concurrently across threads. Only `Spline.Linear()` is local; the quadratic and natural cubic splines are global. For a B-spline (e.g. as a least-squares fitting basis) request `Spline.BSpline(d)` explicitly; it is **not** thread-safe for concurrent evaluation on a shared curve.

`Fit.Bootstrap()` supports linear interpolation (`Spline.Linear()` or `Spline.BSpline(1)`): bootstrapping requires that each new knot leaves earlier segments unchanged. Fit other methods to the whole quote set at once with `fit(spline, quotes)`, which reprices every quote.

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

knot_rates(zrc_lin)                                 # read-only view of the rates, as Continuous rates
zrc_up = reconstruct(zrc_lin; rates = rates .+ 0.001)   # every knot 10bp higher
zrc_mc = reconstruct(zrc_lin; spline = Spline.MonotoneConvex())
```

## Before the First Knot

The DataInterpolations-backed methods hold the zero rate flat at the first knot's rate between
`t = 0` and the first knot. `Spline.MonotoneConvex()` instead interpolates from `t = 0`: the
Hagan-West construction treats the origin as a node, with the instantaneous forward there set by
its boundary condition. A bootstrapped curve's knots are exactly its quote maturities, and the
flat short end prices the first quote's interval.

## Extrapolation Beyond the Last Knot

Every knot-based curve defaults to `extrapolation=:flat_forward`: beyond the last knot
`(tₙ, zₙ)` the instantaneous forward is held at a constant anchor `fₙ`, and for `t > tₙ`

```math
z(t) = f_n + (z_n - f_n)\frac{t_n}{t}.
```

This preserves the discount factor at `tₙ` (discount factors are continuous there) and
converges to `fₙ` as the horizon increases. The zero rate is not flat unless `fₙ = zₙ`.
The anchor depends on the curve:

- **DataInterpolations-backed curves** (`Spline.Linear()`, `Quadratic()`, `Cubic()`,
  `PCHIP()`, `Akima()`, `BSpline(d)`, through `Yield.Spline`, `ZeroRateCurve`, or `fit`)
  anchor on the **last discrete forward**, the average continuously compounded forward over
  the last knot interval:

  ```math
  f_n = \frac{z_n t_n - z_{n-1} t_{n-1}}{t_n - t_{n-1}}.
  ```

  (For a single knot, `fₙ = z₁`.) This anchor depends only on the last two knots and does
  not depend on the interpolant, so a quadratic, cubic, or B-spline end piece cannot push
  the tail to an extreme or negative level. The instantaneous forward generally **jumps**
  at `tₙ`, from the interpolant's endpoint forward `zₙ + tₙz′(tₙ⁻)` to `fₙ`.
- **`Spline.MonotoneConvex()`** keeps its native boundary instantaneous forward `f(tₙ)`
  from the Hagan-West construction, so its forward curve is continuous at `tₙ`. That
  forward is built from the last discrete forwards and collared between 0 and twice the
  last discrete forward (when that forward is positive), so it is well behaved by
  construction.

For example, zero rates `[0.02, 0.025, 0.03, 0.035, 0.04]` at `[1, 2, 5, 10, 30]` have a
last discrete forward of 4.25%, which is the tail forward of every DataInterpolations-backed
curve. Anchoring on the interpolant's endpoint derivative instead would make the tail depend on the interpolation method, and can produce negative or extreme long-end rates.

These are general curve-extension choices, not implementations of an accounting basis
or prescribed regulatory extrapolation. For convergence to an independently chosen
ultimate forward rate, see [`Yield.SmithWilson`](@ref), which accepts `ufr` and a
convergence-speed parameter `α`; model selection and calibration remain explicit choices.

Pass the `extrapolation` keyword to `ZeroRateCurve`, `Yield.Spline`, `reconstruct`, or a spline
`fit` call to select the long-end behavior:

| Value | Long-end zero rate | Notes |
|-------|--------------------|-------|
| `:flat_forward` | `fₙ + (zₙ-fₙ)tₙ/t` | **Default.** `fₙ` is the last discrete forward (MonotoneConvex: its boundary forward `f(tₙ)`). |
| `Yield.FlatForwardAt(r)` | `f + (zₙ-f)tₙ/t` | Independent forward assumption `f`, the continuously compounded value of the rate `r`; usually introduces a forward jump. |
| `:flat_zero` | `zₙ` | Holds the zero rate constant; usually introduces a forward jump at `tₙ`. |
| `:linear` | `zₙ + z′(tₙ⁻)(t-tₙ)` | Extends the zero rate at its boundary slope. |
| `:extension` | Final interpolation piece | Restores the former DataInterpolations behavior and can become extreme far beyond the grid. Not available for `Spline.MonotoneConvex()`. |

```julia
curve = Yield.Spline(Spline.Cubic(), tenors, rates; extrapolation=:flat_zero)
zrc = ZeroRateCurve(rates, tenors, Spline.Cubic(); extrapolation=:linear)
quotes = ZCBYield.(Continuous.(rates), tenors)
fitted = fit(Spline.Linear(), quotes, Fit.Bootstrap(); extrapolation=:extension)
```

The selected value is available as `curve.extrapolation` on every knot curve. It is preserved by
fitting, `reconstruct`, and `Accessors.@set`; derived boundary quantities are recomputed when
knots change. Loss-fitting `Spline.MonotoneConvex()` returns a native `Yield.MonotoneConvex`,
with `Yield.instantaneous_forward`, for every supported policy.

`Yield.FlatForwardAt` takes a rate with an explicit compounding convention:

```julia
assumed = ZeroRateCurve(rates, tenors; extrapolation=Yield.FlatForwardAt(Continuous(0.035)))
annual = ZeroRateCurve(rates, tenors; extrapolation=Yield.FlatForwardAt(Periodic(0.035, 1)))
```

A bare number (`Yield.FlatForwardAt(0.035)`) is a `MethodError`: `Yield.Constant(0.035)` reads a
bare number as annual effective, so the convention must be stated. The rate is
stored as its continuously compounded value and stays fixed when knot rates are bumped or
fitted. Appending a synthetic far knot does not generally impose a chosen terminal
instantaneous forward and can also change interior interpolation.

## Fitting More Quotes Than Knots

`fit(spline, quotes)` places one knot at each quote maturity and reprices every quote. To fit a
curve to more quotes than you want knots (the bonds of a credit index, say, whose maturities
repeat and whose prices no single curve reprices), choose the knots and refit a curve that has
them. The refit is a least-squares fit, [`Fit.LeastSquares`](@ref FinanceModels.Fit.LeastSquares):
the knot rates that minimize the squared pricing residuals. (A refit whose knots are at the quote
maturities is the exact fit, `Fit.Exact()`.)

```julia
knots = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0]
c0 = ZeroRateCurve(fill(0.05, length(knots)), knots, Spline.Cubic())
bbb = fit(c0, bond_quotes)                 # e.g. 1000 `Quote(price, Bond.Fixed(...))`s
```

Unweighted residuals are price errors, which grow with a bond's duration, so long bonds dominate
the fit. Weights of `1 / (price * duration)^2` make each residual approximately a yield error:

```julia
w = [1 / (q.price * duration(q))^2 for q in bond_quotes]   # duration from ActuaryUtilities, say
bbb = fit(c0, bond_quotes, Fit.LeastSquares(weights = w))
```

The fit solves the least-squares problem by Levenberg–Marquardt, and returns only a minimum: at
the fitted knot rates, the Hessian of the sum of squares (with the residuals' curvature) is positive
definite, and the Newton step to the minimum is at most `1e-10`. Quotes with fixed cashflows
(`Cashflow` and `Bond.Fixed` instruments, which every FinanceModels quote type is) are priced
through their cashflow matrix: the curve is evaluated once per distinct cashflow date rather than
once per cashflow. The quotes must determine every knot rate (at least as many quotes with a
positive weight as knots, each knot moving some quote's price), or the fit throws. Fits of fixed
cashflows are differentiable with respect to the quote prices, coupons and weights; see
[Sensitivities Through Calibration](@ref).

Choose an interpolation that is smooth in its knot rates, such as `Spline.Linear()` or
`Spline.Cubic()`. Monotone convex, PCHIP and Akima curves switch formulas between regions of knot
rates, so the least-squares objective has kinks, and its minimum can lie on one, where the fit
cannot converge: it then throws a [`FitConvergenceError`](@ref FinanceModels.FitConvergenceError)
that says so. In a sweep of random 80-bond indices, about half the monotone convex and PCHIP fits
did, and none of the linear, quadratic or cubic ones.

## Key Tradeoffs

### Forward Curve Smoothness

The instantaneous forward is `f(t) = z(t) + t * z′(t)`. Linear zero-rate interpolation gives jumps at knots. PCHIP, MonotoneConvex, Akima, and Cubic give continuous interior forwards, although only the smoother choices also provide continuous forward derivatives. Check the required regularity when a model uses derivatives of the forward curve.

Use [`Yield.instantaneous_forward`](@ref FinanceModels.Yield.instantaneous_forward) to evaluate the forward. At a jump it returns the right-hand value; a symmetric finite difference at the knot would instead average the two sides.

```@example interpolation_forwards
using FinanceModels
using Markdown # hide

rates = [0.02, 0.025, 0.03, 0.035, 0.04]
tenors = [1.0, 2.0, 5.0, 10.0, 20.0]
methods = [
    ("Linear", Spline.Linear()),
    ("PCHIP", Spline.PCHIP()),
    ("MonotoneConvex", Spline.MonotoneConvex()),
    ("Akima", Spline.Akima()),
    ("Cubic", Spline.Cubic()),
]
points = [1.99, 2.0, 2.01, 4.99, 5.0, 5.01]
forwards = map(methods) do (name, method)
    curve = ZeroRateCurve(rates, tenors, method)
    [rate(Yield.instantaneous_forward(curve, t)) for t in points]
end
# The table reports percentages. # hide
header = "| Method | " * join(["f($t)" for t in points], " | ") * " |\n" # hide
rows = ["| " * name * " | " * join(string.(round.(100 .* fs; digits=2)) .* "%", " | ") * " |" for ((name, _), fs) in zip(methods, forwards)] # hide
Markdown.parse(header * "|---|" * repeat("---|", length(points)) * "\n" * join(rows, "\n")) # hide
```

MonotoneConvex preserves nonnegative continuous forwards when the input discrete forwards are nonnegative, under its default extrapolation. Other extrapolation policies can change the tail behavior. See [Hagan & West (2006)](https://doi.org/10.1080/13504860500396032).

### Knot-Rate Sensitivity Locality

Changing an interpolation knot can move rates beyond its adjacent intervals. The example below computes `∂z(t)/∂r₃`, the sensitivity to the five-year continuous zero-rate input. It rebuilds the curve under each perturbation; it is **not** ActuaryUtilities' `KeyRates` shock.

```@example interpolation_locality
using FinanceModels, ForwardDiff
using Markdown # hide

rates = [0.02, 0.03, 0.035, 0.04, 0.045]
tenors = [1.0, 2.0, 5.0, 10.0, 20.0]
methods = [
    ("Linear", Spline.Linear()),
    ("PCHIP", Spline.PCHIP()),
    ("MonotoneConvex", Spline.MonotoneConvex()),
    ("Akima", Spline.Akima()),
    ("Cubic", Spline.Cubic()),
]
points = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
responses = map(methods) do (name, method)
    [ForwardDiff.gradient(r -> rate(zero(ZeroRateCurve(r, tenors, method), t)), rates)[3]
     for t in points]
end
header = "| t | " * join(first.(methods), " | ") * " |\n" # hide
rows = ["| $t | " * join([string(round(g[i]; digits=3)) for g in responses], " | ") * " |" for (i, t) in enumerate(points)] # hide
Markdown.parse(header * "|---|" * repeat("---|", length(methods)) * "\n" * join(rows, "\n")) # hide
```

Linear interpolation is local to the adjacent intervals inside the grid. The other methods spread the effect according to their slope construction. The relative magnitudes in this example depend on the input curve; they are not universal rankings. At each unperturbed knot the sensitivity is zero, because interpolation still passes through that knot.

Beyond the final knot, locality depends on the extrapolation policy as well. With
`:flat_forward`, for fixed tenors,

```math
\frac{\partial z(t)}{\partial r_i}
= \frac{t_n}{t}\,\mathbf{1}_{i=n}
+ \left(1-\frac{t_n}{t}\right)\frac{\partial f_n}{\partial r_i}.
```

The anchor's sensitivities therefore persist throughout the tail. For the
DataInterpolations-backed curves the anchor is the last discrete forward, so only the last
two knots contribute, whatever the interpolant:
`∂fₙ/∂rₙ = tₙ/(tₙ - tₙ₋₁)` and `∂fₙ/∂rₙ₋₁ = -tₙ₋₁/(tₙ - tₙ₋₁)`. A narrow final interval
amplifies both. MonotoneConvex's boundary forward depends on the last few discrete
forwards. With `FlatForwardAt(f)` held fixed, only the final knot rate contributes to the
extrapolated zero rate, with weight `tₙ/t`. These statements concern zero-rate
sensitivities; cashflow PV sensitivities also depend on maturity and discounting.

### Performance

Benchmark the operation you actually need: repeated evaluation, rebuilding curves from perturbed knots, and differentiating through calibration have different costs. Precompile the operation before timing it, and keep the same inputs and derivative order when comparing methods.

In ActuaryUtilities 6, fixed-cashflow `sensitivities` uses analytic derivatives of its continuous-zero shocks. Callback valuations use ForwardDiff. `sensitivities(...)` defaults to first order; request `sensitivities(SecondOrder(), ...)` for convexity. Neither a `KeyRates` grid nor a callback automatically rebuilds the interpolation knots. To measure knot-rate risk, reconstruct the curve inside the valuation, as above.

## Choosing a Method

- **`Spline.MonotoneConvex()`** is the default for its treatment of forwards. Its nonnegative-forward guarantee depends on nonnegative input discrete forwards and the extrapolation policy.
- **`Spline.PCHIP()`** preserves monotonic zero-rate data with local slope construction.
- **`Spline.Linear()`** gives localized knot effects within the grid when forward jumps are acceptable.
- **`Spline.Akima()`** offers another local slope construction; inspect the resulting forwards on your data.
- **`Spline.Cubic()`** gives smoother interior forwards, with global knot effects and possible overshoot.

For every method, inspect the short end, the knots, and the extrapolation tail. A curve suitable for pricing the calibration instruments may still need a different interpolation or extrapolation choice for another use.
