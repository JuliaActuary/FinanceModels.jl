# FinanceModels.jl release notes

## v7.0.0 (unreleased)

### Valuation contexts (requires FinanceCore 3; changed numbers)

`present_value(ctx, contract)` is the single valuation entry point. Its first argument, the
valuation context, is a yield curve or rate, a model with a closed form for the contract, or
`Models(model, store)` / `Models(model; index)`, which also holds the models a contract reads by key.
`Models(curve)` is the single-curve case, `index = curve`.

- **Portfolios are correct.** FinanceCore passed each element's index of a collection as a third
  argument, which FinanceModels read as a valuation time: a 3-year zero-coupon bond at a
  continuous 3% was worth 0.9139 alone but 0.9418 inside `[bond]`, and `[a, b]` and `[b, a]`
  differed. A collection of contracts is now worth the sum of its contracts' values.
- **Closed forms compose.** A `Composite` is worth the sum of its parts, so options, caps, floors
  and swaptions (which have no cashflow projection) can be valued inside a `Composite` or a
  portfolio; they threw. A closed form is defined on the contract, `present_value(ctx,
  c::MyContract)`, reading `discount(ctx, t)`, `ctx[key]` and `valuation_model(ctx)`; it gives the
  same value under `Models(model, store)` as under `model`. Wrappers that act on the cashflow stream
  (`Forward`, `FX.Converted`, transducers) need a projection.
- **No valuation-time argument.** `present_value(model, contract, cur_time)` is removed. For a
  deterministic curve, the value as of `t` of the cashflows at or after `t` is a short reduction;
  see the migration guide. Cashflows before time 0 now accumulate; the default `cur_time = 0`
  dropped them.
- **Removed:** `Projection(contract; index)` (use `Models(model; index)`), `model_requirements`,
  and valuing a `Projection` (`present_value(model, ::Projection)`). Value the contract under a
  context, and use `collect(Projection(contract, ctx))` for its cashflows.
- **One reporting currency per context.** An `FX.Forwards` model discounts in its quote currency,
  and `FX.Forward` is one quote-currency cashflow. A base-currency `FX.BasisSwapLeg` is valued on a
  base-currency curve (`present_value(m.foreign, leg)`), converted with `FX.Converted`, or under a
  context that reports in its base currency. Under any other context (an FX model in another
  currency, `Models` over a plain curve, a `Rate`) it throws: it was valued in base-currency units,
  so a `Composite` would have added EUR to USD. An `FX.Converted` under an FX model's context (or
  `Models` over one) must convert into its reporting currency: a GBP→EUR conversion under a EURUSD
  model was discounted on the USD curve. Inside `FX.Converted`, a leg must pay in the pair's base
  currency and a nested `FX.Converted` must convert into it; otherwise it throws instead of being
  converted again. A plain bond carries no currency and is valued in the context's currency. An `FX.Forward`, which
  prices on the context's own FX model, is not supported there. Fitting an `FX.Forwards` model
  still values par basis-swap quotes on its foreign curve.

### Public extension hooks

The hooks for extending FinanceModels are public (Julia 1.11 `public`; FinanceModels now depends on
Compat), not exported. Two are renamed from 6.x:

- `FinanceModels.default_variables(model)`, was `__default_optic`: the parameters `fit` varies;
- `FinanceModels.default_optimizer(model)`, was `__default_optim`: `fit`'s default optimizer.

Two are new: `Yield.force_of_interest(curve, t)`, a curve's closed-form instantaneous forward as a
number, and `FinanceModels.closed_form(model, contract)`, the model kernel behind the values of
options, caps, floors and swaptions. A custom contract whose formula depends on the pricing model
defines `present_value(ctx, c) = FinanceModels.closed_form(valuation_model(ctx), c)`.
`FinanceModels.simulation_steps` is public too.

### Stable CIR bond prices (changed numbers)

The textbook Cox–Ingersoll–Ross price broke down at small volatility and long maturities: it raises
a number within σ² of 1 to the power 2ab/σ² and forms e^{γτ}. A 10-year discount factor was 5.6e57
at σ = 1e-10 and `Inf` at σ = 1e-12, and 10,000 years at σ = 0.01 gave `NaN`. One log-price kernel
now gives point, conditional and interval factors:

- for either sign of the mean reversion (a negative one at 10,000 years was `NaN` too);
- for volatilities whose square underflows, and at σ = 0, where it gives the deterministic price,
  so the price is continuous as σ → 0 (at τ = ∞ with b = 0 it was `NaN` at σ = 0);
- at τ = ∞ without mean reversion, where the rate is absorbed at 0 and P(∞) = exp(−2r/γ);
- in Float32, since overflow is detected rather than assumed from Float64's range;
- with derivatives through a = σ = 0 (they were `NaN`), and unchanged when the time unit and the
  parameters are rescaled together.

−log P agrees with an 8192-bit evaluation of the textbook formula to about 1e-13 relative wherever
it is finite in Float64; an explosive a < 0 whose price leaves Float64's range can still give
`NaN`. Intervals are differences of −log P, so they stay finite where both point factors underflow.

### Forward-starting contracts project on their own clock (changed numbers)

`Forward(s, c)` shifted `c`'s payment times by `s` but projected `c` against the unshifted
models, so a floating instrument fixed its coupons on the index from time 0: a one-year floater
starting at 2 paid 1.0202 instead of 1.0305 on an upward curve. The instrument now projects
against the models seen from `s` (a yield curve as `ForwardStarting(curve, s)`, an FX model with
its forward rate at `s` as spot), and its cashflows are then shifted by `s`. A model is rebased only
when the instrument reads it, so a plain discount rate as the context, a `NamedTuple` store and
unused store entries work as before; a flat `Rate` index is the same from any start. Fixed
instruments are unchanged.

Wrapper contracts also nest in any order: `Forward(Forward(…))`, `Forward(Composite(…))`,
`Forward(FX.Converted(…))` and `FX.Converted(Forward(…))` threw a `MethodError`. Projecting
them obeys the laws of time translation: a forward at 0 is the contract, nested starts add,
forwarding a composite forwards its parts, and conversion commutes with forwarding. Contracts
defined by either documented extension point, `asfoldable` or `__foldl__`, fold inside any wrapper.

### Transducers over contracts keep their order (changed numbers)

A chain of transducers over a contract applied in reverse inside a projection: `bond |> double |>
add1` projected as [2.1, 4.1] instead of [1.1, 3.1] (`collect(bond |> double |> add1)` was right),
so its present value was wrong too. The chain now applies in the order written, inside every
wrapper, and stateful and early-terminating transducers (`Take`, `Scan`) work over a contract inside
a projection, a portfolio or a `Composite`; they threw. A transducer that emits when it completes
(`Partition(n; flush = true)`, `PartitionBy`) works under an outer `Take`; it threw. A projection that
emits nothing (`Take(0)`, a `Filter` that drops everything, an empty portfolio) collects to an empty
vector, also with transducers applied to it (`collect(Projection(c) |> Take(0))`); it threw.

### Simulated paths are defined only on their simulated grid

A `RatePath` from `simulate` covers times from 0 to `horizon`, or to the first grid point beyond an
unaligned `horizon`. Evaluating it outside that range
(`discount`, interval discounts, `zero`, `forward`, `short_rate`, the present value of a later
cashflow, or `pv_mc` with an explicit `horizon` shorter than the contract) throws
DataInterpolations' `RightExtrapolationError` or `LeftExtrapolationError`. It used to extend the
last simulated step silently, so a one-year simulation priced a payment at year 10. A path's
instantaneous rate is right-continuous: at a grid time `short_rate` and
`Yield.instantaneous_forward` both give the slope of the step that starts there, as a `Continuous`
rate, and the last step's at the path's end. `RatePath` takes a `DataInterpolations.LinearInterpolation` only, the
interpolant `simulate` builds: its rate is a step's slope, which another interpolant would get
wrong. Its grid must start at t = 0 with the value 0, or construction throws an `ArgumentError`.
With UnicodePlots loaded, a path displays up to its last grid time (at most 30); it threw.

`simulate` takes its number of steps from `FinanceModels.simulation_steps(horizon, timestep)`, which
returns `(; nsteps, aligned)`. A horizon within `8eps` of a whole number of steps is aligned, so
roundoff adds no step: `horizon = 0.07, timestep = 0.01` took 8 steps and now takes 7, which changes
those paths. An aligned path ends at the horizon itself: `horizon = 0.3, timestep = 0.1` ended at
`0.30000000000000004`, so a payment just past the horizon was valued. An unaligned horizon takes
the first grid point beyond it, and `aligned` is `false`.

### Short-end limits and instantaneous forwards

`zero(curve, 0)` is the zero rate's limit, the short rate, for every curve with a discount factor
of 1 at time 0: it was `NaN` for Vasicek, CIR, Hull–White, Smith–Wilson, `ForwardStarting` and
custom curves. `Yield.instantaneous_forward(curve, t)` is defined for every curve and returns a
`Continuous` rate: in closed form
for knot curves, Nelson–Siegel(–Svensson), Vasicek, CIR, Hull–White, `Constant`,
`ForwardStarting`, and composite and scaled curves, and by differentiating the log-discount with
ForwardDiff for any other. A yield shift's forward at 0 is its zero rate there. A knot curve's
forward is the right-hand derivative of −log D: at a knot, the forward of the interpolant piece
that starts there, and from the last knot on the tail's. Hull–White takes its drift from its
curve's forward, so its drift and conditional bond prices at knot times match 6.x. Its conditional
bond price also uses its curve's log interval, so it no longer returns `NaN` once point discount
factors underflow, and Float32 models stay Float32. Nelson–Siegel(–Svensson) zero rates near t = 0
come from the loadings' Taylor series: at t = 1e-16 the closed form gave 2.5% for a curve whose
short rate is 2%.

### Explicit rate conventions: typed results and typed inputs

Rate-valued results are `Rate`s, so their convention travels with them:

- `knot_rates(curve)` and a knot curve's `rates` property return a read-only vector of
  `Continuous` rates, a view of the curve's numeric knots: indexing allocates nothing, and `copy`
  and `collect` give a `Vector` of rates. `rate.(knot_rates(curve))` gives the numbers, for
  example to seed a ForwardDiff gradient.
- `Yield.instantaneous_forward(curve, t)` and `short_rate(path, t)` return a `Continuous` rate. A
  curve with a closed-form forward defines `Yield.force_of_interest(curve, t)`, which returns
  the number; a custom curve still needs only `discount`.

FinanceCore reads a bare number as annual effective, so the continuous numbers these functions
returned gave wrong values without an error: on a curve with zero rates of 2%, 3%, 3.5% and 4% at 1,
2, 5 and 10 years, `discount(Yield.Constant(instantaneous_forward(c, 5)), 1)` was 0.96077 instead
of 0.95999, and `ZCBYield(knot_rates(c)[3], 5)` priced at 0.84197 instead of 0.83946. The same code
now gives the right values.

Inputs that took a number also take a `Rate`; numbers keep their meaning:

- **Knot rates** (`ZeroRateCurve`, `Yield.Spline`, `Yield.MonotoneConvex`, `reconstruct`,
  `@set curve.rates`) take numbers (continuously compounded), `Rate`s, or a mixture. A `Rate` is
  converted to its continuous value, so `reconstruct(curve; rates = knot_rates(curve) .+
  Continuous.(shift))` shifts in a stated convention.
- **Continuous inputs**: `Yield.SmithWilson`'s `ufr`, the short rate `r_t` of the Vasicek, CIR and
  Hull–White conditional prices, and the `r` and `q` of `eurocall` and `europut` take a number
  (continuously compounded) or a `Rate`, converted. A `Periodic` `ufr` gave a wrong `zero`, and
  `discount` threw.
- **Coupons, margins and interest-rate strikes** (`Bond.Fixed`, `Bond.Floating`,
  `FX.ParBasisSwap`, and the `Option.Cap`, `Option.Floor` and `Option.Swaption` strikes) take a
  `Periodic` rate of the contract's frequency and use its nominal rate: `Bond.Fixed(Periodic(0.05,
  2), Periodic(2), 10)` is `Bond.Fixed(0.05, Periodic(2), 10)`. A rate of another frequency throws
  an `ArgumentError`, because its nominal and its converted rate give different coupons (0.50 per
  100 face for a 5% annual rate on a ten-year semiannual bond at 4%); a `Continuous` rate is a
  `MethodError`. A typed coupon can differ from the number in the last bits, since a `Rate` stores
  its continuous equivalent.
- **Custom curves' zero rates**: composition and the long-run tail read a curve's `zero` in its own
  convention, and `TenorShift`/`ProjectedShift` rules receive it as a `Continuous` rate; a `zero` that
  returned `Periodic(0.06, 1)` was read as 6% continuous, and passed to a rule as periodic.

Model fields (volatilities, mean-reversion speeds, curve coefficients, the stored `ufr`), quote
coordinates (`implied_quote`) and sensitivities stay numbers. Vasicek's and CIR's `initial` short
rate is a `Rate`, as in 6.x: a number is stored as `Continuous`. The models guide has a table of the
conventions. No value changes: the typed results hold the same numbers, bit for bit.

### Cashflow matrices keep the amounts' type; Smith–Wilson coupon sensitivities

`FinanceModels.cashflows_timepoints` builds its matrix in the amounts' floating-point type:
`BigFloat` amounts are no longer narrowed to `Float64`, and dual-number amounts no longer throw. So
`fit(Yield.SmithWilson(…), quotes)` can be differentiated with ForwardDiff with respect to
cashflow amounts such as `Bond.Fixed` coupons, not only quote prices. The matrix is built in one
pass over the cashflows (about twice as fast for 30 semiannual swap quotes), and payments at times
`-0.0` and `0.0` share a row (each was counted in both).

### `maturity` for forwards and European options

`maturity` is defined for `Forward` (its start plus the instrument's maturity), `Option.EuroCall`
and `Option.EuroPut`. For a forward-starting contract this makes `pv_mc`'s default horizon work;
`pv_mc` still cannot value a European option, which has no scenario cashflow projection.

### Empty values keep the model's number type

Present values of contracts with no cashflows, and of caps and floors with no caplets, are a zero
in the model's number type (`BigFloat`, ForwardDiff dual numbers) instead of a `Float64` `0.0`
(#288). A contract's present value starts its sum from that typed zero, so a Float32 model values
Float32 cashflows in Float32 (they were widened to Float64). A `Bond.Fixed` still projects Float64
amounts, so its value is Float64. Float64 results are unchanged.

### Optimization 5

FinanceModels requires Optimization 5 (was 4.4), OptimizationOptimJL 0.4.6 and AccessibleModels
0.1.14. `fit`'s `optimizer` keyword still takes any OptimizationOptimJL solver. With
OptimizationOptimJL 0.4.9 or later (Optim 2), a loss fit can stop at a different point within the
solver's default tolerance (gradient norm `1e-8`), so fitted knot rates can differ from Optim 1
results by about `1e-9`. Pass `solve_kwargs = (; g_tol = 1e-12)` for a fit that reprices its quotes
to rounding.

### `implied_quote`

`implied_quote(curve, family, maturity)` returns the quote at which a quote
constructor (`CMTYield`, `OISYield`, `ZCBYield`, `ZCBPrice`, or a closure such as
`(r, t) -> ParYield(r, t; frequency = 1)`) reprices on a curve. Its first-order
ForwardDiff derivatives with respect to curve parameters are exact: they come
from the implicit function theorem at the solution rather than from solver
iterations. FinanceModels now depends on ForwardDiff directly.

### Correct swaption Greeks under automatic differentiation

Jamshidian's decomposition prices a swaption as bond options struck at `P(T₀,Tᵢ;r*)`,
where the critical rate `r*` depends on the model parameters. `r*` was found by a
bisection that returned a plain `Float64`, so ForwardDiff Greeks of Vasicek and
Hull–White swaptions silently dropped every `∂Kᵢ/∂r*·dr*/dθ` term: on a 1×5 payer
swaption, Vasicek's `∂/∂a` had the wrong sign and vega was 9% too high. `r*` now carries
its implicit-function derivatives; prices are unchanged up to root-finder precision (#290).

### Vasicek and Hull–White near zero mean reversion (changed numbers)

Vasicek bond prices used a truncated Taylor expansion below |aτ| = 0.02 and, above it, a closed
form whose variance terms cancel: `-log P` was off by up to 4.5e-7 (a = 0.001, τ = 19.9), a
Vasicek swaption by up to 8e-8, and the bond price's derivative in `a` by 6e-5 relative at
a = 0.001. Hull–White's closed forms switched to their a = 0 limits below |a| = 1e-12, so their
derivatives in `a` vanished there (a swaption's `∂/∂a` at a = 0 was 0) and they lost digits just
above the switch. Both models now evaluate these factors to within a few units of Float64
rounding for every mean reversion, zero and negative included, with exact ForwardDiff
derivatives. Prices change by the former errors; away from small aτ, only in the last bits.
Vasicek's discount factor at τ = ∞ with mean reversion a > 0 and no drift term (b = σ²/(2a²)) is its
limit, exp(−r/a − σ²/(4a³)); it was `NaN`.
With explosive mean reversion (aτ < −1), Vasicek's bond price and forward group the terms that grow
like e^{|a|τ}: they cancelled when the short rate was near its long-run level, so a constant 3% rate
(σ = 0, r = b, a = −0.1) priced its 500-year bond at 1 instead of e^{−15}, and its forward was wrong
the same way. Prices and forwards for aτ ≥ −1 are bitwise unchanged.

### Differentiable spline fits

Spline `fit`s (`Fit.Loss` with any interpolation method, `Fit.Bootstrap()` with
`Spline.Linear()`, and `FX.Forwards` with a spline foreign curve) now propagate ForwardDiff
dual numbers in quote prices, rates, and cashflow amounts, and in a `Yield.FlatForwardAt`
extrapolation forward. The fitted knot rates carry the exact first-order derivatives of the
calibration from the implicit function theorem; their values are the primal fit, bitwise.
Dual maturities, nested dual numbers, dual numbers from two differentiations, loss fits that
do not reprice their quotes, and other models' fits throw an `ArgumentError` (#290).

An optimizer-backed fit now throws an `ArgumentError` when a caller's dual number reaches its
parameters or loss through data it does not differentiate: a contract type `fit` does not strip
(an option strike, a `Forward`-wrapped amount), or any input of a non-spline fit. It throws at the
first evaluation where the dual number appears. Before, some reached the optimizer, and ForwardDiff
returned the derivative of the solver's iterations (`0.0` under `NelderMead`). Each running fit has
its own ForwardDiff tags, so the check also works for a fit nested in another fit's loss, even when
its quotes sit in an untyped vector (`Any[]`, `Quote[]`) or depend on the enclosing fit only away
from the starting point. Bootstrap steps and implicit roots check their residuals the same way.

### Derivatives at interpolation kinks

`Spline.MonotoneConvex()`, `Spline.PCHIP()`, and `Spline.Akima()` are only piecewise smooth in
their knot rates. Where they switch formula (a flat stretch of the curve, equal adjacent
forwards, a straight run of Akima knots), ForwardDiff previously returned a derivative that
matched neither an up nor a down bump: on a flat MonotoneConvex curve the knot sensitivities
could be off by a factor of two, flat PCHIP gave `NaN`, and a flat MonotoneConvex curve on a
dense grid gave derivatives driven by rounding noise.

- MonotoneConvex now reports, for each partial, the limit of a centered bump in its direction.
  On a flat stretch of the curve these partials need not sum to the parallel-shift derivative,
  so they do not aggregate like a gradient. Primal values are unchanged.
- PCHIP and Akima throw an `ArgumentError` when dual knot rates move a kink, including the point
  where Akima switches to its fallback slope and its value jumps.
- A differentiated `fit` throws when the fitted curve lies on a kink or within the fit's
  precision of one.

Away from kinks, derivatives are unchanged. See "Sensitivities Through Calibration".

### Notional-independent solver checks

`implied_quote`, the swaption critical rate, and differentiated fits judged a vanishing slope or
an ill-conditioned calibration on an absolute scale, so a quote family with a tiny notional was
refused. A bootstrap likewise solved each quote to an absolute tolerance, so a quote with a
notional of `1e-10` was fitted only to about `1e-5` in its zero rate, and `implied_quote`
accepted a root on an absolute residual (a zero-coupon quote with a notional of `1e-10` was off
by about `1e-5`; at `1e-14` it returned the starting guess). Every check and solve is now
relative to the size of the quote, and `implied_quote` and the swaption critical rate verify the
root they accept.

### Solver settings for `fit`

Optimizer-backed `fit`s accept `solve_kwargs`, passed to `Optimization.solve` (for example
`solve_kwargs = (; maxiters = 10_000, g_tol = 1e-12)`). A differentiated loss fit that stops
too far from an exact fit now has a way to tighten it.

### Refitting a knot curve

`fit(curve, quotes)` for a knot curve (`Yield.AbstractInterpolatedZeroCurve`) is now the spline
fit on the curve's knots: it keeps the curve's tenors, interpolation method and extrapolation
policy, starts from the same knot rates as `fit(spline, quotes)` (not the curve's own), uses the
spline's default optimizer, and carries implicit-function derivatives of dual quotes when there
is one knot per quote. It previously ran the generic optic fit, whose candidates went through the
public constructor, so a flat PCHIP or Akima starting curve threw. It no longer accepts
`variables`, and `FinanceModels.KnotRatesOptic` is removed. An `FX.Forwards` with a knot-curve
foreign curve refits it through the implied foreign quotes.

### PCHIP and Akima fits on evenly spaced maturities

`Spline.PCHIP()` and `Spline.Akima()` loss fits started from knot rates that were linear in the
knot number. On evenly spaced maturities those are collinear, where Akima switches formula and
its derivatives are `NaN`, so a coupon-bearing Akima fit on maturities `1:8` failed to converge.
Both now start from a curve that is strictly increasing and concave in the maturities, away
from every formula switch.

### `discount(curve, Inf)` under a zero tail forward

With a tail forward of zero (for example `Yield.FlatForwardAt(Continuous(0.0))`),
`discount(curve, Inf)` returned `NaN` from `exp(-0 * Inf)`. It now returns the limit: the
discount factor at the last knot. It returns 0 for a positive tail forward and `Inf` for a
negative one. Differentiating it where the long-run forward is zero throws, since any bump moves
the limit to 0 or `Inf`. Under `extrapolation = :extension`, `discount(curve, Inf)` throws a
`DomainError`: that policy's polynomial continuation is evaluated at finite times only, and
returned `NaN` at `Inf` even for a flat extension.

### Bootstrap requires linear interpolation

`Fit.Bootstrap()` now accepts only `Spline.Linear()` and throws an `ArgumentError`
for other strategies. Smoother interpolants let later knots reshape earlier segments,
so bootstrapped coupon quotes silently stopped repricing (residuals up to 0.40% of
par). Fit smoother curves to all quotes at once with `Fit.Loss`. Bootstrap also
validates its inputs and reprices every quote on the returned curve. See the
migration guide.

### `ZeroRateCurve` returns the curve it builds; one knot-curve interface

`ZeroRateCurve` is now a function returning a `Yield.MonotoneConvex` (for the default
`Spline.MonotoneConvex()`) or a `Yield.Spline`. Both subtype the new
`Yield.AbstractInterpolatedZeroCurve`, as do spline `fit` and bootstrap results. Read knots with
`knot_rates` (`Continuous` rates) and `knot_tenors` and build a changed curve with `reconstruct(curve; rates, tenors,
spline, extrapolation)`, which also accepts dual numbers for knot-rate gradients. Knot curves copy
and validate their inputs through one shared construction (unsorted, duplicate, negative or
non-finite tenors, non-finite rates, and too few knots throw an `ArgumentError`), keep read-only
knot vectors, compare structurally (`==`/`isequal`/`hash` on knots, method and extrapolation
policy), throw a `DomainError` for negative times, and print as the `ZeroRateCurve(...)` call that
rebuilds them. Fitting a knot curve rebuilds it once per optimizer candidate instead of once per
knot. `Yield.build_model` and the `Yield.MonotoneConvex()` fit placeholder are removed:
`Spline.MonotoneConvex()` is the one monotone convex selector (#272), fitted through the same path
as every other interpolation method. `Spline.PolynomialSpline` accepts only orders 1 to 3 and
`Spline.BSpline` requires degree 1 or more.

Each construction form has one signature: `ZeroRateCurve(rates, tenors, spline =
Spline.MonotoneConvex(); extrapolation)` takes the method positionally, as in 6.x, and the sampling
form keeps its `ZeroRateCurve(curve, tenors; spline, extrapolation)` keyword. `Yield.Spline`
accepts only the DataInterpolations methods, so `Yield.Spline(Spline.MonotoneConvex(), …)` is a
`MethodError`, as in 6.x; `ZeroRateCurve` and `Yield.MonotoneConvex` build that curve.

### Flat zero rate before the first knot (changed numbers)

DataInterpolations-backed curves now hold the first knot's zero rate between `t = 0` and the
first knot instead of extending the first interpolation piece. **Zero rates before the first knot
change**; values from the first knot on do not, and `Spline.MonotoneConvex()` does not change.
Bootstrapped curves are numerically unchanged, but their knots are now exactly the quote
maturities (the extra `t = 0` knot is gone), so a curve has one knot rate per quote.

### Failed optimizer fits throw `FitConvergenceError`

Every optimizer-backed `fit` now checks the solver's return code and throws a
`FitConvergenceError` (with the solver's `retcode`) instead of returning the
unfitted starting model.

### Flat-forward extrapolation beyond the last knot (changed numbers)

DataInterpolations-backed curves (`Spline.Linear()`, `Quadratic()`, `Cubic()`,
`PCHIP()`, `Akima()`, `BSpline(d)`) no longer continue their final interpolation
piece beyond the last knot. By default they now hold the last discrete forward,
`(zₙtₙ - zₙ₋₁tₙ₋₁)/(tₙ - tₙ₋₁)`, constant there, which keeps discount factors
continuous at the last knot and stops the end polynomial from driving the tail
(continuing it, or even holding its end slope's forward, can give extreme or
negative long-end forwards).
**Zero rates and present values beyond the last knot change** for every such
curve and fit; values up to the last knot do not. For a linear bootstrap of
`ZCBYield.([0.02, 0.025, 0.031, 0.036], [1, 2, 5, 10])` the 30-year continuous
zero rate moves from 5.47% to 3.86% (−161 bp; the 30-year discount factor rises
from 0.194 to 0.314). `Spline.MonotoneConvex()`, the `ZeroRateCurve` default,
keeps its boundary forward and is unchanged.

Pass `extrapolation = :extension` to `Yield.Spline`, `ZeroRateCurve` or `fit` to
keep the previous values. The other policies are
`:flat_zero`, `:linear` and `Yield.FlatForwardAt(rate)`, which takes a
`FinanceCore.Rate` such as `Continuous(0.035)` (a bare number throws). See the
migration guide.

### Quote conventions

- `OISYield` builds annual-pay par swaps beyond one year (was quarterly), matching
  SOFR, €STR, and SONIA overnight index swaps.
- `ParSwapYield` and `InterestRateSwap` require an explicit `frequency`.
- `ParYield` throws an `ArgumentError` when an explicit `frequency` conflicts with a
  `Periodic` rate's own compounding; previously the keyword was ignored.
- Every `frequency` input (`par`, `ParYield`, `ParSwapYield`, `InterestRateSwap`, `FX.ParBasisSwap`,
  `Bond.Fixed`, `Bond.Floating`, `Option.Cap`, `Option.Floor`, `Option.Swaption`) takes an integer
  or a `Periodic`. `par` threw for a `Periodic`, `FX.ParBasisSwap` for an integer, and the bonds took
  no integer. Contracts store a `Periodic`. A non-integer number, which `Cap`, `Floor` and
  `Swaption` accepted, and a `Continuous` bond frequency are a `MethodError`.

See the migration guide for upgrade steps.

### `TransformedYield` removed

The `Yield.TransformedYield` alias, deprecated in 6.1, is removed. Use `Yield.TenorShift`.

### Constructor errors

`Yield.NelsonSiegel`, `Yield.NelsonSiegelSvensson`, `Yield.CairnsPritchard` and
`Yield.CairnsPritchardExtended` throw a `DomainError` whose value is the offending parameters; the
message was passed as the value. `Yield.SmithWilson` with `u` and `qb` of different lengths throws an
`ArgumentError` (was a `DomainError`).

### `CompositeYield` accepts only `+` and `-`

`curve1 + curve2` and `curve1 - curve2` multiply and divide the curves' discount factors, so
every interval factor of the result is the product or quotient of the components' and the
composition commutes with `ForwardStarting`. `Yield.CompositeYield(a, b, op)` with any other
`op` (for example `max`, or `*` of two zero rates) now throws a `MethodError`: such an operation
builds a new curve from zero rates measured from time 0 rather than composing the two curves.
For a pointwise transformation of one curve's zero rates, use a `TenorShift`
(`curve + ((z, t) -> ...)`), such as `curve + ((z, t) -> max(z, Continuous(0.0)))` to floor the
zero rate. Ordering rates with `max` needs FinanceCore 2.6, which FinanceModels now requires.

### Interval factors that don't underflow

For the built-in curves, `discount(curve, from, to)`, `accumulation(curve, from, to)` and
`forward(curve, from, to)` are now computed from the cumulative log-discount L(t) = −log D(t) at
the two endpoints, rather than as a ratio of discount factors. The main effects:

- **Far-tail intervals are finite.** For example `discount(Yield.NelsonSiegel(1.0, 0.05, -0.02, 0.01), 20000, 20001)`
  was `NaN` (0/0) and is now `exp(-0.05)`. `ForwardStarting` far into the tail is finite too.
- **Intervals from time 0 are unchanged.** `discount(curve, 0, t)` equals `discount(curve, t)` bit for
  bit, so this change leaves every contract's `present_value` unchanged. Other intervals can move by a few units
  in the last place, and `zero`/`forward` of `SmithWilson`, `ForwardStarting` and the short-rate
  models by a few more (they no longer round-trip through the discount factor).
- **Empty intervals.** `discount(curve, t, t)` is exactly 1, also at `t = Inf`.
- **Composite curves at infinity.** `discount(a + b, Inf)` (and `-`, scaling and `ForwardStarting`)
  combines the components' long-run behavior before taking the limit, so long-run forwards that
  cancel leave their finite intercept: a flat-forward `ZeroRateCurve` minus a `Constant` at its
  tail forward tends to a finite factor rather than `NaN`. Components with a closed-form tail are
  `ZeroRateCurve`/`Yield.Spline`, `MonotoneConvex` and `Constant`; any other curve contributes its
  zero rate at infinity, as the sum of zero rates did before. A zero rate that tends to ±Inf says
  the log-discount grows faster than linearly but not how fast, so against a `:linear` tail of the
  other sign the limit is unsupported and returns `NaN` rather than a guess.
- **Smith–Wilson** intervals are exp(−ufr·(to − from))·(1 + s(to))/(1 + s(from)). That is exact for either sign
  of the discount factor (a fit to arbitrary prices can make it negative) and finite in the far
  tail; `discount(sw, t)` is unchanged.
- **Custom curves** still need only `discount(curve, t)`, and their intervals stay the ratio
  D(to)/D(from). `forward(curve, 0, t)` is now consistent with that for a curve with D(0) ≠ 1
  (previously it assumed D(0) = 1).
- **Composite and scaled curves** take their components' intervals, combined: each component keeps
  its own form (a custom curve's ratio, Smith–Wilson's far-tail-stable ratio), and an infinite
  endpoint comes from the combined tail.
- **`forward(curve, from, to)`** is the interval's log-discount per unit time, taken by the same
  rule as `discount(curve, from, to)`. On a Smith–Wilson curve, or a custom curve with negative
  discount factors, it therefore exists where both factors are negative (it was a `DomainError`),
  and on composite, scaled and
  `ForwardStarting` curves it can move by a few units in the last place.

### Time derivatives at `t = 0`

`ForwardDiff` derivatives with respect to time at exactly `t = 0`, such as
`ForwardDiff.derivative(t -> discount(curve, t), 0.0)` or a derivative in the start of
`discount(curve, from, to)`, were `NaN` for `Yield.MonotoneConvex` (the `ZeroRateCurve` default),
`NelsonSiegel` and `NelsonSiegelSvensson`, and for composites, scaled curves and yield shifts over
them. Their zero rates had a removable 0/0 at 0, which a dual time skips past. They are now exact:

- `MonotoneConvex` computes its cumulative log-discount directly instead of multiplying its zero
  rate back by `t`, which moves some discount factors in the last bit. Over its first interval its
  zero rate uses a divided form of the Hagan–West integral, which is also more accurate for very
  short maturities (the previous form lost digits as `t → 0`).
- `NelsonSiegel` and `NelsonSiegelSvensson` take their loadings near `t = 0` from Taylor series
  (see "Short-end limits").

A curve without its own `zero` (Smith–Wilson, Vasicek, CIR, `ForwardStarting`, custom curves) has
only L(t)/t. Its zero rate at 0 is the limit, the short rate, but has no time derivative there: the
zero rate, and a yield shift over it, throw a `DomainError` for a time derivative at 0 instead of
returning `NaN`. Its discount factor has one.

## v6.4.0

### `Spline.BSpline` fitted values changed on non-uniform tenor grids

The only change in this release is the move from DataInterpolations 8 to 9 (#293).
That bump silently changed what `Spline.BSpline(d)` produces whenever the knot
tenors are not equally spaced, which includes every standard CMT/Treasury grid.
Fits on equally spaced tenors are bit-identical to v6.3.0.

**Mechanism.** DataInterpolations 8 evaluated `BSplineInterpolation` as
`S(p(t))`: the maturity `t` was first mapped piecewise-linearly onto an equally
spaced parameter grid (one parameter value per knot), and the degree-`d` spline
was built in that parameter space. FinanceModels passed the `:Uniform`
parameter-vector option. In effect, `BSpline(d)` was a degree-`d` spline in
*knot index*, warped back to maturity. Where two neighbouring tenor intervals
differ in width the warp has a slope break, so the zero curve had a derivative
kink and the instantaneous forward rate jumped at that knot. On the CMT grid
below the v6.3.0 `BSpline(4)` forward curve jumps by roughly −63 bps at the 10y
knot; `BSpline(3)` by −71 bps.

DataInterpolations 9 removed the parameter vector. The spline is now built and
evaluated directly in maturity with de Boor averaged knots, so `BSpline(d)` is
what its name says: a C^(d−1) piecewise polynomial in `t`. Forwards are
continuous at every knot. The old curve cannot be reproduced; the concept no
longer exists upstream.

**Consequence.** A true degree-4 interpolant across tenors spanning 1 month to
30 years oscillates between the long knots. With the representative CMT quotes

```julia
tenors = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]
quotes = [0.0530,0.0528,0.0525,0.0515,0.0490,0.0450,0.0430,0.0420,0.0425,0.0435,0.0465,0.0455]
c = fit(Spline.BSpline(4), CMTYield.(quotes, tenors), Fit.Bootstrap())
```

the fitted curve reprices its own 20y quote with an error of −20.8 bps
(v6.3.0: +1.0 bps) and the instantaneous forward at the 20y knot is negative.
`BSpline(3)` on the same grid improves, from 1.0 bps to 0.2 bps at 20y.
`BSpline(2)`, `Linear`, `Cubic`, `PCHIP`, `Akima` and `MonotoneConvex` are
unaffected or slightly better.

**What to do.**

- If you pinned curve values from a `BSpline` fit on unequally spaced tenors,
  re-baseline; the change is expected.
- Prefer `Spline.BSpline(3)` (or `Cubic`, `PCHIP`, `MonotoneConvex`) for
  Treasury-style grids. Degree 4 and above on strongly non-uniform tenors will
  oscillate; that is inherent to a global interpolating spline and not a bug.
