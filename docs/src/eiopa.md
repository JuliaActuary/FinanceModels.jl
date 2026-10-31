# [EIOPA risk-free curves: background and valuation tutorial](@id eiopa-rfr)

`FinanceModels.EIOPA` implements the Solvency II risk-free rate methodology
applicable from **30 January 2027**. This tutorial explains curve construction,
starts from published basic rates, revalues an asset/liability portfolio to
calculate its volatility adjustment, and then introduces stresses and the
alternative market-observation input path.

All numerical inputs below are synthetic. They demonstrate the 2027 method and
are not a monthly EIOPA release. Legacy Smith–Wilson publications use a different
extrapolation method and cannot be imported as 2027 inputs.

## Why these curves exist

The best estimate of liabilities (BEL) discounts projected payments using a
term structure of risk-free rates. For fixed cashflows, its value is

```math
\mathrm{BEL}=\sum_t C_t D(t), \qquad D(t)=(1+z_t)^{-t},
```

where ``z_t`` is an annual-effective zero rate. A 60-year benefit depends on the
60-year discount factor, rather than on a coupon rate or short-term policy rate.

Long maturities often have too little reliable trading to determine the whole
curve directly. EIOPA therefore combines observable rates with an extrapolation
toward a long-term assumption. An undertaking may also apply an approved
volatility adjustment (VA), linked to asset spreads and the relative spread
sensitivity of its assets and liabilities. The upward/downward interest-rate
stresses provide scenarios for the interest-rate solvency capital requirement
(SCR) under the standard formula; the API supplies scenario discounts rather
than a complete SCR calculation.

| Term | Meaning and role |
|:-----|:-----------------|
| RFR | Risk-free rate term structure used for discounting. |
| DLT | Deep, liquid and transparent: selected market observations used to build the basic curve. |
| CRA | Credit risk adjustment deducted from raw market observations before constructing the basic curve. |
| FSP | First smoothing point: the maturity where extrapolation starts. |
| LLFR | Last liquid forward rate: the weighted forward assumption that starts the extrapolation. |
| UFR | Ultimate forward rate: the long-term forward assumption, rather than the zero rate at every long maturity. |
| ``\alpha`` | Positive convergence parameter controlling how quickly forwards approach UFR. |
| RCS | Risk-corrected spread used as an input to VA. Supplied separately from the basic curve. |
| CSSR | Credit spread sensitivity ratio, based on asset and liability revaluations and bounded to [0,1]. |
| PVBP | Present-value sensitivity to spread changes. Workbook diagnostics use per-unit decimal-rate amounts; multiply by 0.0001 for per-bp amounts. |
| VA | Undertaking-specific volatility adjustment applied through FSP and then extrapolated. |

### How the tail and VA work

Up to FSP, constant-forward interpolation makes cumulative log-discounts linear
between observations. Beyond FSP, the instantaneous continuously compounded
forward at horizon ``h=t-\mathrm{FSP}`` is

```math
f(t)=u+(\ell-u)e^{-\alpha h},
\qquad u=\log(1+\mathrm{UFR}),\quad \ell=\mathrm{LLFR}.
```

Integrating this forward preserves the discount factor at FSP. Forwards can jump
there, and then approach UFR. Zero rates approach UFR more gradually because they
also reflect accumulated discounting before FSP.

A decimal annual VA ``v`` adds ``\log(1+v)`` to forwards through FSP. At those
annual nodes the adjusted zero rate is ``(1+z_t)(1+v)-1``. It also changes LLFR by
``w_{\mathrm{FSP}}\log(1+v)``, where ``w_{\mathrm{FSP}}`` is the first LLFR weight.
UFR stays fixed, so the liability VA effect fades after FSP.

Assets have a different revaluation rule: Article 51a(2)(c) increases their
spreads at **all maturities**. Their notional spread shock does not fade at FSP.

## Choose an input path

Both types work with the same `EIOPA.calculate` and `EIOPA.publication` functions.

| Input type | Supply | FinanceModels does |
|:-----------|:-------|:-------------------|
| `EIOPA.PublishedRFR` | Published basic annual zero rates for years 1:FSP, official LLFR, FSP weight, UFR and ``\alpha`` | Import controls, reconstruct the tail, and calculate your CSSR/VA. |
| `EIOPA.RFRInputs` | Selected raw swap/government observations, CRA, FSP, all LLFR weights, UFR and ``\alpha`` | Deduct CRA, bootstrap/interpolate, calculate LLFR, and then calculate CSSR/VA. |

`publication` specifies **output rounding**; it accepts either input type.
The package does not download files or infer currency, valuation date or eligibility.
For an actual publication, select matching data and controls from
[EIOPA's RFR publications](https://www.eiopa.europa.eu/tools-and-data/risk-free-interest-rate-term-structures_en),
convert percentages to decimal rates, and retain the source date/version.

## A complete undertaking valuation from published basic rates

### 1. Import the basic prefix and official controls

```@example eiopa
using FinanceModels

# Synthetic annual-effective basic rates: 2.00%, 2.05%, ..., 2.95%.
basic_rates = [0.0200 + 0.0005 * (t - 1) for t in 1:20]
inputs = EIOPA.PublishedRFR(basic_rates;
    fsp=20, llfr=Continuous(0.032), w_fsp=0.40,
    ufr=Periodic(0.033, 1), α=0.11)
basic_result = EIOPA.calculate(inputs)
(discount_at_60=discount(basic_result.curves.basic, 60),
    bootstrapped=basic_result.bootstrapped)
```

Here `0.032` means a **continuous** LLFR of 3.2%, while `Periodic(0.033, 1)` means
an **annual-effective** UFR of 3.3%. The compounding distinction matters. The FSP
weight is a fraction: 40% is `0.40`. Each vector element corresponds to its year,
so select an actual file's prefix as `annual_basic_rates[1:fsp]`.

`bootstrapped === nothing` confirms that no bootstrap was performed. No CRA is
deducted again and no LLFR is estimated from rounded rates. The tail uses the
supplied FSP node and controls; its last published digit can differ from an
original basic tail constructed before the FSP node was rounded. Preserve the
official full basic series when exact reporting of that series is required.

### 2. Build the notional curves for revaluation

The notional VA is 85% of currency RCS, rounded to whole basis points on the
publication path. With RCS of 100bp it is 85bp. The example assumes an eligible
undertaking choosing an undertaking-specific adjustment (USA) of 80%.

```@example eiopa
rcs = 0.01  # 100bp; 1bp is 0.0001
usa = 0.80
notional = EIOPA.publication(inputs;
    va_inputs=EIOPA.VAInputs(rcs_currency=rcs, usa=usa))

# Turn published annual series into discount curves on the supplied grid.
as_curve(rates) = ZeroRateCurve(Periodic.(rates, 1),
    notional.maturities, Spline.ConstantForward())
basic = as_curve(notional.rates.basic)
liability_notional = as_curve(notional.rates.cssr)

# The asset's own curve already includes its existing 100bp spread.
asset_base = basic + Yield.Constant(Periodic(0.01, 1))
asset_notional = EIOPA.asset_spread_shock(asset_base,
    notional.diagnostics.notional_va)
notional.diagnostics.notional_va * 10000
```

For an actual release, build the baseline `basic` from the file's full published
basic series for the asset/liability revaluation. `notional.rates.basic` above
reconstructs the synthetic example's basic tail; the rounded FSP anchor can cause
last-digit differences from an original published basic series. CSSR and VA
construction still use the basic prefix and official controls in `inputs`.

`liability_notional` is the fading CSSR curve. `asset_notional` applies the parallel
shock to the asset's own curve, preserving its existing spread. For instruments
whose cashflows depend on spreads, apply the shock inside their valuation models
as well. These annual series cover years 1–150; do not use the interpolating
wrapper's default tail as a substitute for the EIOPA tail beyond that grid.

The notional VA uses unadjusted currency RCS even when USA is selected, matching
workbook `CSSR Calculation!H10`. USA affects final currency VA. TD §9.3.5 redefines
currency RCS using USA without explicitly resolving its use in the notional step;
this API follows the workbook convention.

### 3. Revalue the portfolio and calculate CSSR

Use a 30-year fixed-income asset paying 5 per year plus 100 principal, and
liabilities paying 6 per year for 60 years. Amounts use the same units.

```@example eiopa
assets = [Cashflow(5.0, t) for t in 1:30]
push!(assets, Cashflow(100.0, 30))
liabilities = [Cashflow(6.0, t) for t in 1:60]

valuations = (
    assets=present_value(asset_base, assets),
    assets_notional=present_value(asset_notional, assets),
    liabilities=present_value(basic, liabilities),
    liabilities_notional=present_value(liability_notional, liabilities))
va_inputs = EIOPA.VAInputs(rcs_currency=rcs, usa=usa, valuations=valuations)
adjusted = EIOPA.publication(inputs; va_inputs=va_inputs)
(cssr=adjusted.diagnostics.cssr,
    assets_pvbp=adjusted.diagnostics.assets_pvbp,
    liabilities_pvbp=adjusted.diagnostics.liabilities_pvbp,
    total_va_bp=10000 * adjusted.diagnostics.total_va)
```

For valuation inputs, PVBP is a finite-difference sensitivity: the value reduction
divided by the notional decimal VA. Positive PVBP means value falls when spreads
increase. CSSR clamps the ratio
of asset to liability PVBP to [0,1]. Article 51a sets CSSR to one when liability
PVBP is nonpositive; TD §9.3.2 omits that statutory exception.

Alternatively supply `sensitivities=(assets=..., liabilities=...)` in common units,
or a precomputed `cssr`. Only one CSSR source is accepted. If notional VA is zero,
supply sensitivities or CSSR rather than dividing value differences by zero.

### 4. Value liabilities with the undertaking's VA

```@example eiopa
liability_va = as_curve(adjusted.rates.va)
(basic_bel=valuations.liabilities,
    adjusted_bel=present_value(liability_va, liabilities),
    basic_discount_60=discount(basic, 60),
    adjusted_discount_60=discount(liability_va, 60))
```

A positive VA raises discount rates and lowers the value of these positive fixed
liability payments. Currency VA before rounding is `0.85 * CSSR * USA * RCS`.
The final VA differs from the 85bp notional revaluation shock because it
incorporates both CSSR and USA.

### 5. See why asset and liability shocks differ

The plot uses continuous curves with the calculated CSSR to show their shapes;
the publication path additionally rounds annual outputs. For the orange line,
the same basic starting curve isolates the effect of the parallel asset shock.
Plotting requires `CairoMakie` in addition to `FinanceModels`.

```@example eiopa
using CairoMakie
shape = EIOPA.calculate(inputs; va_inputs=EIOPA.VAInputs(
    rcs_currency=rcs, usa=usa, cssr=adjusted.diagnostics.cssr))
parallel = EIOPA.asset_spread_shock(shape.curves.basic,
    shape.diagnostics.notional_va)
curves = (shape.curves.basic, shape.curves.cssr, shape.curves.va, parallel)
labels = ("Basic", "Notional liability CSSR", "Undertaking liability VA", "Parallel asset shock")
colors = ("#2764A5", "#8856A7", "#16857B", "#D97724")
ts = collect(0.5:0.25:150)
fig = Figure(size=(1000, 460), fontsize=15)
axz = Axis(fig[1, 1], xlabel="Maturity (years)", ylabel="Annual zero rate (%)")
axf = Axis(fig[1, 2], xlabel="Maturity (years)", ylabel="Continuous forward (%)")
for (curve, label, color) in zip(curves, labels, colors)
    lines!(axz, ts, [100 * rate(Periodic(1)(zero(curve, t))) for t in ts]; label, color)
    lines!(axf, ts, [100 * rate(Yield.instantaneous_forward(curve, t)) for t in ts]; color)
end
for (ax, u) in ((axz, 3.3), (axf, 100log1p(0.033)))
    vlines!(ax, [inputs.fsp]; color=:gray50, linestyle=:dash)
    hlines!(ax, [u]; color=:gray50, linestyle=:dot)
end
Legend(fig[2, 1:2], axz; nbanks=2, framevisible=false, tellheight=true)
save("eiopa-curves.svg", fig)
save("eiopa-curves.png", fig) # hide
nothing # hide
```

![Annual zero rates and instantaneous forwards: liability VA fades toward the common UFR, while the parallel asset spread shock persists. Dashed vertical lines mark FSP; dotted horizontal lines mark UFR.](eiopa-curves.svg)

Liability forwards converge to the same UFR. The parallel asset spread shock
persists. Annual nodes can introduce forward jumps; at FSP the right-hand forward
equals the supplied LLFR, including its applicable VA change.

## Macro VA and eligibility

For macro VA set `use_macro_va=true`, with `rcs_country` and `debt_to_assets`.
For USA supply `usa` (at most 1.05). They are mutually exclusive. Euro-area
undertakings without USA must apply macro VA (TD §9.3.4), except the peer-country
cases prescribed to have zero macro VA. The default `false` does not infer
currency, domicile or eligibility; select the applicable treatment explicitly.

```@example eiopa
macro_inputs = EIOPA.VAInputs(rcs_currency=rcs, rcs_country=0.02,
    debt_to_assets=0.50, use_macro_va=true, cssr=adjusted.diagnostics.cssr)
macro_result = EIOPA.publication(inputs; va_inputs=macro_inputs)
(country_weight=macro_result.diagnostics.country_weight,
    macro_va=macro_result.diagnostics.macro_va,
    total_va=macro_result.diagnostics.total_va)
```

Spreads and VAs are decimal annual rates; negative values are supported.
Diagnostics expose basic/adjusted LLFRs, CSSR, asset/liability PVBP, country weight
(`country_weight`, ``\omega``), and notional/currency/macro/total VA. PVBP supplied
through sensitivities retains the caller's common units.

## Constructing the basic curve from market observations

Use `RFRInputs` when you have selected DLT market observations and want to
perform bootstrapping and LLFR calculation yourself. This independent synthetic
example starts with swap coupons:

```@example eiopa
market_inputs = EIOPA.RFRInputs([0.02, 0.025, 0.03, 0.032], [1, 5, 10, 20];
    frequency=1, cra_bp=10, fsp=10, llfr_weights=[0.0, 0.0, 0.4, 0.6],
    ufr=Periodic(0.033, 1), α=0.11)
market_result = EIOPA.calculate(market_inputs)
(llfr=market_result.diagnostics.llfr, discount_60=discount(market_result.curves.basic, 60))
```

Supply complete selected observations at increasing positive integer tenors and
aligned nonnegative LLFR weights summing to one (zero before FSP). Swap rates are
nominal annual coupons with `frequency` payments per year. Government observations
are annual-effective zero rates (`instrument=:government`). CRA is in basis points
and is deducted for either instrument. Currency adjustments must already be
reflected in observations. CRA/RCS estimation and previous-day fallback are
outside this API.

Use the assessed FSP explicitly. TD Table 2 includes assessed FSPs with zero
weight; these are accepted, and VA then contributes nothing to LLFR. The workbook
infers FSP from its first positive weight, a shortcut appropriate only when the
assessed FSP has positive weight. Supply ``\alpha`` rather than calibrating it:
0.40 for SEK, 0.11 otherwise, or the applicable approved phasing-in parameter.

`market_result.bootstrapped` retains all DLT maturities. The final basic curve
retains it through FSP; later observations influence LLFR only and need not be
repriced by the final extrapolated curve.

## Interest-rate stress scenarios

### From raw market observations

```@example eiopa
shocks = EIOPA.ShockParameters()
stressed = EIOPA.publication(market_inputs; shock_parameters=shocks)
(up_60=stressed.rates.up[60], down_60=stressed.rates.down[60])
```

Supply `va_inputs` from the matching portfolio revaluation to obtain `up_va` and `down_va`
too. Both scenarios retain the undertaking's total VA. Shocks apply to unrounded
bootstrapped annual rates, recompute LLFR, and move annual-effective UFR by ±15bp.
Downward floors apply to final rates, after VA and rounding. Adjusted downward
LLFR diagnostics refer to extrapolation before that final floor.

### From imported stressed bases

A published basic prefix and LLFR alone cannot determine stressed LLFRs: the
underlying DLT observations beyond FSP are missing. Supply separate up/down bases
with **unfloored** annual zero rates through FSP, their own LLFRs and already-shifted
UFRs. Obtain those from a matching source or your stress calculation. A downward
series whose floor already binds has lost the information needed to add VA first.

The following synthetic bases illustrate the input contract:

```@example eiopa
up_inputs = EIOPA.PublishedRFR([r * 1.20 + 0.01 for r in basic_rates];
    fsp=20, llfr=Continuous(0.043), w_fsp=0.40,
    ufr=Periodic(0.0345, 1), α=0.11)
down_inputs = EIOPA.PublishedRFR([r * 0.80 - 0.01 for r in basic_rates];
    fsp=20, llfr=Continuous(0.018), w_fsp=0.40,
    ufr=Periodic(0.0315, 1), α=0.11)
imported_stress = EIOPA.publication(inputs; va_inputs=va_inputs,
    shock_parameters=shocks, shocked_inputs=(up=up_inputs, down=down_inputs))
(up_va_60=imported_stress.rates.up_va[60], down_va_60=imported_stress.rates.down_va[60])
```

Here relative/absolute shocks are already in the inputs and are not applied again.
`shock_parameters` supplies final downward floors and checks UFR shifts. FSP,
its weight and ``\alpha`` must match the basic inputs. Missing stressed bases produce
an error. The synthetic transformations above are illustrative, rather than the
prescribed maturity-dependent shock table.

Article 46(1c)(a) and the workbook calculate the first shocked LLFR term over the
preceding DLT interval into FSP. This equals the one-year forward into FSP on the
basic constant-forward curve; after maturity-dependent shocks the intervals can
differ. The market-input path uses that DLT-interval convention and retains every
annual shocked node before FSP.

## Full precision, publication rounding and sensitivities

| Function | Output and purpose |
|:---------|:-------------------|
| `calculate` | Continuous discount-curve objects with unrounded VA, for analytical valuation and sensitivities. |
| `publication` | Annual-effective numeric rates with workbook-compatible staged rounding, for regulatory reconciliation. |

Regulatory technical provisions use published basic rates and whole-bp VA
(TD §§10.1.2–10.1.3). Rounding VA alone can change the adjustment by up to half a
basis point. `publication` rounds basic annual rates to five decimals before CSSR/VA,
notional and total VA to whole bp, and final annual rates to five decimals. It then
applies downward floors without another rounding pass, so a binding floor can
retain more digits. Decimal ties follow Excel's 15-significant-digit arithmetic
and round away from zero. Annual outputs default to years 1–150;
`maturities` accepts another increasing positive integer grid.

Use publication curves and rounded notional VA for revaluations supplied to
`publication`, as in the tutorial. For `calculate`, do the revaluation pass with
its precise basic/CSSR curves and notional VA. Publication cannot recover
precision already lost in imported data; the analytical path also starts from
exactly the rates and controls you supply.

First-order AD is supported away from clipping/floor boundaries; publication
rejects dual-number inputs. Maturities and FSP remain fixed. This example holds
the measured CSSR fixed and differentiates the analytical VA curve with respect
to currency RCS:

```@example eiopa
using ForwardDiff
value_at_rcs(s) = discount(EIOPA.calculate(inputs; va_inputs=EIOPA.VAInputs(
    rcs_currency=s, usa=usa, cssr=adjusted.diagnostics.cssr)).curves.va, 60)
ForwardDiff.derivative(value_at_rcs, rcs)
```

`calculate` returns `bootstrapped`, `curves` and `diagnostics`. Curve names are
`basic`, `cssr`, `va`, `up`, `down`, `up_va`, `down_va`; unavailable optional
outputs are `nothing`. `publication` returns `maturities`, corresponding `rates`
and diagnostics. Precise shocked curves interpolate pre-FSP annual nodes with
constant forwards. Custom `ShockParameters` interpolate parameters linearly and
hold endpoint values flat outside their grid.

## Reusable curve components

`Spline.ConstantForward()` works with `ZeroRateCurve` and
`fit(Spline.ConstantForward(), quotes, Fit.Bootstrap())`. It interpolates cumulative
log-discounts linearly, whereas `Spline.Linear()` interpolates zero rates. It
supports knot/reconstruction interfaces and usual tail policies except `:extension`.

`Yield.ExponentialForward(base; start, initial_forward, ufr, α)` preserves a curve
through `start`, then converges its instantaneous forward exponentially to UFR.
Both rate assumptions have explicit compounding. The start log-discount is cached;
`Accessors.@set` reconstruction validates assumptions and refreshes that cache.

## Sources and version

- [EIOPA-BoS-26/198, May 2026](https://www.eiopa.europa.eu/document/download/7da8a1f3-e395-4789-bdc1-76984f8510b2_en): §8 and Annex D for basic curves; §§9–10 for VA, including published-input controls in §10.1.
- [EIOPA workbook, 19 May 2026](https://www.eiopa.europa.eu/document/download/01e0930c-d78d-4806-8e83-a27609ea57c5_en): reference calculations and publication rounding.
- [EU 2026/269](https://eur-lex.europa.eu/eli/reg_del/2026/269/oj/eng): amendments to Articles 46, 51a, 165 and 166.

The default shock table is pinned to that regulation. UFR, FSP, selected observations,
LLFR controls, convergence parameters and risk-corrected spreads remain caller-supplied.
The regulation takes precedence over omissions in the technical document.
Offline tests reconcile the official workbook's seven annual curves and full-precision
intermediates; separate tests compare imported inputs with the market-input path.

## API reference

```@autodocs
Modules = [FinanceModels.EIOPA]
```
