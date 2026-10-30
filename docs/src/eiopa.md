# EIOPA risk-free term structures

`FinanceModels.EIOPA` implements the new methodology applicable from 30 January
2027: constant-forward bootstrapping, LLFR/UFR extrapolation, undertaking-specific
volatility adjustments (VA), and standard-formula interest-rate shocks.

## Sources and version

- [EIOPA-BoS-26/198, May 2026](https://www.eiopa.europa.eu/document/download/7da8a1f3-e395-4789-bdc1-76984f8510b2_en): §8 and Annex D for basic curves; §§9–10 for VA.
- [EIOPA workbook, 19 May 2026](https://www.eiopa.europa.eu/document/download/01e0930c-d78d-4806-8e83-a27609ea57c5_en): reference calculations and publication rounding.
- [EU 2026/269](https://eur-lex.europa.eu/eli/reg_del/2026/269/oj/eng): amendments to Articles 46, 51a, 165 and 166.

The default shock table is versioned to that regulation. UFR, FSP, DLT observations,
LLFR weights, convergence parameters and risk-corrected spreads are supplied by
the caller. No data downloading, CRA/RCS estimation or previous-day fallback is
performed. Currency adjustments must already be reflected in observations; CRA
is deducted by this API. The caller determines eligibility for VA settings.

The regulation takes precedence over omissions in the technical document.
Article 51a sets CSSR to one when liability PVBP is nonpositive; §9.3.2 of the
technical document omits that exception. Both calculation paths implement it.

## Basic curve

```@example eiopa
using FinanceModels
inputs = EIOPA.RFRInputs([0.02, 0.025, 0.03, 0.032], [1, 5, 10, 20];
    frequency=1, cra_bp=10, fsp=10, llfr_weights=[0.0, 0.0, 0.4, 0.6],
    ufr=Periodic(0.033, 1), α=0.11)
result = EIOPA.calculate(inputs)
discount(result.curves.basic, 60)
```

Supply complete, selected DLT observations at increasing positive integer tenors,
an observed FSP, and aligned weights summing to one (zero before FSP). Swap inputs
are nominal annual coupons, with `frequency` payments per year. Government inputs
are annual-effective zero rates (`instrument=:government`). CRA is in basis points;
UFR is a `Rate` with explicit compounding. Alpha is supplied, not calibrated: use
0.40 for SEK, 0.11 otherwise, or the applicable approved phasing-in parameter.

Supply FSP from the applicable DLT assessment. It is independent of whether its
LLFR weight is positive: the technical document's Table 2 includes assessed FSPs
with zero weight. Such inputs are accepted, and VA then has no effect on LLFR.
The demonstration workbook infers FSP from its first positive weight; that shortcut
is appropriate only when the assessed FSP has a positive weight.

`result.bootstrapped` retains all DLT maturities. The final basic curve retains it
through FSP; later observations influence LLFR only and need not be repriced by
the final extrapolated curve. Discounts are continuous at FSP; forwards may jump.

## CSSR and volatility adjustment

First construct the notional CSSR curve for liability revaluation, and apply the
same notional VA as a parallel spread shock to each asset's own valuation curve:

```@example eiopa
notional = EIOPA.calculate(inputs;
    va_inputs=EIOPA.VAInputs(rcs_currency=0.01))
asset_base = notional.curves.basic + Yield.Constant(Periodic(0.01, 1))
asset_notional = EIOPA.asset_spread_shock(asset_base,
    notional.diagnostics.notional_va)
discount(asset_notional, 60)
```

Article 51a(2)(c) increases asset spreads at **all maturities**, including beyond
FSP. Preserve each asset's existing spread and valuation assumptions; apply the
shock within the asset model when cashflows depend on spreads. Liabilities use
`notional.curves.cssr`, whose VA impact fades after FSP. The following example
revalues fixed asset and liability cashflows at year 60:

```@example eiopa
va = EIOPA.VAInputs(rcs_currency=0.01,
    valuations=(assets=100 * discount(asset_base, 60),
        assets_notional=100 * discount(asset_notional, 60),
        liabilities=200 * discount(notional.curves.basic, 60),
        liabilities_notional=200 * discount(notional.curves.cssr, 60)))
adjusted = EIOPA.calculate(inputs; va_inputs=va)
adjusted.diagnostics.total_va
```

Alternatively supply `sensitivities=(assets=..., liabilities=...)` in common PVBP
units, or a precomputed `cssr`. Only one CSSR source is accepted. When notional VA
is zero, supply sensitivities/CSSR rather than dividing valuation differences by
zero. Spreads and VA are decimal annual rates; negative VA is supported.

For macro VA set `use_macro_va=true`, with `rcs_country` and `debt_to_assets`.
For an undertaking-specific adjustment supply `usa` (at most 1.05). These are
mutually exclusive. Euro-area undertakings without USA must apply macro VA
(TD §9.3.4), except the peer-country cases prescribed to have zero macro VA.
The default `false` does not infer currency, domicile or eligibility; the caller
must select the applicable treatment.

The notional VA is 85% of the unadjusted published currency RCS, even when USA is
selected, matching workbook `CSSR Calculation!H10`. USA affects the final currency
VA. TD §9.3.5 redefines currency RCS using USA without explicitly resolving its
use in the notional step; this API follows the workbook convention.

VA adds `log1p(VA)` to pre-FSP forwards and changes only the first weighted LLFR
term. The UFR remains fixed; shifting the entire extrapolated curve differs.

## Precise curves and publication outputs

```@example eiopa
shocks = EIOPA.ShockParameters()
stressed = EIOPA.calculate(inputs; va_inputs=va, shock_parameters=shocks)
published = EIOPA.publication(inputs; va_inputs=va, shock_parameters=shocks)
(published.maturities[60], published.rates.basic[60], published.rates.up[60])
```

`calculate` is the full-precision analytical path for sensitivities and AD,
including unrounded VA. Regulatory technical provisions use published basic
rates and whole-basis-point VA (TD §§10.1.2–10.1.3); use `publication` for that
convention. Unrounded VA can differ by up to half a basis point before FSP.
`calculate` returns `bootstrapped`, `curves` and `diagnostics`.
Curve names are `basic`, `cssr`, `va`, `up`, `down`, `up_va` and
`down_va`; unavailable optional outputs are `nothing`.

Diagnostics include all basic/adjusted LLFRs, notional/currency/macro/total VA,
CSSR, the country weight (`country_weight`, ω), and asset/liability PVBP.
For valuation inputs PVBP is the value reduction divided by the notional decimal
VA; for sensitivity inputs the supplied common units are retained. Adjusted
downward LLFR diagnostics refer to the extrapolation before the final floor.

`publication` rounds basic annual zero rates to five decimal places before the
CSSR/VA stage and notional/total VA to whole basis points. Shocks use unrounded
bootstrapped annual rates, their own LLFRs and annual-effective UFR shifts of ±15
basis points, retaining the undertaking's total VA. Final rates are rounded to
five decimals, then downward floors apply without another rounding pass. A binding
floor can retain more digits. Publication outputs contain numeric annual-effective
`rates`, defaulting to years 1–150; `maturities` accepts another increasing positive
integer grid. Decimal ties follow Excel's 15-significant-digit arithmetic and
round away from zero.

For valuations supplied to `publication`, construct liability curves from its
published annual rates and use its rounded notional VA for the parallel asset
shock. For example:

```@example eiopa
published_notional = EIOPA.publication(inputs;
    va_inputs=EIOPA.VAInputs(rcs_currency=0.01))
published_cssr = ZeroRateCurve(Periodic.(published_notional.rates.cssr, 1),
    published_notional.maturities, Spline.ConstantForward())
published_asset_notional = EIOPA.asset_spread_shock(asset_base,
    published_notional.diagnostics.notional_va)
discount(published_cssr, 60)
```

Use precise CSSR curves and notional VA for valuations supplied to `calculate`.
First-order AD is supported in precise calculations
away from clipping/floor boundaries; publication rejects dual-number inputs.
Maturities and FSP are fixed. Precise shocked curves use constant-forward
interpolation between annual pre-FSP nodes. Custom `ShockParameters` interpolate
parameters linearly between their tenors and hold endpoints flat outside the grid.

Article 46(1c)(a) and the workbook calculate the first shocked LLFR term from the
preceding DLT observation to FSP. On the basic constant-forward curve that equals
the one-year forward into FSP in TD §8.5.6; after maturity-dependent shocks those
intervals can differ. Both paths use the workbook's DLT-interval convention for
shocked LLFRs, while retaining every annual shocked node in the pre-FSP curve.

## Reusable curve components

`Spline.ConstantForward()` works with `ZeroRateCurve` and
`fit(Spline.ConstantForward(), quotes, Fit.Bootstrap())`. It interpolates cumulative
log-discounts linearly, unlike `Spline.Linear()` (linear zero rates). It supports
the knot/reconstruction interfaces and usual tail policies except `:extension`.

`Yield.ExponentialForward(base; start, initial_forward, ufr, α)` preserves a curve
through `start`, then converges its instantaneous forward exponentially to UFR.
Both assumptions are explicit `Rate`s. Existing curve defaults remain unchanged.
The FSP log-discount is cached at construction. Reconstruction with `Accessors.@set`
validates updated assumptions and refreshes that cache.

## API reference

```@autodocs
Modules = [FinanceModels.EIOPA]
```
