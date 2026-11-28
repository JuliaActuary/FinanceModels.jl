# [Fitting a Curve to a Bond Portfolio](@id bond-portfolio-fit)

Suppose you have prices for 1,000 fixed-rate bonds in a credit index and want one curve to
value another bond from the same market. Bonds with the same maturity can trade at different
yields. Fit a small set of knot rates to all the prices, then check how much the choice of
weights matters for the bond you want to value.

## Build the quotes

Each observation pairs a full price per unit of face value with the contract that pays the
cashflows. Here every bond pays semiannually on the same half-year grid, starting six months
from valuation. A coupon of `0.05` pays `0.025` per half year. For prices quoted per 100 of
face, divide by 100. If the source reports a clean price, add accrued interest in the same
units first and use the remaining payment schedule as of the price date.

The example generates its own prices so you can run it without a data file. Each bond is
priced on a common curve plus a bond-specific spread with a standard deviation of 30 basis
points. For observed data, replace that pricing step with `Quote(full_price, bond)`.

```@example bond_portfolio
using FinanceModels, Random

rng = Xoshiro(42)
bonds = [
    Bond.Fixed(rand(rng, 4:14) / 200, Periodic(2), rand(rng, 1:40) / 2)
    for _ in 1:1_000
]
reference = Yield.NelsonSiegel(2.0, 0.045, -0.018, 0.012)
quotes = map(bonds) do bond
    spread = Yield.Constant(Continuous(0.003 * randn(rng)))
    Quote(pv(reference + spread, bond), bond)
end

(bonds = length(quotes), distinct_maturities = length(unique(maturity.(quotes))))
```

## Fit eight knot rates to all 1,000 prices

Supply a curve with its own knots. Passing a spline alone would place a knot at each quote
maturity, and this portfolio has repeated maturities. The grid below spans the observed
half-year to 20-year range, with closer spacing at the short end. Linear interpolation makes
each zero rate between two knots depend only on those neighbors. The fit uses all the
coupon and principal payments between knots.

```@example bond_portfolio
knots = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0]
grid = ZeroRateCurve(fill(0.04, length(knots)), knots, Spline.Linear())
price_fit = fit(grid, quotes, Fit.LeastSquares())

collect(zip(knot_tenors(price_fit), round.(100 .* rate.(knot_rates(price_fit)); digits = 3)))
```

Each pair gives a tenor in years and its fitted continuously compounded zero rate in
percent. The input curve supplies the grid and interpolation method. FinanceModels starts
the refit from its standard knot rates. It minimizes the sum of squared differences between
model and market prices. A one-point price miss per 100 of face counts equally for a
two-year and a twenty-year bond.

## Give comparable yield misses comparable weight

A one-basis-point yield change generally moves the longer bond's price more. To target yield
errors, measure how much each bond's price changes for one basis point at its own market
yield. Use the same nominal, semiannually compounded yield convention for every bond.

```@example bond_portfolio
ytm(q) = rate(convert(Periodic(2), internal_rate_of_return(q)))
market_yields = ytm.(quotes)

price_per_bp = map(quotes, market_yields) do q, y
    pv_at_lower_yield = pv(Yield.Constant(Periodic(y - 0.0001, 2)), q.instrument)
    pv_at_higher_yield = pv(Yield.Constant(Periodic(y + 0.0001, 2)), q.instrument)
    (pv_at_lower_yield - pv_at_higher_yield) / 2
end
weights = 1 ./ price_per_bp .^ 2
yield_fit = fit(grid, quotes, Fit.LeastSquares(; weights))
nothing # hide
```

Dividing a small pricing error by `price_per_bp` expresses its magnitude in approximate
yield basis points. In the squared-error objective, each squared price error therefore gets
weight `1 / price_per_bp^2`. Compute the weights once from the market observations so they
stay fixed during calibration.

## Check both fits in the units you care about

Reprice the original contracts and solve for their yields at the model prices. Compare the
price error, its first-order yield equivalent, and the actual yield difference. Each column
reports a root-mean-square error over the same 1,000 observations.

```@example bond_portfolio
rms(x) = sqrt(sum(abs2, x) / length(x))
function fit_errors(curve, quotes, market_yields, price_per_bp)
    model_prices = [pv(curve, q.instrument) for q in quotes]
    price_gaps = model_prices .- [q.price for q in quotes]
    model_yields = [ytm(Quote(p, q.instrument)) for (p, q) in zip(model_prices, quotes)]
    return (
        price_per_100 = 100 * rms(price_gaps),
        approx_yield_bp = rms(price_gaps ./ price_per_bp),
        yield_bp = 10_000 * rms(model_yields .- market_yields),
    )
end

price_errors = fit_errors(price_fit, quotes, market_yields, price_per_bp)
yield_errors = fit_errors(yield_fit, quotes, market_yields, price_per_bp)
using Markdown # hide
Markdown.parse("| Fit | Price points per 100 | Approximate yield (bp) | Actual yield (bp) |\n|:---|---:|---:|---:|\n" * # hide
    join(["| $name | $(round(e.price_per_100; digits=3)) | $(round(e.approx_yield_bp; digits=2)) | $(round(e.yield_bp; digits=2)) |" # hide
          for (name, e) in (("Equal price weights", price_errors), ("Inverse squared price sensitivity", yield_errors))], "\n")) # hide
```

Equal price weights give the smaller price error. The second fit gives short bonds more
influence and reduces the approximate yield error. Both fits miss actual yields by about
30 basis points, close to the spread dispersion built into the sample. The small difference
between the last two columns comes from the curvature of price as a function of yield.

To inspect the bond with the largest price miss under the yield-weighted fit, keep the
sign of the difference. A positive gap means the curve values the bond above its market price.

```@example bond_portfolio
gaps = [100 * (pv(yield_fit, q.instrument) - q.price) for q in quotes]
i = argmax(abs.(gaps))
(bond = i, maturity_years = maturity(quotes[i]), market_per_100 = round(100 * quotes[i].price; digits = 3),
 gap_per_100 = round(gaps[i]; digits = 3))
```

On observed data, check the coupon, remaining payments, and quoted price for that bond.
If many bonds in the same maturity range miss in the same direction, inspect that part
of the curve before adding another knot.

## Value another bond

For a new seven-and-a-half-year bond with a 4.75% coupon, the two curves give these full prices
per 100 of face.

```@example bond_portfolio
new_bond = Bond.Fixed(0.0475, Periodic(2), 7.5)
(price_fit = round(100 * pv(price_fit, new_bond); digits = 3),
 yield_fit = round(100 * pv(yield_fit, new_bond); digits = 3))
```

The difference shows how the calibration objective carries through to a price you would
actually use. To measure how that price responds to changes in individual market quotes,
see [Sensitivities Through Calibration](@ref).
