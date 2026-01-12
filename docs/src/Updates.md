# FinanceModels.jl

## Design

- **Contracts** represent insturments that ultimately obligate a payment of cashflows, which may or may not be scenario dependant.
- **Quotes** are observed or reference prices that may be used to `fit` models.
- **Models** are the combination of **assumptions** and **logic** that can then be used to realize the assumed cashflows that arise from a contract.

## Motivation

FinanceModels.jl is the evolution of Yields.jl. Yields.jl was originally designed for very nice usage of term structures of yield curves, but three aspects held it back:

1. The design was very oriented towards interest rates, and it was awkward to stick, e.g. volatility models into a package called Yields.jl
2. The API for contructing curves was inconsistent because there are different ways to construct a given curve and the inputs to constructing a simple bootstrapped curve with a spline through given yields vs a best-fit of a variety of instrumnets was simply a different paradigm.
3. There was a lack of ability to even express some types of contracts that are useful for model-fitting or modeling in general.

## TODOs

- `bond.frequency.frequency` is awkward
- Core contracts:
  - Composite contact (e.g. Fixed + Float -> Swap)
  - Forward contact
  - Derivatives?
  - distinguish between clean and dirty prices
- Projections
  - Everythign is currently coerced to a F64/F64 Cashflow, but would like to be flexible with amount and timepoints
- How to integrate Dates?
- Core methods:
  - port Yields.jl methods

- Ergonomics
    -

- Package design:
  - promote `pv` to FinanceCore given it's utility here
  - promote `Cashflow` up to FC
