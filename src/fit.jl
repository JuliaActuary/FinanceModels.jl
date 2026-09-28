module Fit

    abstract type FitMethod end


    """
        Fit.Loss(function)

    `function` should be a loss measure, such as `x->x^2` or `x->abs(x)`. This is used by the optimization algorithm in `fit` to determine optimal parameters as defined by this loss function.

    A subtype of FitMethod.

    # Examples
    ```julia-repl
    julia> mod0 = Yield.Constant();

    julia> quotes = ZCBPrice([0.9, 0.8, 0.7,0.6]);

    julia> fit(mod0,quotes,Fit.Loss(x->x^2))
    FinanceModels.Yield.Constant{Rate{Float64, Periodic}}(Periodic(0.12822921882254446, 1))
    ```

    (With `UnicodePlots` loaded, fitted yield models display as a zero-rate chart instead.)
    """
    struct Loss{T} <: FitMethod
        fn::T
    end

    """
        Bootstrap()

    A singleton type passed to `fit` to bootstrap a spline curve one quote at a time.
    Each step solves for the zero rate at the next quote maturity to match its price.

    Supports `Spline.Linear()` (equivalently `Spline.PolynomialSpline(1)` or
    `Spline.BSpline(1)`). Bootstrapping requires that adding a knot on the right
    leaves every earlier curve segment unchanged, so earlier quotes stay exactly
    priced. Linear interpolation of zero rates has this property; quadratic,
    cubic, higher-order B-spline, PCHIP, and Akima interpolation do not, because
    a later knot changes the shape of earlier segments. Fit those strategies
    across the complete quote set with `Fit.Loss(x -> x^2)`, and fit a monotone
    convex curve with `fit(Spline.MonotoneConvex(), quotes)`.

    After solving, every quote is repriced on the returned curve; a residual
    beyond root-finder precision throws an `ArgumentError` naming the quote.

    A subtype of FitMethod.

    # Examples

    ```julia
    quotes = ZCBPrice([0.99, 0.97, 0.94])
    curve = fit(Spline.Linear(), quotes, Fit.Bootstrap())
    discount(curve, 2) ≈ 0.97 # true
    ```
    """
    struct Bootstrap <: FitMethod
        # spline method
    end


end

"""
    __default_optic(model)

 Returns the variables to optimize over for the given model. This is an optic/lens specifying which parameters of the model can vary. See extended help for more.
An optic argument is a tuple of optic => interval pairs specifying which model parameters to optimize and their bounds.

# Examples

We might have a model as follows where we want `fit` to optize parameters `a` and `b`:

```julia
struct MyModel <:FinanceModels.AbstractModel
        a 
        b 
end

__default_optic(m::MyModel) = (
    @optic(_.a) => 0.0 .. 100.0,
    @optic(_.b) => -10.0 .. 10.0,
)
```

# Extended help

An arbitrarily complex model may be the object we intend to fit - how does `fit` know what free variables are able to be solved for within the given model?
`variables` is a tuple of optic => interval pairs. What does this mean?
- An optic (or "lens") is a way to define an accessor to a given object. Example:

```julia-repl
julia> using Accessors, AccessibleModels, IntervalSets

julia> obj = (a = "AA", b = "BB");

julia> lens = @optic _.a
(@optic _.a)

julia> lens(obj)
"AA"
```
An optic argument is a tuple of optic => interval pairs. For example, we might have a model as follows where we want 
`fit` to optize parameters `a` and `b`:

```julia
struct MyModel <:FinanceModels.AbstractModel
        a 
        b 
end

__default_optic(m::MyModel) = (
    @optic(_.a) => 0.0 .. 100.0,
    @optic(_.b) => -10.0 .. 10.0,
)
```
In this way, fit know which arbitrary parameters in a given object may be modified. Technically, we are not modifying the immutable `MyModel`, but instead efficiently creating a new instance. This is enabled by [AccessibleModels.jl](https://github.com/JuliaAPlavin/AccessibleModels.jl).

Note that not all optimization algorithms want a bounded interval. In that case, simply leave off the paired range. The prior example would then become:

```julia
__default_optic(m::MyModel) = (
    (@optic(_.a),),
    (@optic(_.b),),
)
```
```

    

"""
__default_optic(m::Yield.Constant) = ((@optic(_.rate.continuous_value) => -1.0 .. 1.0),)
"""
    KnotRatesOptic()

Batch optic over every knot rate of a `Yield.AbstractInterpolatedZeroCurve`. `getall` returns the
rates as a tuple; `setall` replaces them all with **one** [`reconstruct`](@ref) of the curve (one
interpolant build per optimizer candidate). This is the default fitting variable for those
curves: a tuple of per-knot optics would rebuild the curve once per knot per candidate, i.e.
O(n²) work and allocation in the number of knots.
"""
struct KnotRatesOptic end
Accessors.OpticStyle(::Type{KnotRatesOptic}) = Accessors.ModifyBased()
Accessors.getall(m::Yield.AbstractInterpolatedZeroCurve, ::KnotRatesOptic) = Tuple(Yield.knot_rates(m))
Accessors.setall(m::Yield.AbstractInterpolatedZeroCurve, ::KnotRatesOptic, vals) = Yield.reconstruct(m; rates = vals)
Accessors.modify(f, m::Yield.AbstractInterpolatedZeroCurve, o::KnotRatesOptic) =
    Accessors.setall(m, o, map(f, Accessors.getall(m, o)))

__default_optic(m::Yield.AbstractInterpolatedZeroCurve) = ((KnotRatesOptic() => -1.0 .. 1.0),)

__default_optic(m::Yield.NelsonSiegel) = (
    @optic(_.τ₁) => 0.0 .. 100.0,
    @optic(_.β₀) => -10.0 .. 10.0,
    @optic(_.β₁) => -10.0 .. 10.0,
    @optic(_.β₂) => -10.0 .. 10.0,
)
__default_optic(m::Yield.NelsonSiegelSvensson) = (
    @optic(_.τ₁) => 0.0 .. 100.0,
    @optic(_.τ₂) => 0.0 .. 100.0,
    @optic(_.β₀) => -10.0 .. 10.0,
    @optic(_.β₁) => -10.0 .. 10.0,
    @optic(_.β₂) => -10.0 .. 10.0,
    @optic(_.β₃) => -10.0 .. 10.0,
)
__default_optic(m::Yield.CairnsPritchard) = (
    @optic(_.c₁) => 0.001 .. 10.0,
    @optic(_.c₂) => 0.001 .. 10.0,
    @optic(_.b₀) => -1.0 .. 1.0,
    @optic(_.b₁) => -10.0 .. 10.0,
    @optic(_.b₂) => -10.0 .. 10.0,
)
__default_optic(m::Yield.CairnsPritchardExtended) = (
    @optic(_.c₁) => 0.001 .. 10.0,
    @optic(_.c₂) => 0.001 .. 10.0,
    @optic(_.c₃) => 0.001 .. 10.0,
    @optic(_.b₀) => -1.0 .. 1.0,
    @optic(_.b₁) => -10.0 .. 10.0,
    @optic(_.b₂) => -10.0 .. 10.0,
    @optic(_.b₃) => -10.0 .. 10.0,
)
__default_optic(m::Equity.BlackScholesMerton{T, U, V}) where {T, U, V <: Volatility.Constant} = ((@optic(_.σ.σ) => 0.0 .. 10.0),)
__default_optic(m::Volatility.Constant) = ((@optic(_.σ) => 0.0 .. 10.0),)
__default_optic(m::ShortRate.Vasicek) = (
    @optic(_.a) => 0.0 .. 5.0,
    @optic(_.b) => -0.1 .. 0.5,
    @optic(_.σ) => 0.0 .. 1.0,
    @optic(_.initial.continuous_value) => -0.05 .. 0.2,
)
__default_optic(m::ShortRate.CoxIngersollRoss) = (
    @optic(_.a) => 0.0 .. 5.0,
    @optic(_.b) => 0.0 .. 0.5,
    @optic(_.σ) => 0.0 .. 1.0,
    @optic(_.initial.continuous_value) => 0.0 .. 0.2,
)
__default_optic(m::ShortRate.HullWhite) = (
    @optic(_.a) => 0.0 .. 5.0,
    @optic(_.σ) => 0.0 .. 1.0,
)
# FX.Forwards: the free variables live on the base-currency (`foreign`) curve — spot and
# the domestic curve are calibration inputs — so compose the foreign curve's own optics
# through the `foreign` field.
__default_optic(m::FX.Forwards) = map(o -> __fx_foreign_optic(o), __default_optic(m.foreign))
__fx_foreign_optic(o::Base.Pair) = Accessors.opcompose(@optic(_.foreign), o.first) => o.second
__fx_foreign_optic(o::Tuple) = (Accessors.opcompose(@optic(_.foreign), only(o)),)
__fx_foreign_optic(o) = Accessors.opcompose(@optic(_.foreign), o)


__default_optim(m) = OptimizationOptimJL.LBFGS()
__default_optim(m::T) where {T <: Spline.SplineCurve} = OptimizationOptimJL.Newton()
__default_optim(::Spline.MonotoneConvex) = OptimizationOptimJL.LBFGS()

__default_loss(m) = Fit.Loss(x -> x^2)

# One AD backend for every `fit` loss function. `SecondOrder` serves both first-order
# optimizers (LBFGS/Fminbox use the gradient, via the inner backend) and second-order
# ones (Newton/IPNewton use the Hessian), so OptimizationBase never auto-promotes a
# first-order declaration and warns. First-order paths never instantiate the Hessian,
# so declaring it costs them nothing.
const __FIT_ADTYPE = DifferentiationInterface.SecondOrder(AutoForwardDiff(), AutoForwardDiff())

# The fitting loss: `loss_method.fn` summed over the price residuals of `quotes` under `model`.
__quote_loss(model, loss_method, quotes) = mapreduce(+, quotes) do q
    loss_method.fn(present_value(model, q.instrument) - q.price)
end

# A knot-curve fit's trial curve (optimizer candidate, bootstrap step, calibration Jacobian): knot
# rates `z` over tenors validated up front, built without a copy or a finite-rate check, so a
# non-finite candidate mid-search is a bad residual rather than an exception.
__trial_curve(spline, z, tenors, extrapolation) =
    Yield.__build(spline, Yield.KnotGrid(Yield.Unchecked(), z, tenors); extrapolation)

"""
    FitConvergenceError(retcode, msg)

Thrown by [`fit`](@ref FinanceModels.fit) when the optimizer does not report a successful solve. `retcode` is
the solver's `SciMLBase.ReturnCode` (for example `ReturnCode.MaxIters` or
`ReturnCode.Failure`) and `msg` describes the failed fit.

A failed solve can leave the parameters at the starting guess, so `fit` throws rather than
return an unfitted model. Address the cause and fit again: a different starting model, a
different `optimizer`, tighter or longer `solve_kwargs` (for example
`solve_kwargs = (; maxiters = 10_000)`), or quotes that the model can fit. Choose the model
deliberately; switching to another interpolation method on failure changes the curve.
"""
struct FitConvergenceError{R} <: Exception
    retcode::R
    msg::String
end

function Base.showerror(io::IO, e::FitConvergenceError)
    return print(io, "FitConvergenceError: ", e.msg, " (optimizer return code ", e.retcode, ")")
end

# Minimize `loss(x, quotes)` from `x0` and return the solution. The quotes are the problem's data
# (capturing them in `loss` made a Nelson-Siegel fit about 8% slower). Every optimizer-backed fit
# either returns a successful solution or throws: a failed solve may leave the solution equal to
# `x0`, and a model built from it would make optimizer failure look like a valid (but unfitted)
# result.
function __minimize(loss, x0, quotes, optimizer, solve_kwargs; lb = nothing, ub = nothing)
    prob = Optimization.OptimizationProblem(Optimization.OptimizationFunction(loss, __FIT_ADTYPE), x0, quotes; lb, ub)
    sol = Optimization.solve(prob, optimizer; solve_kwargs...)
    Optimization.SciMLBase.successful_retcode(sol.retcode) || throw(
        FitConvergenceError(sol.retcode, "model fitting did not converge")
    )
    return sol.u
end

"""
    fit(
        model, 
        quotes, 
        method=Fit.Loss(x -> x^2);
        variables=__default_optic(model), 
        optimizer=__default_optim(model),
        solve_kwargs=(;)
        )

Fit a model to a collection of quotes using a loss function and optimization method.

## Arguments
- `model`: The initial model to fit, which is generally an instantiated but un-optimized model.
- `quotes`: A collection of quotes to fit the model to.
- `method::F=Fit.Loss(x -> x^2)`: The loss function to use for fitting the model. Defaults to the squared loss function. 
  - `method` can also be `Bootstrap()` with `Spline.Linear()`. Other interpolation strategies require a full-curve `Fit.Loss`.
- `variables=__default_optic(model)`: The variables to optimize over. This is a tuple of optic => interval pairs specifying which parameters of the model can vary. See extended help for more.
- `optimizer=__default_optim(model)`: The optimization algorithm to use. The default optimization for a given model is `LBFGS()` from Optim.jl (via OptimizationOptimJL), a quasi-Newton method with automatic differentiation via ForwardDiff. See extended help for more on customizing the solver.
- `solve_kwargs=(;)`: Keyword arguments passed to `Optimization.solve` with the optimizer, such as
  `(; maxiters = 10_000, abstol = 1e-12)` or Optim.jl's `g_tol`. Use them to tighten a loss fit.
- `extrapolation=:flat_forward`: For `Spline.SplineCurve` fits (including
  `Spline.MonotoneConvex()`), the long-end policy beyond the last knot (see
  [`Yield.Spline`](@ref FinanceModels.Yield.Spline)). Also accepts `:flat_zero`, `:linear`,
  `Yield.FlatForwardAt(rate)`, and (except for MonotoneConvex) `:extension`.

Fitting a `Spline.SplineCurve` places one knot at each quote maturity and returns the same
curve type as [`ZeroRateCurve`](@ref): a `Yield.MonotoneConvex` for `Spline.MonotoneConvex()`
(whose default optimizer is `LBFGS()`) and a `Yield.Spline` otherwise (default `Newton()`).
`Fit.Bootstrap()` accepts `Spline.Linear()` only. Fitting an existing knot curve varies its knot
rates and preserves its tenors, method and `extrapolation`. With fewer quotes than a polynomial or
B-spline needs for its order, the order is reduced to one less than the number of knots (see
[`Spline.PolynomialSpline`](@ref FinanceModels.Spline.PolynomialSpline)): `fit(Spline.Cubic(), quotes)`
with two quotes returns a linear curve.

The optimization routine will then attempt to modify parameters of `model` to best fit the quoted prices of the contracts underlying the `quotes` by calling `present_value(model,contract)`. The optimization will minimize the loss function specified within `Fit.Loss(...)`. 

Different types of quotes are appropriate for different kinds of models. For example, if you try to value a set of equity `Option.EuroCall`s with a `Yield.Constant`, you will get an error because the `present_value(m<:Yield.Constant,o<:Option.EuroCall)` is not defined.

## Returns
- The fitted model.

## Differentiating through a fit

Spline fits (`fit(spline, quotes)` and `fit(Spline.Linear(), quotes, Fit.Bootstrap())`,
including `FX.Forwards` with a spline foreign curve) are differentiable with ForwardDiff: when
quote prices, rates, or cashflow amounts are dual numbers, the fitted knot rates carry the exact
first-order derivatives of the calibration. They come from the implicit function theorem at
the fitted curve, not from the solver's iterations, and their values are the primal fit exactly.

```julia
tenors = [1.0, 2.0, 5.0, 10.0]
value(rates) = pv(fit(Spline.Linear(), OISYield.(rates, tenors), Fit.Bootstrap()), cashflows)
ForwardDiff.gradient(value, [0.03, 0.032, 0.035, 0.037])   # ∂value/∂quoted rates
```

- The derivatives are first order, with respect to the quotes passed to `fit`. Nested
  (higher-order) dual numbers, and dual maturities or cashflow times, throw.
- The fit must reprice its quotes. Bootstrap does; a loss fit is differentiated only if the
  largest absolute component of one Newton correction of its knot rates towards the exact fit
  is at most `1e-6`. That correction is a local estimate of the fit's error, not a guaranteed
  distance to the exact solution.
- A fitted `Spline.MonotoneConvex()`, `Spline.PCHIP()`, or `Spline.Akima()` curve that lies on a
  kink of its interpolation (flat quotes, for example) has no derivative there and throws.
- Other models' `fit`s throw an `ArgumentError` when given dual quotes.

See [Sensitivities Through Calibration](@ref) for details.

# Examples
```julia-repl
julia> model = Yield.Constant();

julia> quotes = ZCBPrice([0.9, 0.8, 0.7,0.6]);

julia> fit(model,quotes)
FinanceModels.Yield.Constant{Rate{Float64, Periodic}}(Periodic(0.12822921882254446, 1))
```

(With `UnicodePlots` loaded, fitted yield models display as a zero-rate chart instead.)

# Extended help

## Customizing the Solver

The default solver is `LBFGS()` from Optim.jl (via OptimizationOptimJL). This is a quasi-Newton method that uses automatic differentiation (ForwardDiff) to compute gradients efficiently.
 - Any solver from OptimizationOptimJL can be used, e.g. `fit(...; optimizer=OptimizationOptimJL.Newton())` or `fit(...; optimizer=OptimizationOptimJL.NelderMead())`.
 - Solver settings go in `solve_kwargs`, e.g. `fit(...; solve_kwargs = (; maxiters = 10_000, g_tol = 1e-12))`.
 - More documentation is available from the upstream packages:
   - [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/)
   - [Optimization.jl](https://docs.sciml.ai/Optimization/stable/)
   - [AccessibleModels.jl](https://github.com/JuliaAPlavin/AccessibleModels.jl)

## Defining the variables

An arbitrarily complex model may be the object we intend to fit - how does `fit` know what free variables are able to be solved for within the given model?
`variables` is a tuple of optic => interval pairs. What does this mean?
- An optic (or "lens") is a way to define an accessor to a given object. Example:

```julia-repl
julia> using Accessors, AccessibleModels, IntervalSets

julia> obj = (a = "AA", b = "BB");

julia> lens = @optic _.a
(@optic _.a)

julia> lens(obj)
"AA"
```
An optic argument is a tuple of optic => interval pairs. For example, we might have a model as follows where we want 
`fit` to optimize parameters `a` and `b`:

```julia
struct MyModel <:FinanceModels.AbstractModel
     a 
     b 
end

__default_optic(m::MyModel) = (
    @optic(_.a) => 0.0 .. 100.0,
    @optic(_.b) => -10.0 .. 10.0,
)
```
In this way, fit know which arbitrary parameters in a given object may be modified. Technically, we are not modifying the immutable `MyModel`, but instead efficiently creating a new instance. This is enabled by [AccessibleModels.jl](https://github.com/JuliaAPlavin/AccessibleModels.jl).

Note that not all optimization algorithms want a bounded interval. In that case, simply leave off the paired range. The prior example would then become:

```julia
__default_optic(m::MyModel) = (
    (@optic(_.a),),
    (@optic(_.b),),
)
```
```


## Additional Examples

See the tutorials in the package documentation for FinanceModels.jl or the docstrings of FinanceModels.jl's available model types.
"""
function fit(
        mod0,
        quotes,
        method::F = __default_loss(mod0);
        variables = __default_optic(mod0),
        optimizer = __default_optim(mod0),
        solve_kwargs = (;)
    ) where
    {F <: Fit.Loss}
    __quotes_carry_ad(quotes) && throw(
        ArgumentError(
            "fit does not differentiate through the calibration of $(nameof(typeof(mod0))) models; " *
                "spline fits such as fit(Spline.Linear(), quotes) do. To differentiate a valuation with " *
                "respect to a fitted curve's knot rates, use reconstruct(curve; rates = dual_rates)."
        )
    )
    # AccessibleModels parameterizes the model: the transformed vector the optimizer moves, its
    # bounds, and the model rebuilt from a candidate.
    params = AccessibleModel(mod0, variables)
    build(x) = AccessibleModels.from_transformed(x, params)
    x0 = collect(AccessibleModels.transformed_vec(params))
    bounds = AccessibleModels.transformed_bounds(params)
    lb = haskey(bounds, :lb) ? collect(bounds.lb) : nothing
    ub = haskey(bounds, :ub) ? collect(bounds.ub) : nothing
    # Ensure x0 is strictly interior to avoid Fminbox boundary warnings
    # and to avoid degenerate starting points (e.g., σ=0 for Black-Scholes)
    if lb !== nothing && ub !== nothing
        for i in eachindex(x0)
            if x0[i] <= lb[i] || x0[i] >= ub[i]
                x0[i] = (lb[i] + ub[i]) / 2
            end
        end
    end
    # `__FIT_ADTYPE` is SecondOrder so a bounds-compatible second-order optimizer
    # (e.g. `IPNewton()`) finds a Hessian; the default `Fminbox(LBFGS())` uses only the
    # gradient (via the inner `AutoForwardDiff`) and is unaffected.
    loss(x, qs) = convert(eltype(x), __quote_loss(build(x), method, qs))
    return build(__minimize(loss, x0, quotes, optimizer, solve_kwargs; lb, ub))

end

# Flat knot rates are singular under ForwardDiff for interpolants whose slope formula
# contains divided differences (notably PCHIP and Akima), so each method starts from a
# small slope: the DataInterpolations methods near a flat 5%, MonotoneConvex from 1% to 5%
# (a single knot starts at the midpoint).
# PCHIP and Akima also switch formulas where adjacent secant slopes are equal, which a seed
# linear in the knot number gives on an evenly spaced grid, and where PCHIP's end slopes are 0
# or three times the end secant. They start from a curve that is strictly increasing and
# strictly concave in the tenors themselves, which stays away from every switch.
__curve_fit_seed(n, lo, hi) = n <= 1 ? fill((lo + hi) / 2, n) : collect(range(lo, hi; length = n))
__knot_fit_seed(::Spline.SplineCurve, tenors) = __curve_fit_seed(length(tenors), 0.049, 0.051)
__knot_fit_seed(::Spline.MonotoneConvex, tenors) = __curve_fit_seed(length(tenors), 0.01, 0.05)
__knot_fit_seed(::Union{Spline.PCHIP, Spline.Akima}, tenors) = 0.05 .- 0.01 .* exp.(-tenors ./ (1 + maximum(tenors; init = 0)))

function fit(
        mod0::T, quotes; optimizer = __default_optim(mod0), solve_kwargs = (;),
        extrapolation = :flat_forward
    ) where {T <: Spline.SplineCurve}
    return fit(mod0, quotes, __default_loss(mod0); optimizer, solve_kwargs, extrapolation)
end

# Every interpolation method, including `Spline.MonotoneConvex()`, fits through this one
# method: knots at the sorted quote maturities, one trial curve per candidate built by the same
# `Yield.__build` as `ZeroRateCurve`, and a validated `ZeroRateCurve` as the result.
function fit(
        mod0::T, quotes, method::F; optimizer = __default_optim(mod0), solve_kwargs = (;),
        extrapolation = :flat_forward
    ) where {T <: Spline.SplineCurve, F <: Fit.Loss}
    # The solve runs on primal quotes; dual numbers in the quotes or the extrapolation policy
    # are propagated afterwards by the implicit function theorem (see `__implicit_knot_curve`).
    quotes, primal_quotes = __calibration_quotes(quotes)
    primal_extrapolation = __primal_extrapolation(extrapolation)
    # Validate the policy and the knot grid once, up front (duplicate maturities, too few
    # knots for the interpolant, `:extension` with MonotoneConvex, …) with the same errors as
    # direct construction; trial curves reuse them.
    tenors = sort!(maturity.(primal_quotes))
    grid0 = Yield.KnotGrid(__knot_fit_seed(mod0, tenors), tenors, mod0; who = "fit($(mod0))")
    __check_primal_quotes(Yield.__build(mod0, grid0; extrapolation = primal_extrapolation), primal_quotes)
    loss(u, qs) = __quote_loss(__trial_curve(mod0, u, grid0.tenors, primal_extrapolation), method, qs)
    rates = __minimize(loss, grid0.rates, primal_quotes, optimizer, solve_kwargs)
    curve = Yield.ZeroRateCurve(rates, grid0.tenors, mod0; extrapolation = primal_extrapolation)   # public result: validated (finite rates)
    return __implicit_knot_curve(curve, quotes, primal_quotes, extrapolation)
end

# FX.Forwards with a spline placeholder as the foreign curve: given spot, the domestic
# curve, and market quotes, each quote reduces to an equivalent quote on the
# base-currency discount curve — closed-form implied zero-coupon quotes for outright
# forwards (DF_f(t) = (K·DF_d(t) + price)/spot, see `FX.implied_zcb_quotes`) and par
# cashflow strips for basis swaps (see `FX.ParBasisSwap`) — so curve construction
# reduces to fitting the spline through those quotes; no optimizer runs over the FX
# model itself, and one quote set may mix both instrument types. The two dispatch
# methods are deliberately separate (not a `Union`) so neither is ambiguous against the
# generic optic-based `fit(mod0, quotes, ::Fit.Loss)`. Keywords (`optimizer`, `solve_kwargs`,
# `extrapolation`) go to the foreign-curve fit unchanged.
function __fit_fx_via_implied(mod0, quotes, method; kwargs...)
    implied = map(q -> FX.__implied_foreign_quote(mod0, q), quotes)
    return @set mod0.foreign = fit(mod0.foreign, implied, method; kwargs...)
end
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.Bootstrap; kwargs...) where {P, S, D, F <: Spline.SplineCurve} = __fit_fx_via_implied(mod0, quotes, method; kwargs...)
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.Loss; kwargs...) where {P, S, D, F <: Spline.SplineCurve} = __fit_fx_via_implied(mod0, quotes, method; kwargs...)

# Appending a knot must preserve every earlier segment, not merely the earlier
# knot values. In particular, "local" smooth splines need not have this property.
__supports_bootstrap(::Spline.SplineCurve) = false
__supports_bootstrap(s::Union{Spline.PolynomialSpline, Spline.BSpline}) = s.order == 1

__bootstrap_alternative(::Spline.SplineCurve) = "fit this strategy across all quotes with Fit.Loss(x -> x^2)"
__bootstrap_alternative(::Spline.MonotoneConvex) = "fit a monotone convex curve to all quotes with fit(Spline.MonotoneConvex(), quotes)"

function fit(
        mod0::T, quotes, method::Fit.Bootstrap;
        extrapolation = :flat_forward
    ) where {T <: Spline.SplineCurve}
    __supports_bootstrap(mod0) || throw(
        ArgumentError(
            "Fit.Bootstrap does not support $mod0: adding a knot would change earlier " *
                "curve segments and reprice earlier quotes. Use Spline.Linear(), or " *
                __bootstrap_alternative(mod0) * "."
        )
    )
    # The solve runs on primal quotes; dual numbers in the quotes or the extrapolation policy
    # are propagated afterwards by the implicit function theorem (see `__implicit_knot_curve`).
    quotes, primal_quotes = __calibration_quotes(quotes)
    primal_extrapolation = __primal_extrapolation(extrapolation)
    order = sortperm(primal_quotes; by = maturity)
    quotes, primal_quotes = quotes[order], primal_quotes[order]
    n = length(primal_quotes)
    times = [float(maturity(q)) for q in primal_quotes]
    # A quote maturing at t = 0 does not depend on the knot rate there, so the solve would
    # return its seed. The other grid rules (a quote at all, finite and distinct maturities)
    # are `KnotGrid`'s: the returned curve's knots are exactly the quote maturities, validated
    # once, up front, with the same errors as direct construction.
    all(>(0), times) || throw(ArgumentError("bootstrap quote maturities must be positive; got $times"))
    grid0 = Yield.KnotGrid(zeros(n), times, mod0; who = "fit($(mod0), Bootstrap)")
    __check_primal_quotes(Yield.__build(mod0, grid0; extrapolation = primal_extrapolation), primal_quotes)
    zs = zeros(n)
    scales = zeros(n)

    for i in eachindex(primal_quotes)
        q = primal_quotes[i]
        # The i-th continuous zero rate is the only unknown — earlier knots are
        # already solved — so each step is a scalar root-find (exact repricing)
        # rather than an optimizer pass that rebuilds the interpolant per
        # AD-traced evaluation. The trial curve is the returned curve truncated at
        # this knot, and (for local interpolants) later knots leave it unchanged up
        # to this maturity, so the final curve reprices every quote to root precision.
        f = function (z)
            zs[i] = z
            c = __trial_curve(mod0, zs[1:i], times[1:i], primal_extrapolation)
            return present_value(c, q.instrument) - q.price
        end
        seed = i == 1 ? 0.0 : zs[i - 1] # seed with the previous zero rate
        # Solve the residual relative to the quote's size (its price, or its value at the seed
        # for a zero-price quote), so the root's precision does not depend on its notional.
        s = max(abs(q.price), abs(f(seed) + q.price))
        scales[i] = iszero(s) ? one(s) : s
        g(z) = f(z) / scales[i]
        # a generous continuous-zero-rate range for the bracketed fallback
        zs[i] = __solve_primal_root(g, seed, (-1.0, 1.0))
    end
    curve = Yield.ZeroRateCurve(zs, times, mod0; extrapolation = primal_extrapolation)
    # Every quote's cashflows must end at its maturity for later knots to leave
    # it priced. Check the returned curve so a contract that pays beyond its
    # maturity cannot drift silently.
    for (q, scale) in zip(primal_quotes, scales)
        residual = present_value(curve, q.instrument) - q.price
        abs(residual) <= 1.0e-8 * scale || throw(
            ArgumentError(
                "bootstrap could not reprice the quote maturing at $(maturity(q)) " *
                    "(residual $residual); its cashflows may extend beyond its maturity"
            )
        )
    end
    return __implicit_knot_curve(curve, quotes, primal_quotes, extrapolation)
end

function fit(mod0::Yield.SmithWilson, quotes)
    cm, ts = cashflows_timepoints(quotes)
    prices = [q.price for q in quotes]

    return Yield.SmithWilson(ts, cm, prices; ufr = mod0.ufr, α = mod0.α)

end
