module Fit

    abstract type FitMethod end


    """
        Fit.Loss(f)
        Fit.Loss()

    Fit a model by minimizing the sum of `f` applied to its pricing residuals (each quote's value
    under the model less its price), such as `abs` or `x -> x^2`. `Fit.Loss()` uses f(r) = r².

    This is the default method for [`fit`](@ref FinanceModels.fit). `optimizer` and `solve_kwargs`
    configure its optimizer. Knot-curve fits start from FinanceModels' fixed seed rates,
    rather than a supplied curve's rates, and seek a local minimum; they need not reprice the
    quotes. Spline.MonotoneConvex uses LBFGS; other splines use Newton.

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
    Loss() = Loss(__squared)
    __squared(x) = x^2

    """
        Fit.Exact(; maxiters = 100)

    Fit a knot curve that reprices every quote: with one knot per quote, the quotes' repricing
    conditions are a square system in the knot rates, which `fit` solves by Newton's method with a
    line search, for at most `maxiters` iterations. Select it explicitly with
    `fit(spline, quotes, Fit.Exact())`, which places one knot at each quote maturity, or with
    `fit(curve, quotes, Fit.Exact())` to keep an existing curve's knots. The default remains `Fit.Loss()`.

    The fitted curve reprices each quote to within `sqrt(eps())` of its size: the gross value of its
    projected cashflows on the curve (the sum of their absolute amounts, discounted). Custom
    pricing-only contracts retain a scale based on their price and value. If the solve cannot meet that tolerance, the fit throws a
    [`FitConvergenceError`](@ref FinanceModels.FitConvergenceError). This can reflect quotes that no
    curve of its interpolation method reprices (a price that implies a non-positive discount factor,
    for example) or numerical nonconvergence; `Fit.Loss()` minimizes the squared residuals instead. A
    refit with `Fit.Exact()` needs as many quotes as knots.

    A subtype of FitMethod.
    """
    struct Exact <: FitMethod
        maxiters::Int
        # Constructed by keyword only, so that options can be added without changing existing calls.
        Exact(; maxiters::Integer = 100) = new(maxiters)
    end

    """
        Bootstrap()

    A singleton type passed to `fit` to bootstrap a spline curve one quote at a time.
    Each step solves for the zero rate at the next quote maturity to match its price.
    Built-in cashflow contracts measure residuals against gross discounted projected payments.
    Custom pricing-only contracts retain their price/value scale.

    Supports `Spline.Linear()` (equivalently `Spline.PolynomialSpline(1)` or
    `Spline.BSpline(1)`). Bootstrapping requires that adding a knot on the right
    leaves every earlier curve segment unchanged, so earlier quotes stay exactly
    priced. Linear interpolation of zero rates has this property; quadratic,
    cubic, higher-order B-spline, PCHIP, and Akima interpolation do not, because
    a later knot changes the shape of earlier segments. Fit those strategies
    to the complete quote set at once with `fit(spline, quotes, Fit.Exact())`, which also
    reprices every quote.

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

    """
        Fit.LeastSquares(; weights = nothing, maxiters = 200)

    Fit a knot curve (a `Spline.SplineCurve` or an existing knot curve such as a `ZeroRateCurve`) by
    weighted least squares: the knot rates that minimize

        Σᵢ wᵢ (present_value(curve, qᵢ.instrument) - qᵢ.price)²

    over the quotes `qᵢ`, with `wᵢ = 1` when `weights` is `nothing`. Weights are a vector with one
    finite, non-negative number per quote; a weight of zero leaves its quote out. For example, with `D`
    a bond's modified duration, weights of `1 / (price * D)^2` make each weighted residual
    approximately a yield error, so that long bonds do not dominate the fit.

    Select this method explicitly with `fit(ZeroRateCurve(rates, knots, spline), quotes, Fit.LeastSquares())`,
    for example to fit ten knots to an index of a thousand bonds. The default remains `Fit.Loss()`.
    Quote maturities may repeat. There must be at least as many quotes with a positive weight as knots, each knot moving some quote's price. The solver starts
    from FinanceModels' standard initial rates, rather than the supplied curve's rates, and seeks a
    local minimum, which need not be the best fit of all curves.

    The minimization is Levenberg–Marquardt. `maxiters` limits its iterations; final polishing may take
    additional Gauss–Newton and Newton steps. A fit is accepted only where the Hessian of half the sum
    of squares (with the residuals' curvature) is positive definite, and the Newton correction, a local
    estimate of the remaining error in the knot rates, is at most `1e-10`. A fit that reprices to
    rounding can finish directly. Otherwise the fit throws a
    [`FitConvergenceError`](@ref FinanceModels.FitConvergenceError).

    With knots at the quote maturities (`fit(spline, quotes, Fit.LeastSquares())`), a curve that
    reprices the quotes is a global minimum. Convergence is not guaranteed: Monotone Convex,
    PCHIP and Akima can stop at interpolation switches even on square systems that Exact solves.
    Where the fit converges without repricing, it returns a local minimum; Exact would throw.

    Quotes must have fixed cashflows (`Cashflow`, `Bond.Fixed` and compositions of them), which
    covers every built-in quote family. Other contracts require Fit.Loss(). ForwardDiff supports
    first-order derivatives with respect to quote prices, cashflow amounts and the weights.

    A subtype of FitMethod.

    # Examples

    ```julia
    knots = [1.0, 2.0, 5.0, 10.0, 20.0]
    curve = fit(ZeroRateCurve(fill(0.04, length(knots)), knots, Spline.Cubic()), bond_quotes, Fit.LeastSquares())
    weighted = fit(ZeroRateCurve(fill(0.04, length(knots)), knots, Spline.Cubic()), bond_quotes,
        Fit.LeastSquares(weights = w))
    ```
    """
    struct LeastSquares{W} <: FitMethod
        weights::W
        maxiters::Int
        # Constructed by keyword only, so that options can be added without changing existing calls.
        LeastSquares(; weights = nothing, maxiters::Integer = 200) = new{typeof(weights)}(weights, maxiters)
    end


end

"""
    default_variables(model)

The parameters `fit` optimizes for `model`: a tuple of optic => interval pairs (or 1-tuples of an
optic, for an unbounded parameter). For example, to optimize `a` and `b` of a custom model:

```julia
struct MyModel <: FinanceModels.AbstractModel
    a
    b
end

FinanceModels.default_variables(m::MyModel) = (
    @optic(_.a) => 0.0 .. 100.0,
    @optic(_.b) => -10.0 .. 10.0,
)
```

See [`fit`](@ref FinanceModels.fit), "Defining the variables".
"""
default_variables(m::Yield.Constant) = ((@optic(_.rate.continuous_value) => -1.0 .. 1.0),)
default_variables(m::Yield.NelsonSiegel) = (
    @optic(_.τ₁) => 0.0 .. 100.0,
    @optic(_.β₀) => -10.0 .. 10.0,
    @optic(_.β₁) => -10.0 .. 10.0,
    @optic(_.β₂) => -10.0 .. 10.0,
)
default_variables(m::Yield.NelsonSiegelSvensson) = (
    @optic(_.τ₁) => 0.0 .. 100.0,
    @optic(_.τ₂) => 0.0 .. 100.0,
    @optic(_.β₀) => -10.0 .. 10.0,
    @optic(_.β₁) => -10.0 .. 10.0,
    @optic(_.β₂) => -10.0 .. 10.0,
    @optic(_.β₃) => -10.0 .. 10.0,
)
default_variables(m::Yield.CairnsPritchard) = (
    @optic(_.c₁) => 0.001 .. 10.0,
    @optic(_.c₂) => 0.001 .. 10.0,
    @optic(_.b₀) => -1.0 .. 1.0,
    @optic(_.b₁) => -10.0 .. 10.0,
    @optic(_.b₂) => -10.0 .. 10.0,
)
default_variables(m::Yield.CairnsPritchardExtended) = (
    @optic(_.c₁) => 0.001 .. 10.0,
    @optic(_.c₂) => 0.001 .. 10.0,
    @optic(_.c₃) => 0.001 .. 10.0,
    @optic(_.b₀) => -1.0 .. 1.0,
    @optic(_.b₁) => -10.0 .. 10.0,
    @optic(_.b₂) => -10.0 .. 10.0,
    @optic(_.b₃) => -10.0 .. 10.0,
)
default_variables(m::Equity.BlackScholesMerton{T, U, V}) where {T, U, V <: Volatility.Constant} = ((@optic(_.σ.σ) => 0.0 .. 10.0),)
default_variables(m::Volatility.Constant) = ((@optic(_.σ) => 0.0 .. 10.0),)
default_variables(m::ShortRate.Vasicek) = (
    @optic(_.a) => 0.0 .. 5.0,
    @optic(_.b) => -0.1 .. 0.5,
    @optic(_.σ) => 0.0 .. 1.0,
    @optic(_.initial.continuous_value) => -0.05 .. 0.2,
)
default_variables(m::ShortRate.CoxIngersollRoss) = (
    @optic(_.a) => 0.0 .. 5.0,
    @optic(_.b) => 0.0 .. 0.5,
    @optic(_.σ) => 0.0 .. 1.0,
    @optic(_.initial.continuous_value) => 0.0 .. 0.2,
)
default_variables(m::ShortRate.HullWhite) = (
    @optic(_.a) => 0.0 .. 5.0,
    @optic(_.σ) => 0.0 .. 1.0,
)
# FX.Forwards: the free variables live on the base-currency (`foreign`) curve — spot and
# the domestic curve are calibration inputs — so compose the foreign curve's own optics
# through the `foreign` field.
default_variables(m::FX.Forwards) = map(o -> __fx_foreign_optic(o), default_variables(m.foreign))
__fx_foreign_optic(o::Base.Pair) = Accessors.opcompose(@optic(_.foreign), o.first) => o.second
__fx_foreign_optic(o::Tuple) = (Accessors.opcompose(@optic(_.foreign), only(o)),)
__fx_foreign_optic(o) = Accessors.opcompose(@optic(_.foreign), o)

"""
    default_optimizer(model)

The optimizer [`fit`](@ref FinanceModels.fit) uses for `model` when none is given:
`OptimizationOptimJL.LBFGS()`, and `OptimizationOptimJL.Newton()` for a spline other than
`Spline.MonotoneConvex()`. These optimizers are used by `Fit.Loss()`, the default fitting method.
`Fit.Exact()` runs no optimizer. To change a custom model's default:

```julia
FinanceModels.default_optimizer(m::MyModel) = OptimizationOptimJL.NelderMead()
```
"""
default_optimizer(m) = OptimizationOptimJL.LBFGS()
default_optimizer(m::T) where {T <: Spline.SplineCurve} = OptimizationOptimJL.Newton()
default_optimizer(::Spline.MonotoneConvex) = OptimizationOptimJL.LBFGS()

__default_loss(m) = Fit.Loss()

# The AD backend of one solve. `SecondOrder` serves both first-order optimizers (LBFGS/Fminbox use
# the gradient, via the inner backend) and second-order ones (Newton/IPNewton use the Hessian), so
# OptimizationBase never auto-promotes a first-order declaration and warns; first-order paths never
# build the Hessian. ForwardDiff tells dual numbers apart only by their tags' types, so solves that
# overlap (a fit inside another fit's loss, or fits on several threads) need tags of different types:
# each running solve leases an owner number `O` and puts it in its tags.
struct __FitTagKey{O, N} end

# The owner numbers of the running solves. A solve leases the lowest free one and returns it when it
# ends, so sequential fits reuse owner 1 and its compiled methods, and only overlapping solves add one.
const __FIT_OWNERS = Bool[]
const __FIT_OWNERS_LOCK = ReentrantLock()
function __lease_fit_owner()
    return lock(__FIT_OWNERS_LOCK) do
        i = findfirst(!, __FIT_OWNERS)
        i === nothing ? (push!(__FIT_OWNERS, true); length(__FIT_OWNERS)) : (__FIT_OWNERS[i] = true; i)
    end
end
__return_fit_owner(i) = lock(() -> (__FIT_OWNERS[i] = false), __FIT_OWNERS_LOCK)

# The fitting loss: `loss_method.fn` summed over the price residuals of `quotes` under `model`.
__quote_loss(model, loss_method, quotes) = mapreduce(+, quotes) do q
    loss_method.fn(__quote_value(model, q) - q.price)
end

# A quote's value under the model being fitted. A par basis-swap quote is a base-currency quote: its
# leg is valued on the FX model's foreign curve, the native route (`FX.__implied_foreign_quote` checks
# its pair), as the spline fits of an FX model do.
__quote_value(model, q) = present_value(model, q.instrument)
__quote_value(m::FX.Forwards, q::Quote{<:Any, <:FX.BasisSwapLeg}) = present_value(m.foreign, FX.__implied_foreign_quote(m, q).instrument)

# A knot-curve fit's trial curve (optimizer candidate, bootstrap step, calibration Jacobian): knot
# rates `z` over tenors validated up front, built without a copy or a finite-rate check, so a
# non-finite candidate mid-search is a bad residual rather than an exception. The caller resolves

"""
    FitConvergenceError(retcode, msg)

Thrown by [`fit`](@ref FinanceModels.fit) when the optimizer does not report a successful solve,
or when a spline fit's exact solve does not reprice its quotes. `retcode` is the solver's
`SciMLBase.ReturnCode` (for example `ReturnCode.MaxIters`, `ReturnCode.Stalled` or
`ReturnCode.Failure`) and `msg` describes the failed fit.

A failed solve can leave the parameters at the starting guess, so `fit` throws rather than
return an unfitted model. Address the cause and fit again: a different starting model, a
different `optimizer`, tighter or longer solver settings (for example `Fit.Loss()` with
`solve_kwargs = (; maxiters = 10_000)`, or `Fit.Exact(maxiters = 500)`), or quotes that the model
can fit. A spline fit's exact solve fails when no curve of its interpolation method reprices the
quotes (a price that implies a non-positive discount factor, for example); use
[`Fit.Loss`](@ref FinanceModels.Fit.Loss) to minimize a loss instead.

Choose the model deliberately; switching to another interpolation method on failure changes the
curve.
"""
struct FitConvergenceError{R} <: Exception
    retcode::R
    msg::String
end

function Base.showerror(io::IO, e::FitConvergenceError)
    return print(io, "FitConvergenceError: ", e.msg, " (return code ", e.retcode, ")")
end

# Minimize `loss(x, quotes)` from `x0` and return the solution. The quotes are the problem's data
# (capturing them in `loss` made a Nelson-Siegel fit about 8% slower). Every optimizer-backed fit
# either returns a successful solution or throws: a failed solve may leave the solution equal to
# `x0`, and a model built from it would make optimizer failure look like a valid (but unfitted)
# result. Every evaluation runs through `__own_duals`; SciMLBase gives `x` the type of a dual number
# it finds in the quotes, so that case throws at the first evaluation. `who` and `hint` name the fit
# in the error. The solve holds its owner number until it returns or throws.
function __minimize(loss, x0, quotes, optimizer, solve_kwargs; lb = nothing, ub = nothing, who, hint)
    owner = __lease_fit_owner()
    try
        return __minimize(Val(owner), loss, x0, quotes, optimizer, solve_kwargs; lb, ub, who, hint)
    finally
        __return_fit_owner(owner)
    end
end
function __minimize(::Val{O}, loss, x0, quotes, optimizer, solve_kwargs; lb, ub, who, hint) where {O}
    outer = ForwardDiff.Tag(__FitTagKey{O, :outer}(), eltype(x0))
    inner = ForwardDiff.Tag(__FitTagKey{O, :inner}(), eltype(x0))  # created second, so it nests inside `outer`
    own = Tuple{typeof(outer), typeof(inner)}
    adtype = DifferentiationInterface.SecondOrder(AutoForwardDiff(; tag = outer), AutoForwardDiff(; tag = inner))
    objective(x, qs) = convert(eltype(x), __own_duals(x -> loss(x, qs), x, eltype(x), own, who, hint))
    prob = Optimization.OptimizationProblem(Optimization.OptimizationFunction(objective, adtype), x0, quotes; lb, ub)
    sol = Optimization.solve(prob, optimizer; solve_kwargs...)
    Optimization.SciMLBase.successful_retcode(sol.retcode) || throw(
        FitConvergenceError(sol.retcode, "model fitting did not converge")
    )
    return sol.u
end

"""
    fit(model, quotes[, method]; kwargs...)

Fit a model to a collection of quotes.

## Arguments
- `model`: The initial model to fit, which is generally an instantiated but un-optimized model.
- `quotes`: A collection of quotes to fit the model to.
- `method`: how to fit:
  - [`Fit.Loss(f)`](@ref FinanceModels.Fit.Loss): minimize the sum of `f` over the pricing residuals
    with an optimizer; `Fit.Loss()` is the squared loss;
  - [`Fit.Exact()`](@ref FinanceModels.Fit.Exact): reprice every quote, for a knot curve with one knot
    per quote;
  - [`Fit.LeastSquares()`](@ref FinanceModels.Fit.LeastSquares): minimize weighted squared residuals of fixed-cashflow quotes.
  - [`Fit.Bootstrap()`](@ref FinanceModels.Fit.Bootstrap): bootstrap a supported local interpolant one quote at a time.

  Spline fits, knot-curve refits, and generic optimizer fits default to `Fit.Loss()`.
- `variables=default_variables(model)`: The variables to optimize over. This is a tuple of optic => interval pairs specifying which parameters of the model can vary. See extended help for more. Spline and knot-curve fits vary the knot rates and do not take `variables`.
- `optimizer=default_optimizer(model)`: The optimization algorithm to use. The default optimization for a given model is `LBFGS()` from Optim.jl (via OptimizationOptimJL), a quasi-Newton method with automatic differentiation via ForwardDiff. See extended help for more on customizing the solver.
- `solve_kwargs=(;)`: Keyword arguments passed to `Optimization.solve` with the optimizer, such as
  `(; maxiters = 10_000, abstol = 1e-12)` or Optim.jl's `g_tol`. Use them to tighten a loss fit. The
  `optimizer` and `solve_kwargs` configure a [`Fit.Loss`](@ref FinanceModels.Fit.Loss) fit.
- `extrapolation=:flat_forward`: For `Spline.SplineCurve` fits (including
  `Spline.MonotoneConvex()`), the long-end policy beyond the last knot (see
  [`Yield.Spline`](@ref FinanceModels.Yield.Spline)). Also accepts `:flat_zero`, `:linear`,
  `Yield.FlatForwardAt(rate)`, and (except for MonotoneConvex) `:extension`.

Fitting a `Spline.SplineCurve` places one knot at each distinct quote maturity and returns the
same curve type as [`ZeroRateCurve`](@ref). Refitting an existing knot curve preserves its knots,
interpolation method and extrapolation policy. See each [`Fit`](@ref FinanceModels.Fit.Loss)
method above for its solver, acceptance criteria and limitations.
With fewer quotes than a polynomial or B-spline needs for its order, the order is reduced to one
less than the number of knots (see
[`Spline.PolynomialSpline`](@ref FinanceModels.Spline.PolynomialSpline)): `fit(Spline.Cubic(), quotes)`
with two quotes returns a linear curve.

The optimization routine will then attempt to modify parameters of `model` to best fit the quoted prices of the contracts underlying the `quotes` by calling `present_value(model,contract)`. The optimization will minimize the loss function specified within `Fit.Loss(...)`. 

Different types of quotes are appropriate for different kinds of models. For example, if you try to value a set of equity `Option.EuroCall`s with a `Yield.Constant`, you will get an error because the `present_value(m<:Yield.Constant,o<:Option.EuroCall)` is not defined.

## Returns
- The fitted model.

## Differentiating through a fit

Spline fits (`fit(spline, quotes)` and `fit(Spline.Linear(), quotes, Fit.Bootstrap())`,
including `FX.Forwards` with a spline foreign curve) and fits of an existing knot curve with one
knot per quote are differentiable with ForwardDiff: when
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
- A refitted knot curve needs one knot per quote; with another number of knots the implicit
  solve is not square and throws a `DimensionMismatch`.
- The fit must reprice its quotes. Bootstrap and `fit(spline, quotes, Fit.Exact())` do; a loss fit is differentiated only if the
  largest absolute component of one Newton correction of its knot rates towards the exact fit
  is at most `1e-6`. That correction is a local estimate of the fit's error, not a guaranteed
  distance to the exact solution.
- A fitted `Spline.MonotoneConvex()`, `Spline.PCHIP()`, or `Spline.Akima()` curve that lies on a
  kink of its interpolation (flat quotes, for example) has no derivative there and throws.
- Other models' `fit`s differentiate none of their inputs. Every solve checks its parameters
  and its loss (or residual) at each evaluation: a caller's dual number that reaches them, from
  the quotes, a contract type the fit does not strip, or a fixed model field, throws an
  `ArgumentError` instead of differentiating the solver's iterations.

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

FinanceModels.default_variables(m::MyModel) = (
    @optic(_.a) => 0.0 .. 100.0,
    @optic(_.b) => -10.0 .. 10.0,
)
```
In this way, `fit` knows which arbitrary parameters in a given object may be modified. Technically, we are not modifying the immutable `MyModel`, but instead efficiently creating a new instance. This is enabled by [AccessibleModels.jl](https://github.com/JuliaAPlavin/AccessibleModels.jl).

Note that not all optimization algorithms want a bounded interval. In that case, simply leave off the paired range. The prior example would then become:

```julia
FinanceModels.default_variables(m::MyModel) = (
    (@optic(_.a),),
    (@optic(_.b),),
)
```


## Additional Examples

See the tutorials in the package documentation for FinanceModels.jl or the docstrings of FinanceModels.jl's available model types.
"""
function fit(
        mod0,
        quotes,
        method::F = __default_loss(mod0);
        variables = default_variables(mod0),
        optimizer = default_optimizer(mod0),
        solve_kwargs = (;)
    ) where
    {F <: Fit.Loss}
    quotes = collect(quotes)   # read by every loss evaluation
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
    loss(x, qs) = __quote_loss(build(x), method, qs)
    return build(
        __minimize(
            loss, x0, quotes, optimizer, solve_kwargs; lb, ub,
            who = "fit($(nameof(typeof(mod0))))",
            hint = "This fit does not differentiate its inputs; spline fits such as fit(Spline.Linear(), quotes) " *
                "differentiate quote prices. To differentiate a valuation with respect to a fitted curve's " *
                "knot rates, use reconstruct(curve; rates = dual_rates).",
        )
    )

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
        mod0::T, quotes, method::F; optimizer = default_optimizer(mod0), solve_kwargs = (;),
        extrapolation = :flat_forward
    ) where {T <: Spline.SplineCurve, F <: Fit.Loss}
    quotes, primal_quotes = __calibration_quotes(quotes)
    tenors = sort!([maturity(q) for q in primal_quotes])
    return __fit_knot_rates(mod0, tenors, quotes, primal_quotes, method, extrapolation; optimizer, solve_kwargs)
end

function fit(mod0::T, quotes, method::Fit.Exact; extrapolation = :flat_forward) where {T <: Spline.SplineCurve}
    quotes, primal_quotes = __calibration_quotes(quotes)
    tenors = sort!([maturity(q) for q in primal_quotes])
    return __fit_knot_rates(mod0, tenors, quotes, primal_quotes, method, extrapolation)
end

# A least-squares spline fit has one knot per quote maturity. A curve that reprices the quotes is a global
# minimum; where none does, the fit returns a local one, where `Fit.Exact()` throws.
function fit(mod0::T, quotes, method::Fit.LeastSquares; extrapolation = :flat_forward) where {T <: Spline.SplineCurve}
    quotes, primal_quotes = __calibration_quotes(quotes)
    tenors = sort!(maturity.(primal_quotes))
    return __fit_least_squares(mod0, tenors, quotes, primal_quotes, method, extrapolation)
end

# Refitting a knot curve keeps its knot tenors, interpolation method and extrapolation policy, and
# starts from the spline fit's seed rather than the curve's own rates: the price loss can be
# numerically flat at those (a high rate at a long tenor discounts to ~0), and the solver would
# accept the unfitted start.
function fit(
        c::Yield.AbstractInterpolatedZeroCurve, quotes, method::Fit.Loss;
        optimizer = default_optimizer(c.spline), solve_kwargs = (;)
    )
    quotes, primal_quotes = __calibration_quotes(quotes)
    return __fit_knot_rates(c.spline, Yield.knot_tenors(c), quotes, primal_quotes, method, c.extrapolation; optimizer, solve_kwargs)
end

function fit(c::Yield.AbstractInterpolatedZeroCurve, quotes, method::Fit.Exact)
    quotes, primal_quotes = __calibration_quotes(quotes)
    return __fit_knot_rates(c.spline, Yield.knot_tenors(c), quotes, primal_quotes, method, c.extrapolation)
end

# A least-squares refit minimizes the weighted squared pricing residuals on the curve's knots, wherever
# they are (`__fit_least_squares`).
function fit(c::Yield.AbstractInterpolatedZeroCurve, quotes, method::Fit.LeastSquares)
    quotes, primal_quotes = __calibration_quotes(quotes)
    return __fit_least_squares(c.spline, Yield.knot_tenors(c), quotes, primal_quotes, method, c.extrapolation)
end

# Every fit of knot rates, exact or by loss, for every interpolation method including
# `Spline.MonotoneConvex()`: one trial curve per candidate built by the same `Yield.__build` as
# `ZeroRateCurve`, and a validated `ZeroRateCurve` as the result. The callers collect the quotes
# once, with their primal copies (`__calibration_quotes`). `method` is a `Fit.Loss`, or `Fit.Exact`
# to reprice the quotes.
function __fit_knot_rates(spline, tenors, quotes, primal_quotes, method, extrapolation; kwargs...)
    # The solve runs on primal quotes; dual numbers in the quotes or the extrapolation policy
    # are propagated afterwards by the implicit function theorem (see `__implicit_knot_curve`).
    primal_extrapolation = __primal_extrapolation(extrapolation)
    # Validate the policy and the knot grid once, up front (duplicate maturities, too few
    # knots for the interpolant, `:extension` with MonotoneConvex, …) with the same errors as
    # direct construction; trial curves reuse them.
    who = "fit($(spline))"
    grid0 = Yield.KnotGrid(__knot_fit_seed(spline, tenors), tenors, spline; who)
    interpolation = Yield.__interpolation(spline, length(grid0.tenors))
    rates = __solve_knots(method, spline, interpolation, grid0, primal_quotes, primal_extrapolation, who; kwargs...)
    curve = Yield.ZeroRateCurve(rates, grid0.tenors, spline; extrapolation = primal_extrapolation)   # public result: validated (finite rates)
    return __implicit_knot_curve(curve, quotes, primal_quotes, extrapolation)
end

function __solve_knots(method::Fit.Loss, spline, interpolation, grid0, primal_quotes, primal_extrapolation, who; optimizer, solve_kwargs)
    loss(u, qs) = __quote_loss(__trial_curve(spline, interpolation, u, grid0.tenors, primal_extrapolation), method, qs)
    return __minimize(loss, grid0.rates, primal_quotes, optimizer, solve_kwargs; who, hint = __SPLINE_FIT_HINT)
end
function __solve_knots(method::Fit.Exact, spline, interpolation, grid0, primal_quotes, primal_extrapolation, who)
    length(primal_quotes) == length(grid0.tenors) || throw(
        ArgumentError(
            "$who: Fit.Exact() solves the repricing conditions of one quote per knot; got $(length(primal_quotes)) " *
                "quotes for $(length(grid0.tenors)) knots. Fit.Loss() minimizes the squared residuals instead."
        )
    )
    r = __RepricingResiduals(spline, interpolation, grid0.tenors, primal_quotes, primal_extrapolation, who)
    return __solve_knot_rates(r, grid0.rates, method.maxiters)
end

# FX.Forwards with a spline placeholder or a knot curve as the foreign curve: given spot, the
# domestic curve, and market quotes, each quote reduces to an equivalent quote on the
# base-currency discount curve — closed-form implied zero-coupon quotes for outright
# forwards (DF_f(t) = (K·DF_d(t) + price)/spot, see `FX.implied_zcb_quotes`) and par
# cashflow strips for basis swaps (see `FX.ParBasisSwap`) — so curve construction
# reduces to fitting the spline through those quotes; no optimizer runs over the FX
# model itself, and one quote set may mix both instrument types. The two dispatch
# methods are deliberately separate (not a `Union`) so neither is ambiguous against the
# generic optic-based `fit(mod0, quotes, ::Fit.Loss)`. The method and the keywords (`optimizer`, `solve_kwargs`, `extrapolation`) go to the
# foreign-curve fit unchanged.
function __fit_fx_via_implied(mod0, quotes, method; kwargs...)
    implied = map(q -> FX.__implied_foreign_quote(mod0, q), quotes)
    return @set mod0.foreign = fit(mod0.foreign, implied, method; kwargs...)
end
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.Bootstrap; kwargs...) where {P, S, D, F <: Spline.SplineCurve} = __fit_fx_via_implied(mod0, quotes, method; kwargs...)
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.Loss; kwargs...) where {P, S, D, F <: Union{Spline.SplineCurve, Yield.AbstractInterpolatedZeroCurve}} =
    __fit_fx_via_implied(mod0, quotes, method; kwargs...)
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.Exact; kwargs...) where {P, S, D, F <: Union{Spline.SplineCurve, Yield.AbstractInterpolatedZeroCurve}} =
    __fit_fx_via_implied(mod0, quotes, method; kwargs...)
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.LeastSquares; kwargs...) where {P, S, D, F <: Union{Spline.SplineCurve, Yield.AbstractInterpolatedZeroCurve}} =
    __fit_fx_via_implied(mod0, quotes, method; kwargs...)

# Appending a knot must preserve every earlier segment, not merely the earlier
# knot values. In particular, "local" smooth splines need not have this property.
__supports_bootstrap(::Spline.SplineCurve) = false
__supports_bootstrap(s::Union{Spline.PolynomialSpline, Spline.BSpline}) = s.order == 1

__bootstrap_alternative(s::Spline.SplineCurve) = "fit it to all quotes at once with fit($s, quotes)"
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
    who = "fit($(mod0), Bootstrap)"
    grid0 = Yield.KnotGrid(zeros(n), times, mod0; who)
    # Bootstrapped interpolants are linear (`__supports_bootstrap`), and no grid size reduces a
    # linear interpolation, so every step builds its trial curve through the full grid's.
    zs = __bootstrap_knots(mod0, Yield.__interpolation(mod0, n), primal_quotes, times, primal_extrapolation, who)
    curve = Yield.ZeroRateCurve(zs, times, mod0; extrapolation = primal_extrapolation)
    # Every quote's cashflows must end at its maturity for later knots to leave
    # it priced. Check the returned curve so a contract that pays beyond its
    # maturity cannot drift silently. Each residual is measured against its quote's size on the
    # returned curve (`__quote_size`).
    for q in primal_quotes
        residual = present_value(curve, q.instrument) - q.price
        size = __quote_size(q, curve, residual + q.price)
        __relative_residual(residual, size) <= 1.0e-8 || throw(
            ArgumentError(
                "bootstrap could not reprice the quote maturing at $(maturity(q)) " *
                    "(residual $residual); its cashflows may extend beyond its maturity"
            )
        )
    end
    return __implicit_knot_curve(curve, quotes, primal_quotes, extrapolation)
end

# Bootstrap keeps a fixed residual scale during each root search. Updating it between
# searches resolves tiny prices without the overflow/underflow of f(z)/size(z).
function __bootstrap_knots(mod0, interpolation, primal_quotes, times, primal_extrapolation, who)
    zs = zeros(length(primal_quotes))
    for i in eachindex(primal_quotes)
        q = primal_quotes[i]
        trial(z) = __trial_curve(mod0, interpolation, [view(zs, 1:(i - 1)); z], view(times, 1:i), primal_extrapolation)
        f = __PrimalResidual(z -> present_value(trial(z), q.instrument) - q.price, who, __SPLINE_FIT_HINT)
        seed = i == 1 ? 0.0 : zs[i - 1]
        z = seed
        s = __nonzero(__quote_size(q, trial(z), f(z) + q.price))
        for _ in 1:3
            g = let s = s
                z -> f(z) / s
            end
            z = __solve_primal_root(g, z, (-1.0, 1.0))
            s_root = __quote_size(q, trial(z), f(z) + q.price)
            s_root >= s / 2 && break
            s = __nonzero(s_root)
        end
        zs[i] = z
    end
    return zs
end

function fit(mod0::Yield.SmithWilson, quotes)
    quotes = collect(quotes)
    cm, ts = cashflows_timepoints(quotes)
    prices = [q.price for q in quotes]

    return Yield.SmithWilson(ts, cm, prices; ufr = mod0.ufr, α = mod0.α)

end
