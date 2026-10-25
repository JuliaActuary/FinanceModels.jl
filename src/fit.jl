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
    to the complete quote set at once with `fit(spline, quotes)`, which also
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
`Optim.LBFGS()`, and `Optim.Newton()` for a spline other than
`Spline.MonotoneConvex()`. A spline is fitted with its optimizer only for a `Fit.Loss`, an
`optimizer` or `solve_kwargs`: without them, `fit(spline, quotes)` solves for the curve that
reprices every quote, and runs no optimizer. To change a custom model's default:

```julia
using Optim
FinanceModels.default_optimizer(m::MyModel) = Optim.NelderMead()
```
"""
default_optimizer(m) = Optim.LBFGS()
default_optimizer(m::T) where {T <: Spline.SplineCurve} = Optim.Newton()
default_optimizer(::Spline.MonotoneConvex) = Optim.LBFGS()

__default_loss(m) = Fit.Loss(x -> x^2)

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
# `interpolation = Yield.__interpolation(spline, length(tenors))` once, outside its solve, so every
# trial curve has one concrete type.
__trial_curve(spline, interpolation, z, tenors, extrapolation) =
    Yield.__build(spline, interpolation, Yield.KnotGrid(Yield.Unchecked(), z, tenors); extrapolation)

"""
    FitConvergenceError(retcode, msg)

Thrown by [`fit`](@ref FinanceModels.fit) when the optimizer does not report a successful solve,
or when a spline fit's exact solve does not reprice its quotes. `retcode` is a `Symbol`: for an
optimizer, `:MaxIters` (it reached its iteration limit) or `:Failure` (it stopped without
converging); for an exact solve, `:MaxIters`, `:Stalled` (no step reduced the residuals) or
`:InternalLinearSolveFailed` (a singular or non-finite Jacobian). `msg` describes the failed fit.

A failed solve can leave the parameters at the starting guess, so `fit` throws rather than
return an unfitted model. Address the cause and fit again: a different starting model, a
different `optimizer`, tighter or longer `solve_kwargs` (for example
`solve_kwargs = (; iterations = 10_000)`), or quotes that the model can fit. A spline fit's exact
solve fails when no curve of its interpolation method reprices the quotes (a price that implies
a non-positive discount factor, for example); minimizing a loss with `Fit.Loss(x -> x^2)`
returns the closest curve instead. Choose the model deliberately; switching to another
interpolation method on failure changes the curve.
"""
struct FitConvergenceError{R} <: Exception
    retcode::R
    msg::String
end

function Base.showerror(io::IO, e::FitConvergenceError)
    return print(io, "FitConvergenceError: ", e.msg, " (return code ", e.retcode, ")")
end

# Minimize `loss(x, quotes)` from `x0` with Optim and return the solution. The quotes are the problem's
# data (capturing them in `loss` made a Nelson-Siegel fit about 8% slower). Every optimizer-backed fit
# either returns a successful solution or throws: a failed solve may leave the solution equal to
# `x0`, and a model built from it would make optimizer failure look like a valid (but unfitted)
# result. The solve differentiates none of the quotes, so a dual number in them throws up front;
# every evaluation also runs through `__own_duals`. `who` and `hint` name the fit in the error. The
# solve holds its owner number until it returns or throws.
function __minimize(loss, x0, quotes, optimizer, solve_kwargs; lb = nothing, ub = nothing, who, hint)
    __holds_dual(quotes) && throw(__foreign_dual_error(who, hint))
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
    gradient = AutoForwardDiff(; tag = inner)
    hessian = DifferentiationInterface.SecondOrder(AutoForwardDiff(; tag = outer), gradient)
    objective(x, qs) = convert(eltype(x), __own_duals(x -> loss(x, qs), x, eltype(x), own, who, hint))
    opt = __bounded_optimizer(optimizer, lb, ub)
    d = __objective(objective, x0, quotes, opt, gradient, hessian)
    res = __optimize(d, x0, opt, lb, ub, Optim.Options(; solve_kwargs...))
    retcode = Optim.converged(res) ? :Success : Optim.iteration_limit_reached(res) ? :MaxIters : :Failure
    retcode === :Success || throw(FitConvergenceError(retcode, "model fitting did not converge"))
    return Optim.minimizer(res)
end

# The optimizer of a solve with bounds `lb` and `ub`: an optimizer that takes no bounds runs inside
# `Fminbox`, and a `ParticleSwarm` takes them itself.
function __bounded_optimizer(opt, lb, ub)
    lb === nothing && ub === nothing && return opt
    opt isa Optim.AbstractConstrainedOptimizer && return opt
    opt isa Optim.ParticleSwarm && return Optim.ParticleSwarm(; lower = lb, upper = ub, n_particles = opt.n_particles)
    return Optim.Fminbox(opt)
end

# The objective `opt` minimizes, with the derivatives it uses: the gradient from the inner ForwardDiff
# backend, and the Hessian (for `Newton`, `NewtonTrustRegion` and `IPNewton`) from both.
function __objective(objective, x0, quotes, opt, gradient, hessian)
    p = DifferentiationInterface.Constant(quotes)
    f(x) = objective(x, quotes)
    opt isa Optim.ZerothOrderOptimizer && return Optim.NonDifferentiable(f, x0)
    gprep = DifferentiationInterface.prepare_gradient(objective, gradient, x0, p)
    g!(G, x) = DifferentiationInterface.gradient!(objective, G, gprep, gradient, x, p)
    fg!(G, x) = first(DifferentiationInterface.value_and_gradient!(objective, G, gprep, gradient, x, p))
    opt isa Union{Optim.Newton, Optim.NewtonTrustRegion, Optim.IPNewton} || return Optim.OnceDifferentiable(f, g!, fg!, x0)
    hprep = DifferentiationInterface.prepare_hessian(objective, hessian, x0, p)
    h!(H, x) = DifferentiationInterface.hessian!(objective, H, hprep, hessian, x, p)
    return Optim.TwiceDifferentiable(f, g!, fg!, h!, x0, real(zero(eltype(x0))))
end

# Run `opt`: without bounds (or a `ParticleSwarm` holding them), within the bounds (`Fminbox`, `SAMIN`),
# or under box constraints (`IPNewton`).
__optimize(d, x0, opt, lb, ub, options) = Optim.optimize(d, x0, opt, options)
__optimize(d, x0, opt::Optim.AbstractConstrainedOptimizer, lb, ub, options) = Optim.optimize(d, lb, ub, x0, opt, options)
__optimize(d, x0, opt::Optim.ConstrainedOptimizer, lb, ub, options) = Optim.optimize(
    d, Optim.TwiceDifferentiableConstraints(something(lb, eltype(x0)[]), something(ub, eltype(x0)[])), x0, opt, options
)

"""
    fit(
        model, 
        quotes, 
        method=Fit.Loss(x -> x^2);
        variables=default_variables(model),
        optimizer=default_optimizer(model),
        solve_kwargs=(;)
        )

Fit a model to a collection of quotes using a loss function and optimization method.

## Arguments
- `model`: The initial model to fit, which is generally an instantiated but un-optimized model.
- `quotes`: A collection of quotes to fit the model to.
- `method::F=Fit.Loss(x -> x^2)`: The loss function to use for fitting the model. Defaults to the squared loss function.
  - `method` can also be `Bootstrap()` with `Spline.Linear()`. A spline fit without a `method` reprices its quotes exactly (see below).
- `variables=default_variables(model)`: The variables to optimize over. This is a tuple of optic => interval pairs specifying which parameters of the model can vary. See extended help for more. Spline and knot-curve fits vary the knot rates and do not take `variables`.
- `optimizer=default_optimizer(model)`: The optimization algorithm, an Optim.jl optimizer. The default for a given model is `Optim.LBFGS()`, a quasi-Newton method with automatic differentiation via ForwardDiff. See extended help for more on customizing the solver.
- `solve_kwargs=(;)`: Keyword arguments for `Optim.Options`, such as `(; iterations = 10_000)`,
  `g_tol` or `f_abstol`. Use them to tighten a loss fit.
- `extrapolation=:flat_forward`: For `Spline.SplineCurve` fits (including
  `Spline.MonotoneConvex()`), the long-end policy beyond the last knot (see
  [`Yield.Spline`](@ref FinanceModels.Yield.Spline)). Also accepts `:flat_zero`, `:linear`,
  `Yield.FlatForwardAt(rate)`, and (except for MonotoneConvex) `:extension`.

Fitting a `Spline.SplineCurve` places one knot at each quote maturity and returns the same
curve type as [`ZeroRateCurve`](@ref): a `Yield.MonotoneConvex` for `Spline.MonotoneConvex()`
and a `Yield.Spline` otherwise. With one knot per quote, `fit(spline, quotes)` solves for the
knot rates that reprice every quote: Newton's method on the quotes' pricing residuals, with a
line search, from a fixed starting curve near 5%. It runs no optimizer. The returned curve
reprices each quote to within `sqrt(eps())` of its price (in practice to rounding), and a solve
that cannot reprice them, because no curve with this interpolation does, throws a
[`FitConvergenceError`](@ref FinanceModels.FitConvergenceError) rather than return the closest
curve. A `Fit.Loss`, an `optimizer` or `solve_kwargs` fit the spline by minimizing the loss
(squared, if only an optimizer or its settings are given) with an Optim.jl optimizer
instead: `LBFGS()` for `Spline.MonotoneConvex()`, `Newton()` otherwise. With quotes that a
curve reprices, the loss fit converges to the same curve, to the optimizer's tolerance.
`Fit.Bootstrap()` accepts `Spline.Linear()` only. Fitting an existing knot curve
(`fit(curve, quotes)`, a [`Yield.AbstractInterpolatedZeroCurve`](@ref)) fits new knot rates at
its tenors with its interpolation method and `extrapolation`. When its knots are at the quote
maturities, it reprices the quotes as the spline fit does. Otherwise, or with a `Fit.Loss`, an
`optimizer` or `solve_kwargs`, it minimizes the loss, from the same starting rates as a spline
fit (not the curve's own); its default optimizer is the spline's.
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
- The fit must reprice its quotes. Bootstrap and `fit(spline, quotes)` do; a loss fit is differentiated only if the
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

The default solver is `LBFGS()` from Optim.jl. This is a quasi-Newton method that uses automatic differentiation (ForwardDiff) to compute gradients efficiently.
 - Any Optim.jl optimizer can be used, e.g. `using Optim` and `fit(...; optimizer = Optim.NelderMead())`. When the variables are bounded, an optimizer that takes no bounds runs inside `Optim.Fminbox`, which does not take `Newton()`; `IPNewton()` takes the bounds as constraints.
 - Solver settings go in `solve_kwargs`, the keywords of `Optim.Options`, e.g. `fit(...; solve_kwargs = (; iterations = 10_000, g_tol = 1e-12))`. Under `Fminbox`, `iterations` limits each inner solve and `outer_iterations` the number of them.
 - More documentation is available from the upstream packages:
   - [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/)
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

# A spline fit places its knots at the sorted quote maturities. With a loss, it minimizes the loss.
function fit(
        mod0::T, quotes, method::F; optimizer = default_optimizer(mod0), solve_kwargs = (;),
        extrapolation = :flat_forward
    ) where {T <: Spline.SplineCurve, F <: Fit.Loss}
    quotes, primal_quotes = __calibration_quotes(quotes)
    tenors = sort!([maturity(q) for q in primal_quotes])
    return __fit_knot_rates(mod0, tenors, quotes, primal_quotes, method, extrapolation; optimizer, solve_kwargs)
end

# Without a loss, it reprices its quotes: one knot per quote makes the repricing conditions a square
# system, which Newton's method solves (`__solve_knot_rates`). An optimizer or its settings ask for a
# loss fit, of the squared loss.
function fit(
        mod0::T, quotes; optimizer = nothing, solve_kwargs = nothing, extrapolation = :flat_forward
    ) where {T <: Spline.SplineCurve}
    (optimizer === nothing && solve_kwargs === nothing) || return fit(
        mod0, quotes, __default_loss(mod0); extrapolation,
        optimizer = something(optimizer, default_optimizer(mod0)), solve_kwargs = something(solve_kwargs, (;)),
    )
    quotes, primal_quotes = __calibration_quotes(quotes)
    tenors = sort!([maturity(q) for q in primal_quotes])
    return __fit_knot_rates(mod0, tenors, quotes, primal_quotes, __ExactFit(), extrapolation; optimizer, solve_kwargs)
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

# Without a loss, a refit whose knots are at the quote maturities is a spline fit, and reprices its
# quotes as one does. Other knot grids have no exact fit in general, so they minimize the squared loss,
# as an optimizer or its settings ask.
function fit(c::Yield.AbstractInterpolatedZeroCurve, quotes; optimizer = nothing, solve_kwargs = nothing)
    quotes, primal_quotes = __calibration_quotes(quotes)
    tenors = Yield.knot_tenors(c)
    exact = optimizer === nothing && solve_kwargs === nothing &&
        length(tenors) == length(primal_quotes) && sort!([float(maturity(q)) for q in primal_quotes]) == tenors
    exact || return fit(
        c, quotes, __default_loss(c);
        optimizer = something(optimizer, default_optimizer(c.spline)), solve_kwargs = something(solve_kwargs, (;)),
    )
    return __fit_knot_rates(c.spline, tenors, quotes, primal_quotes, __ExactFit(), c.extrapolation; optimizer, solve_kwargs)
end

# Every fit of knot rates, exact or by loss, for every interpolation method including
# `Spline.MonotoneConvex()`: one trial curve per candidate built by the same `Yield.__build` as
# `ZeroRateCurve`, and a validated `ZeroRateCurve` as the result. The callers collect the quotes
# once, with their primal copies (`__calibration_quotes`). `method` is a `Fit.Loss`, or
# `__ExactFit()` to reprice the quotes.
function __fit_knot_rates(spline, tenors, quotes, primal_quotes, method, extrapolation; optimizer, solve_kwargs)
    # The solve runs on primal quotes; dual numbers in the quotes or the extrapolation policy
    # are propagated afterwards by the implicit function theorem (see `__implicit_knot_curve`).
    primal_extrapolation = __primal_extrapolation(extrapolation)
    # Validate the policy and the knot grid once, up front (duplicate maturities, too few
    # knots for the interpolant, `:extension` with MonotoneConvex, …) with the same errors as
    # direct construction; trial curves reuse them.
    who = "fit($(spline))"
    grid0 = Yield.KnotGrid(__knot_fit_seed(spline, tenors), tenors, spline; who)
    interpolation = Yield.__interpolation(spline, length(grid0.tenors))
    rates = __solve_knots(method, spline, interpolation, grid0, primal_quotes, primal_extrapolation, who; optimizer, solve_kwargs)
    curve = Yield.ZeroRateCurve(rates, grid0.tenors, spline; extrapolation = primal_extrapolation)   # public result: validated (finite rates)
    return __implicit_knot_curve(curve, quotes, primal_quotes, extrapolation)
end

struct __ExactFit end

function __solve_knots(method::Fit.Loss, spline, interpolation, grid0, primal_quotes, primal_extrapolation, who; optimizer, solve_kwargs)
    loss(u, qs) = __quote_loss(__trial_curve(spline, interpolation, u, grid0.tenors, primal_extrapolation), method, qs)
    return __minimize(loss, grid0.rates, primal_quotes, optimizer, solve_kwargs; who, hint = __SPLINE_FIT_HINT)
end
__solve_knots(::__ExactFit, spline, interpolation, grid0, primal_quotes, primal_extrapolation, who; optimizer, solve_kwargs) =
    __solve_knot_rates(__RepricingResiduals(spline, interpolation, grid0.tenors, primal_quotes, primal_extrapolation, who), grid0.rates)

# The repricing residuals of an exact knot-curve fit, Rᵢ(z) = present_value(curve(z), qᵢ.instrument) - qᵢ.price,
# over trial curves built as in a loss fit. Like a loss fit's objective, every evaluation checks that the only
# dual numbers it meets are its own Jacobian's (`__own_duals`); the solve runs on Float64 knot rates.
struct __RepricingResiduals{S, I, T, Q, E}
    spline::S
    interpolation::I
    tenors::T
    quotes::Q
    extrapolation::E
    who::String
end
function (r::__RepricingResiduals)(z)
    own = Tuple{ForwardDiff.Tag{typeof(r), Float64}}
    curve = __trial_curve(r.spline, r.interpolation, z, r.tenors, r.extrapolation)
    return map(q -> __own_duals(c -> __quote_value(c, q) - q.price, curve, eltype(z), own, r.who, __SPLINE_FIT_HINT), r.quotes)
end

# An exact fit's iteration limit, the relative residual at which it has converged, the size of a Newton step
# (in knot rates) below which it stops, and the largest relative residual it accepts when it stops: a result
# that does not reprice every quote to `sqrt(eps)`, as `__implicit_root` requires, throws.
const __EXACT_FIT_MAXITERS = 100
const __EXACT_FIT_FTOL = 1.0e-14
const __EXACT_FIT_STEPTOL = 1.0e-14
const __EXACT_FIT_RTOL = sqrt(eps(Float64))

# Newton's method on the repricing residuals `r` from the seed `z0`. Each residual is taken relative to its
# quote's size (its price, or its value at the seed for a zero-price quote, as in the bootstrap), so that the
# line search and the stopping rules do not depend on notionals; the Newton steps do not depend on that
# scaling. A backtracking line search on the sum of squared residuals accepts only steps that decrease it. The
# iteration stops at rounding: residuals within `__EXACT_FIT_FTOL`, a step within `__EXACT_FIT_STEPTOL`, or a
# step that no longer decreases the residuals.
function __solve_knot_rates(r::__RepricingResiduals, z0)
    v = r(z0)
    scale = map(r.quotes, v) do q, vᵢ
        s = max(abs(q.price), abs(vᵢ + q.price))
        iszero(s) ? one(s) : s
    end
    z, f = copy(z0), v ./ scale
    ϕ = sum(abs2, f)
    retcode = :MaxIters
    for _ in 1:__EXACT_FIT_MAXITERS
        maximum(abs, f) <= __EXACT_FIT_FTOL && (retcode = :Success; break)
        J = ForwardDiff.jacobian(r, z) ./ scale
        F = lu(J; check = false)
        (all(isfinite, J) && issuccess(F)) || (retcode = :InternalLinearSolveFailed; break)
        δ = F \ f
        if maximum(abs, δ) <= __EXACT_FIT_STEPTOL
            z = z .- δ
            f = r(z) ./ scale
            retcode = :Success
            break
        end
        α, accepted = 1.0, false
        for _ in 1:30
            zt = z .- α .* δ
            ft = r(zt) ./ scale
            ϕt = sum(abs2, ft)
            if ϕt <= (1 - 1.0e-4 * α) * ϕ   # false for a non-finite trial
                z, f, ϕ, accepted = zt, ft, ϕt, true
                break
            end
            α /= 2
        end
        accepted || (retcode = :Stalled; break)
    end
    worst = argmax(abs.(f))
    abs(f[worst]) <= __EXACT_FIT_RTOL && return z
    stopped = if retcode === :MaxIters
        "did not reprice its quotes within $(__EXACT_FIT_MAXITERS) Newton iterations"
    elseif retcode === :InternalLinearSolveFailed
        "stopped at knot rates $z, where the Jacobian of the quotes' pricing residuals is singular or not finite"
    else
        "cannot reduce its quotes' pricing residuals further"
    end
    throw(
        FitConvergenceError(
            retcode === :Success ? :Stalled : retcode,
            "$(r.who) $stopped: the largest relative residual is $(f[worst]), on the quote maturing at " *
                "$(maturity(r.quotes[worst])). No curve with this interpolation may reprice the quotes (for " *
                "example, a price that implies a non-positive discount factor), or the quote prices may not " *
                "determine every knot rate (a quote whose cashflows do not depend on its own maturity's knot). " *
                "fit(spline, quotes, Fit.Loss(x -> x^2)) minimizes the squared residuals instead."
        )
    )
end

# FX.Forwards with a spline placeholder or a knot curve as the foreign curve: given spot, the
# domestic curve, and market quotes, each quote reduces to an equivalent quote on the
# base-currency discount curve — closed-form implied zero-coupon quotes for outright
# forwards (DF_f(t) = (K·DF_d(t) + price)/spot, see `FX.implied_zcb_quotes`) and par
# cashflow strips for basis swaps (see `FX.ParBasisSwap`) — so curve construction
# reduces to fitting the spline through those quotes; no optimizer runs over the FX
# model itself, and one quote set may mix both instrument types. The two dispatch
# methods are deliberately separate (not a `Union`) so neither is ambiguous against the
# generic optic-based `fit(mod0, quotes, ::Fit.Loss)`. The method (none: a spline placeholder's
# exact fit) and the keywords (`optimizer`, `solve_kwargs`, `extrapolation`) go to the
# foreign-curve fit unchanged.
function __fit_fx_via_implied(mod0, quotes, method...; kwargs...)
    implied = map(q -> FX.__implied_foreign_quote(mod0, q), quotes)
    return @set mod0.foreign = fit(mod0.foreign, implied, method...; kwargs...)
end
fit(mod0::FX.Forwards{P, S, D, F}, quotes; kwargs...) where {P, S, D, F <: Union{Spline.SplineCurve, Yield.AbstractInterpolatedZeroCurve}} =
    __fit_fx_via_implied(mod0, quotes; kwargs...)
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.Bootstrap; kwargs...) where {P, S, D, F <: Spline.SplineCurve} = __fit_fx_via_implied(mod0, quotes, method; kwargs...)
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.Loss; kwargs...) where {P, S, D, F <: Union{Spline.SplineCurve, Yield.AbstractInterpolatedZeroCurve}} =
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
    zs, scales = __bootstrap_knots(mod0, Yield.__interpolation(mod0, n), primal_quotes, times, primal_extrapolation, who)
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

# The bootstrap's knot rates, solved in maturity order, and each quote's residual scale. A function
# barrier on the resolved `interpolation`, so that every step's trial curve has one concrete type.
function __bootstrap_knots(mod0, interpolation, primal_quotes, times, primal_extrapolation, who)
    n = length(primal_quotes)
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
        f = __PrimalResidual(
            function (z)
                zs[i] = z
                c = __trial_curve(mod0, interpolation, zs[1:i], times[1:i], primal_extrapolation)
                return present_value(c, q.instrument) - q.price
            end, who, __SPLINE_FIT_HINT
        )
        seed = i == 1 ? 0.0 : zs[i - 1] # seed with the previous zero rate
        # Solve the residual relative to the quote's size (its price, or its value at the seed
        # for a zero-price quote), so the root's precision does not depend on its notional.
        s = max(abs(q.price), abs(f(seed) + q.price))
        scales[i] = iszero(s) ? one(s) : s
        g(z) = f(z) / scales[i]
        # a generous continuous-zero-rate range for the bracketed fallback
        zs[i] = __solve_primal_root(g, seed, (-1.0, 1.0))
    end
    return zs, scales
end

function fit(mod0::Yield.SmithWilson, quotes)
    quotes = collect(quotes)
    cm, ts = cashflows_timepoints(quotes)
    prices = [q.price for q in quotes]

    return Yield.SmithWilson(ts, cm, prices; ufr = mod0.ufr, α = mod0.α)

end
