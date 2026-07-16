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

    A singleton type which is passed to `fit` in order to bootstrap Splines. The curves are fit such that the spline passes through the zero rates of the curve. 

    A subtype of FitMethod.

    # Examples
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
# One bounded optic per knot rate. The previous `Tuple(@optic(_.rates) .=> -1 .. 1)`
# broadcast over a `ClosedInterval` (which is not iterable) and so threw on the
# first call, breaking the generic `fit(::MonotoneConvex, …)` path entirely.
__default_optic(m::Yield.MonotoneConvex) = ntuple(i -> (@optic(_.rates[i]) => -1.0 .. 1.0), length(m.rates))

# `MonotoneConvex` caches `f`/`fᵈ` derived from `(rates, times)` and has only a
# 2-arg constructor, so ConstructionBase's default positional reconstruction
# (`MonotoneConvex(f, fᵈ, rates, times)`) has no method — breaking `@set`/
# `setproperties` and hence the Accessors-driven `fit` path. Reconstruct from the
# stored rates/times, which recomputes the cache and keeps `f`/`fᵈ` consistent.
Accessors.ConstructionBase.constructorof(::Type{<:Yield.MonotoneConvex}) =
    (f, fᵈ, rates, times) -> Yield.MonotoneConvex(rates, times)
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
__default_optic(m::Equity.BlackScholesMerton{T,U,V}) where {T,U,V<:Volatility.Constant} = ((@optic(_.σ.σ) => 0.0 .. 10.0),)
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

__default_loss(m) = Fit.Loss(x -> x^2)

# One AD backend for every `fit` loss function. `SecondOrder` serves both first-order
# optimizers (LBFGS/Fminbox use the gradient, via the inner backend) and second-order
# ones (Newton/IPNewton use the Hessian), so OptimizationBase never auto-promotes a
# first-order declaration and warns. First-order paths never instantiate the Hessian,
# so declaring it costs them nothing.
const __FIT_ADTYPE = DifferentiationInterface.SecondOrder(AutoForwardDiff(), AutoForwardDiff())

# Build an `OptimizationFunction` whose loss reprices `quotes` under the model returned
# by `build(u)` (parameters → model), summing `loss_method.fn` over the price residuals.
# `build(u)` runs once per loss evaluation, then the model is reused across quotes.
function __reprice_loss(build, loss_method, quotes)
    function loss(u, _p)
        m = build(u)
        return mapreduce(+, quotes) do q
            loss_method.fn(present_value(m, q.instrument) - q.price)
        end
    end
    return Optimization.OptimizationFunction(loss, __FIT_ADTYPE)
end

"""
    fit(
        model, 
        quotes, 
        method=Fit.Loss(x -> x^2);
        variables=__default_optic(model), 
        optimizer=__default_optim(model)
        )

Fit a model to a collection of quotes using a loss function and optimization method.

## Arguments
- `model`: The initial model to fit, which is generally an instantiated but un-optimized model.
- `quotes`: A collection of quotes to fit the model to.
- `method::F=Fit.Loss(x -> x^2)`: The loss function to use for fitting the model. Defaults to the squared loss function. 
  - `method` can also be `Bootstrap()`. If this is the case, `model` should be a spline such as `Spline.Linear()`, `Spline.Cubic()`...
- `variables=__default_optic(model)`: The variables to optimize over. This is a tuple of optic => interval pairs specifying which parameters of the model can vary. See extended help for more.
- `optimizer=__default_optim(model)`: The optimization algorithm to use. The default optimization for a given model is `LBFGS()` from Optim.jl (via OptimizationOptimJL), a quasi-Newton method with automatic differentiation via ForwardDiff. See extended help for more on customizing the solver.

The optimization routine will then attempt to modify parameters of `model` to best fit the quoted prices of the contracts underlying the `quotes` by calling `present_value(model,contract)`. The optimization will minimize the loss function specified within `Fit.Loss(...)`. 

Different types of quotes are appropriate for different kinds of models. For example, if you try to value a set of equity `Option.EuroCall`s with a `Yield.Constant`, you will get an error because the `present_value(m<:Yield.Constant,o<:Option.EuroCall)` is not defined.

## Returns
- The fitted model.

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
        optimizer = __default_optim(mod0)
    ) where
    {F <: Fit.Loss}
    # find the rate that minimizes the loss function w.r.t. the calculated price vs the quotes
    # AccessibleModels uses maximization internally, so we negate the loss to minimize it
    function neg_loss_fn(m, qs)
        loss = mapreduce(+, qs) do q
            p = present_value(m, q.instrument)
            method.fn(p - q.price)
        end
        return -loss
    end
    amodel = AccessibleModel(Base.Fix2(neg_loss_fn, quotes), mod0, variables)
    tf = AccessibleModels.transformed_func(amodel)
    # `__FIT_ADTYPE` is SecondOrder so a bounds-compatible second-order optimizer
    # (e.g. `IPNewton()`) finds a Hessian; the default `Fminbox(LBFGS())` uses only the
    # gradient (via the inner `AutoForwardDiff`) and is unaffected.
    optf = Optimization.OptimizationFunction((x, p) -> convert(eltype(x), -tf(x, p)), __FIT_ADTYPE)
    x0 = collect(AccessibleModels.transformed_vec(amodel))
    bounds = AccessibleModels.transformed_bounds(amodel)
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
    prob = Optimization.OptimizationProblem(optf, x0, AccessibleModels.rawdata(amodel); lb, ub)
    sol = Optimization.solve(prob, optimizer)
    return AccessibleModels.from_transformed(sol.u, amodel)

end

function fit(mod0::Yield.MonotoneConvexUnInit, quotes, method::F=Fit.Loss(x -> x^2);
    optimizer=OptimizationOptimJL.LBFGS()
) where {F<:Fit.Loss}
    # Extract times from quotes (sorted)
    times = sort([maturity(q.instrument) for q in quotes])

    # Create loss function for MonotoneConvex
    loss = __monotone_convex_loss_function(times, method, quotes)

    # Initial guess - use a non-uniform guess to avoid NaN in MonotoneConvex
    # (a flat initial guess causes 0/0 in the g function due to division by zero
    # in sector iv); `range` requires distinct endpoints for length 1
    n = length(times)
    x0 = n == 1 ? [0.03] : collect(range(0.01, 0.05, length=n))

    prob = Optimization.OptimizationProblem(loss, x0)
    sol = Optimization.solve(prob, optimizer)
    return Yield.MonotoneConvex(sol.u, times)
end

__monotone_convex_loss_function(times, loss_method, quotes) =
    __reprice_loss(u -> Yield.MonotoneConvex(u, times), loss_method, quotes)

function fit(mod0::T, quotes, method::F) where {T <: Spline.SplineCurve, F <: Fit.Loss}
    times = sort!(maturity.(quotes))


    optf = __spline_loss_function(mod0, times, method, quotes)
    prob = Optimization.OptimizationProblem(optf, fill(0.05, length(quotes)))
    sol = Optimization.solve(prob, __default_optim(mod0))
    return Yield.Spline(mod0, times, sol.u)

end

# `Spline.MonotoneConvex` is a tag for the native Hagan-West curve, not a
# DataInterpolations spline: the generic spline paths above build `Yield.Spline`
# (no MonotoneConvex method) and bootstrap's incremental, t=0-anchored knots do
# not match Hagan-West's non-local forward construction. Both entry points route
# to the dedicated `Yield.MonotoneConvex` fit, which reprices the whole quote set
# at once (matching `Spline.MonotoneConvex`'s documented "dispatches to
# Yield.MonotoneConvex" behavior).
function fit(::Spline.MonotoneConvex, quotes, method::Fit.Loss; optimizer = OptimizationOptimJL.LBFGS())
    return fit(Yield.MonotoneConvex(), quotes, method; optimizer)
end
fit(::Spline.MonotoneConvex, quotes, ::Fit.Bootstrap) = fit(Yield.MonotoneConvex(), quotes)

# FX.Forwards with a spline placeholder as the foreign curve: given spot, the domestic
# curve, and forward quotes, the implied base-currency discount factors are closed-form
# (DF_f(t) = (K·DF_d(t) + price)/spot — see `FX.implied_zcb_quotes`), so curve
# construction reduces to fitting the spline through the implied zero-coupon quotes; no
# optimizer runs over the FX model itself. The two dispatch methods are deliberately
# separate (not a `Union`) so neither is ambiguous against the generic optic-based
# `fit(mod0, quotes, ::Fit.Loss)`.
function __fit_fx_via_implied(mod0, quotes, method)
    implied = FX.implied_zcb_quotes(mod0, quotes)
    return @set mod0.foreign = fit(mod0.foreign, implied, method)
end
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.Bootstrap) where {P, S, D, F <: Spline.SplineCurve} = __fit_fx_via_implied(mod0, quotes, method)
fit(mod0::FX.Forwards{P, S, D, F}, quotes, method::Fit.Loss) where {P, S, D, F <: Spline.SplineCurve} = __fit_fx_via_implied(mod0, quotes, method)

function fit(mod0::T, quotes, method::Fit.Bootstrap) where {T <: Spline.SplineCurve}
    quotes = sort(collect(quotes); by = maturity)
    n = length(quotes)
    times = [float(maturity(q)) for q in quotes]
    # duplicate knots make the interpolant degenerate and would otherwise
    # surface as a cryptic root-bracketing failure deep in the solve
    allunique(times) || throw(ArgumentError("bootstrap quotes must have distinct maturities; got duplicates among $times"))
    zs = zeros(n)

    for i in eachindex(quotes)
        q = quotes[i]
        # The i-th continuous zero rate is the only unknown — earlier knots are
        # already solved — so each step is a scalar root-find (exact repricing)
        # rather than an optimizer pass that rebuilds the interpolant per
        # AD-traced evaluation. The candidate curve pins z(0) to the first zero
        # rate, exactly as the returned curve does, so that (for local
        # interpolants) the final curve reprices every quote to root precision.
        f = function (z)
            zs[i] = z
            z_at_0 = i == 1 ? z : zs[1]
            c = Yield.Spline(mod0, [zero(eltype(times)); times[1:i]], [z_at_0; zs[1:i]])
            return present_value(c, q.instrument) - q.price
        end
        seed = i == 1 ? 0.0 : zs[i - 1] # seed with the previous zero rate
        zs[i] = try
            Roots.find_zero(f, seed, Roots.Order1())
        catch e
            # only a convergence failure falls back to a bracketed solve over a
            # generous continuous-zero-rate range; genuine errors raised while
            # pricing the quote must surface, not be retried
            e isa Roots.ConvergenceFailed || rethrow()
            Roots.find_zero(f, (-1.0, 1.0), Roots.A42())
        end
    end
    return Yield.Spline(mod0, [zero(eltype(times)); times], [first(zs); zs])
end

function fit(mod0::Yield.SmithWilson, quotes)
    cm, ts = cashflows_timepoints(quotes)
    prices = [q.price for q in quotes]

    return Yield.SmithWilson(ts, cm, prices; ufr = mod0.ufr, α = mod0.α)

end

__spline_loss_function(mod0::T, times, loss_method, quotes) where {T <: Spline.SplineCurve} =
    __reprice_loss(u -> Yield.Spline(mod0, times, u), loss_method, quotes)