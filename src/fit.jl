module Fit

    abstract type FitMethod end


    """
        Fit.Loss(function)

    `function` should be a loss measure, such as `x->x^2` or `x->abs(x)`. This is used by the optimization algorithm in `fit` to determine optimal parameters as defined by this loss function.

    A subtype of FitMethod.

    # Examples
    ```julia
    julia> mod0 = Yield.Constant();

    julia> quotes = ZCBPrice([0.9, 0.8, 0.7,0.6]);

    julia> fit(mod0,quotes,Fit.Loss(x-x^2))

                  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Yield Curve (FinanceModels.Yield.Constant)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀           
                  ┌────────────────────────────────────────────────────────────┐           
         0.120649 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ Zero rates
                  │⠀⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  │⠀⣧⢰⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  │⠀⣿⣾⠀⣀⣸⠀⢸⢳⣇⢀⣀⣀⣀⣀⣀⠀⡀⣀⣀⣀⡀⡀⣀⢀⣀⡀⡀⣀⢀⡀⣀⡀⢀⣀⡀⢀⡀⢀⣀⡀⢀⡀⠀⣀⡀⢀⡀⢀⣀⡀⢀⣀⠀⣀⡀⢀⣀⠀⢀│           
                  │⢠⢻⡟⡆⣿⡟⣦⠚⠀⢸⣾⠛⠛⠘⠛⠘⢲⡗⠛⠃⠛⠓⠓⠛⠚⠛⠑⠓⠛⠃⠓⠛⠑⠚⡟⠓⢻⡗⠚⠀⠓⠚⠑⠒⠃⠓⠚⠑⠚⠀⠓⠃⠘⠒⠃⠓⠃⠘⠒⠃│           
                  │⢸⢸⡇⢹⡏⠁⠉⠀⠀⠈⠉⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
       Continuous │⢸⢸⡇⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  │⢸⠀⠁⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  │⢸⠀⠀⠘⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  │⡎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  │⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  │⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
         0.120649 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
                  └────────────────────────────────────────────────────────────┘           
                  ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀time⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀30⠀  

    ```

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
__default_optic(m::Yield.MonotoneConvexUnInit) = Tuple(@optic(_.rates) .=> -1.0 .. 1.0)
__default_optic(m::Yield.MonotoneConvex) = Tuple(@optic(_.rates) .=> -1.0 .. 1.0)
__default_optic(m::Yield.IntermediateYieldCurve{T}) where {T <: Spline.SplineCurve} = ((@optic(_.ys[end]) => 0.0 .. 1.0),)
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
__default_optic(m::Equity.BlackScholesMerton{T,U,V}) where {T,U,V<:Volatility.Constant} = ((@optic(_.σ.σ) => 0.0 .. 10.0),)
__default_optic(m::Volatility.Constant) = ((@optic(_.σ) => 0.0 .. 10.0),)
__default_optic(m::ShortRate.Vasicek) = (
    @optic(_.a) => 0.0 .. 5.0,
    @optic(_.b) => -0.1 .. 0.5,
    @optic(_.σ) => 0.0 .. 1.0,
)
__default_optic(m::ShortRate.CoxIngersollRoss) = (
    @optic(_.a) => 0.0 .. 5.0,
    @optic(_.b) => 0.0 .. 0.5,
    @optic(_.σ) => 0.0 .. 1.0,
)
__default_optic(m::ShortRate.HullWhite) = (
    @optic(_.a) => 0.0 .. 5.0,
    @optic(_.σ) => 0.0 .. 1.0,
)


__default_optim(m) = OptimizationOptimJL.LBFGS()
__default_optim(m::T) where {T <: Spline.SplineCurve} = OptimizationOptimJL.Newton()

__default_utype(m) = SVector

__default_loss(m) = Fit.Loss(x -> x^2)

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

Different types of quotes are appropriate for different kinds of models. For example, if you try to value a set of equtiy `EuroCall`s with a `Yield.Constant`, you will get an error because the `present_value(m<:Yield.Constant,o<:EuroCall)` is not defined.

## Returns
- The fitted model.

# Examples
```julia
julia> model = Yield.Constant();

julia> quotes = ZCBPrice([0.9, 0.8, 0.7,0.6]);

julia> fit(model,quotes)

              ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Yield Curve (FinanceModels.Yield.Constant)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀           
              ┌────────────────────────────────────────────────────────────┐           
     0.120649 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ Zero rates
              │⠀⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              │⠀⣧⢰⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              │⠀⣿⣾⠀⣀⣸⠀⢸⢳⣇⢀⣀⣀⣀⣀⣀⠀⡀⣀⣀⣀⡀⡀⣀⢀⣀⡀⡀⣀⢀⡀⣀⡀⢀⣀⡀⢀⡀⢀⣀⡀⢀⡀⠀⣀⡀⢀⡀⢀⣀⡀⢀⣀⠀⣀⡀⢀⣀⠀⢀│           
              │⢠⢻⡟⡆⣿⡟⣦⠚⠀⢸⣾⠛⠛⠘⠛⠘⢲⡗⠛⠃⠛⠓⠓⠛⠚⠛⠑⠓⠛⠃⠓⠛⠑⠚⡟⠓⢻⡗⠚⠀⠓⠚⠑⠒⠃⠓⠚⠑⠚⠀⠓⠃⠘⠒⠃⠓⠃⠘⠒⠃│           
              │⢸⢸⡇⢹⡏⠁⠉⠀⠀⠈⠉⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
   Continuous │⢸⢸⡇⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              │⢸⠀⠁⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              │⢸⠀⠀⠘⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              │⡎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              │⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              │⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
     0.120649 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│           
              └────────────────────────────────────────────────────────────┘           
              ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀time⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀30⠀  

```

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
        utype = __default_utype(mod0),
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
    optf = Optimization.OptimizationFunction((x, p) -> convert(eltype(x), -tf(x, p)), AutoForwardDiff())
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
    # (a flat initial guess causes 0/0 in the g function due to division by zero in sector iv)
    n = length(times)
    x0 = collect(range(0.01, 0.05, length=n))

    prob = Optimization.OptimizationProblem(loss, x0)
    sol = Optimization.solve(prob, optimizer)
    return Yield.MonotoneConvex(sol.u, times)
end

function __monotone_convex_loss_function(times, loss_method, quotes)
    function loss(u, p)
        m = Yield.MonotoneConvex(u, times)
        return mapreduce(+, quotes) do q
            pv = present_value(m, q.instrument)
            loss_method.fn(pv - q.price)
        end
    end
    return Optimization.OptimizationFunction(loss, AutoForwardDiff())
end

function fit(mod0::T, quotes, method::F) where {T <: Spline.SplineCurve, F <: Fit.Loss}
    times = sort!(maturity.(quotes))


    optf = __spline_loss_function(mod0, times, __default_loss(mod0), quotes)
    prob = Optimization.OptimizationProblem(optf, fill(0.05, length(quotes)))
    sol = Optimization.solve(prob, __default_optim(mod0))
    return Yield.Spline(mod0, times, sol.u)

end

function fit(mod0::T, quotes, method::Fit.Bootstrap) where {T <: Spline.SplineCurve}
    discount_vector = [0.0]
    times = [maturity(quotes[1])]

    discount_vector[1] = let
        m = fit(Yield.Constant(), [quotes[1]], Fit.Loss(x -> x^2))
        discount(m, times[1])
    end

    for i in eachindex(quotes)[2:end]
        q = quotes[i]
        push!(times, maturity(q))
        push!(discount_vector, 0.0)
        m = Yield.IntermediateYieldCurve(mod0, times, discount_vector)
        discount_vector[i] = let
            m = fit(m, [q], Fit.Loss(x -> x^2))
            discount(m, times[i])
        end

    end
    zero_vec = -log.(clamp.(discount_vector, 0.00001, 1)) ./ times
    return Yield.Spline(mod0, [zero(eltype(times)); times], [first(zero_vec); zero_vec])
    # return Yield.Spline(mod0, times, zero_vec)

end

function fit(mod0::Yield.SmithWilson, quotes)
    cm, ts = cashflows_timepoints(quotes)
    prices = [q.price for q in quotes]

    return Yield.SmithWilson(ts, cm, prices; ufr = mod0.ufr, α = mod0.α)

end

function __loss_single_function(mod0, loss_method, quotes)
    function loss(m, quotes)
        return mapreduce(+, quotes) do q
            p = present_value(m, q.instrument)
            l = loss_method.fn(p - q.price)
        end
    end
    return Base.Fix2(Optimization.OptimizationFunction(loss), quotes) # a function that takes a model and returns the loss
end

function __spline_loss_function(mod0::T, times, loss_method, quotes) where {T <: Spline.SplineCurve}
    function loss(u, p)
        m = Yield.Spline(mod0, times, u)
        return mapreduce(+, quotes) do q
            p = present_value(m, q.instrument)
            loss_method.fn(p - q.price)
        end
    end
    return Optimization.OptimizationFunction(
        loss,
        DifferentiationInterface.SecondOrder(AutoForwardDiff(), AutoForwardDiff())
    )
end