"""
    Fit

A module providing types and methods for fitting financial models to market data.

To use the `fit` function, you need to load the optimization dependencies:
- OptimizationMetaheuristics
- OptimizationOptimJL  
- AccessibleOptimization
- DifferentiationInterface
- Optimization

# Example
```julia
using FinanceModels
using OptimizationMetaheuristics  # Load this to enable fit functionality

model = Yield.Constant()
quotes = ZCBPrice([0.9, 0.8, 0.7, 0.6])
fitted_model = fit(model, quotes)
```
"""
module Fit

    abstract type FitMethod end


    """
        Fit.Loss(function)

    `function` should be a loss measure, such as `x->x^2` or `x->abs(x)`. This is used by the optimization algorithm in `fit` to determine optimal parameters as defined by this loss function.

    A subtype of FitMethod.

    # Examples
    ```julia
    julia> using FinanceModels, OptimizationMetaheuristics
    
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
    ```julia
    using FinanceModels, OptimizationMetaheuristics
    
    quotes = ZCBYield([0.01, 0.02, 0.03])
    model = fit(Spline.Linear(), quotes, Fit.Bootstrap())
    ```
    """
    struct Bootstrap <: FitMethod
        # spline method
    end


end

# Stub for fit function - actual implementation loaded via extension
"""
    fit(model, quotes, method=Fit.Loss(x -> x^2); kwargs...)

Fit a model to market quotes using an optimization method.

**Note:** This function requires loading optimization packages. Please load one of:
- `using OptimizationMetaheuristics` (recommended, default optimizer)
- `using OptimizationOptimJL` (for Newton-based optimizers)

See the docstring after loading the extension for full documentation.

# Example
```julia
using FinanceModels
using OptimizationMetaheuristics  # Required for fit to work

model = Yield.Constant()
quotes = ZCBPrice([0.9, 0.8, 0.7, 0.6])
fitted_model = fit(model, quotes)
```
"""
function fit(args...; kwargs...)
    error("""
    The `fit` function requires optimization packages to be loaded.
    
    Please add one of the following to your code:
        using OptimizationMetaheuristics  # For metaheuristic optimization (default)
        using OptimizationOptimJL          # For gradient-based optimization
    
    Example:
        using FinanceModels
        using OptimizationMetaheuristics
        
        model = Yield.Constant()
        quotes = ZCBPrice([0.9, 0.8, 0.7, 0.6])
        fitted_model = fit(model, quotes)
    """)
end
