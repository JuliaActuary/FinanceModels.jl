module FinanceModelsFitExt

using FinanceModels
using FinanceModels: Fit, Spline, Yield, Equity, Volatility, cashflows_timepoints
using FinanceCore
using FinanceCore: present_value, discount, accumulation, maturity
using Accessors
using IntervalSets
using StaticArrays

import AccessibleOptimization
using AccessibleOptimization: OptArgs, OptProblemSpec, solve
import DifferentiationInterface
using DifferentiationInterface: AutoForwardDiff
import Optimization
import OptimizationMetaheuristics
using OptimizationMetaheuristics: ECA
import OptimizationOptimJL

__default_optic(m::Yield.Constant) = OptArgs(@optic(_.rate.value) => -1.0 .. 1.0)
__default_optic(m::Yield.IntermediateYieldCurve{T}) where {T <: Spline.SplineCurve} = OptArgs(@optic(_.ys[end]) => 0.0 .. 1.0)
__default_optic(m::Yield.NelsonSiegel) = OptArgs(
    [
        @optic(_.τ₁) => 0.0 .. 100.0
        @optic(_.β₀) => -10.0 .. 10.0
        @optic(_.β₁) => -10.0 .. 10.0
        @optic(_.β₂) => -10.0 .. 10.0
    ]...
)
__default_optic(m::Yield.NelsonSiegelSvensson) = OptArgs(
    [
        @optic(_.τ₁) => 0.0 .. 100.0,
        @optic(_.τ₂) => 0.0 .. 100.0,
        @optic(_.β₀) => -10.0 .. 10.0,
        @optic(_.β₁) => -10.0 .. 10.0,
        @optic(_.β₂) => -10.0 .. 10.0,
        @optic(_.β₃) => -10.0 .. 10.0,
    ]...
)
__default_optic(m::Equity.BlackScholesMerton) = __default_optic(m.σ)
__default_optic(m::Volatility.Constant) = OptArgs(@optic(_.σ) => -0.0 .. 10.0)


__default_optim(m) = ECA()
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
- `variables=__default_optic(model)`: The variables to optimize over. This is an optic specifying which parameters of the modle can vary. See extended help for more.
- `optimizer=__default_optim(model)`: The optimization algorithm to use. The default optimization for a given model is ECA from Metahueristics.jl; see extended help for more on customizing the solver including setting the seed.

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

The default solver is `ECA()` from Metahueristics.jl. This is a stochastic global optimizer which will run with a random seed by default.
 - To make the seed static, you can specify the kwarg to `fit` with a customized ECA: e.g. `fit(...;optimizer=ECA(seed=123))`
 - A number of options are available for `ECA()` or you may specify a different solver.
 - More documentation is available from the upstream packages:
   - [Metaheuristics.jl](https://jmejia8.github.io/Metaheuristics.jl/stable/)
   - [Optimization.jl](https://docs.sciml.ai/Optimization/stable/)
   - [AccessibleOptimization.jl](https://gitlab.com/aplavin/AccessibleOptimization.jl)

## Defining the variables

An arbitrarily complex model may be the object we intend to fit - how does `fit` know what free variables are able to be solved for within the given model?
`variables` is a singlular or vector optic argument. What does this mean?
- An optic (or "lens") is a way to define an accessor to a given object. Example:

```julia-repl
julia> using Accessors, AccessibleOptimization, IntervalSets

julia> obj = (a = "AA", b = "BB");

julia> lens = @optic _.a
(@optic _.a)

julia> lens(obj)
"AA"
```
An optic argument is a singular or vector of lenses with an optional range of acceptable parameters. For example, we might have a model as follows where we want 
`fit` to optimize parameters `a` and `b`:

```julia
struct MyModel <:FinanceModels.AbstractModel
     a 
     b 
end

__default_optic(m::MyModel) = OptArgs([
    @optic(_.a) => 0.0 .. 100.0,
    @optic(_.b) => -10.0 .. 10.0,
]...)
```
In this way, fit know which arbitrary parameters in a given object may be modified. Technically, we are not modifying the immutable `MyModel`, but instead efficiently creating a new instance. This is enabled by [AccessibleOptimization.jl](https://gitlab.com/aplavin/AccessibleOptimization.jl).

Note that not all optimization algorithms want a bounded interval. In that case, simply leave off the paired range. The prior example would then become:

```julia
__default_optic(m::MyModel) = OptArgs([
    @optic(_.a),
    @optic(_.b),
]...)
```
```


## Additional Examples

See the tutorials in the package documentation for FinanceModels.jl or the docstrings of FinanceModels.jl's available model types.
"""
function FinanceModels.fit(
        mod0,
        quotes,
        method::F = __default_loss(mod0);
        variables = __default_optic(mod0),
        utype = __default_utype(mod0),
        optimizer = __default_optim(mod0)
    ) where
    {F <: Fit.Loss}
    # find the rate that minimizes the loss function w.r.t. the calculated price vs the quotes
    f = __loss_single_function(mod0, method, quotes)
    # some solvers want a `Vector` instead of `SVector` for utype (override __default_utype(m))
    ops = OptProblemSpec(f, utype, mod0, variables)
    sol = solve(ops, optimizer)
    return sol.uobj

end

function FinanceModels.fit(mod0::T, quotes, method::F) where {T <: Spline.SplineCurve, F <: Fit.Loss}
    times = sort!(maturity.(quotes))


    optf = __spline_loss_function(mod0, times, __default_loss(mod0), quotes)
    prob = Optimization.OptimizationProblem(optf, fill(0.05, length(quotes)))
    sol = solve(prob, __default_optim(mod0))
    return Yield.Spline(mod0, times, sol.u)

end

function FinanceModels.fit(mod0::T, quotes, method::Fit.Bootstrap) where {T <: Spline.SplineCurve}
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

function FinanceModels.fit(mod0::Yield.SmithWilson, quotes)
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

end
