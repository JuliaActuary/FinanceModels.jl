module FinanceModels

using Reexport
import Compat
@reexport using FinanceCore
using FinanceCore: present_value, discount, accumulation
using DifferentiationInterface: AutoForwardDiff
import DifferentiationInterface
import Optimization
import OptimizationOptimJL
using IntervalSets
using AccessibleModels
using Accessors
using LinearAlgebra
using Transducers
import DataInterpolations
using Transducers: @next, complete, __foldl__, asfoldable
import SpecialFunctions
import Roots
import ForwardDiff
using Random

include("utils.jl")
include("implicit.jl")
include("Contract.jl")
include("model/Model.jl")
include("Projection.jl")
include("valuation_contexts.jl")
include("fit_implicit.jl")
include("fit.jl")

export Cashflow, Quote, Forward, CommonEquity, Option, InterestRateSwap

using .Bond: ZCBYield, ZCBPrice, ParSwapYield, ParYield, CMTYield, ForwardYield, OISYield
export Bond, ZCBYield, ZCBPrice, ParSwapYield, ParYield, CMTYield, ForwardYield, OISYield

export Spline

export NullModel, Yield, discount, accumulation, zero, forward

using .Yield: ZeroRateCurve, knot_rates, knot_tenors, reconstruct
export ZeroRateCurve, knot_rates, knot_tenors, reconstruct

using .Yield: par, implied_quote
export par, implied_quote

export Equity, Volatility
export FX
export ShortRate, AbstractStochasticModel, RatePath, simulate, pv_mc, short_rate
export Projection, CashflowProjection, Models, valuation_model
export pv
export Fit, fit, FitConvergenceError

# Public, not exported
Compat.@compat public simulation_steps, simulation_times, default_variables, default_optimizer, closed_form

include("precompile.jl")
end
