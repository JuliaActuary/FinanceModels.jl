module FinanceModels

import Dates
using Reexport
@reexport using FinanceCore
using FinanceCore: present_value, discount, accumulation
using OptimizationOptimJL
using OptimizationMetaheuristics
using StaticArrays
using IntervalSets
using AccessibleOptimization
using Accessors
using LinearAlgebra
using Transducers
import BSplineKit
import UnicodePlots
using Transducers: @next, complete, __foldl__, asfoldable
import Distributions



include("utils.jl")
include("Contract.jl")
include("model/Model.jl")
include("Projection.jl")
include("Fit.jl")

export Cashflow, Quote, Forward, CommonEquity, Option

using .Bond: ZCBYield, ZCBPrice, ParSwapYield, ParYield, CMTYield, ForwardYields, OISYield
export Bond, ZCBYield, ZCBPrice, ParSwapYield, ParYield, CMTYield, ForwardYields, OISYield

export Spline

export NullModel, Yield, discount, accumulation, zero, forward

using .Yield: par
export par

export Equity, Volatility
export Projection, CashflowProjection
export pv
export Fit, fit

include("precompile.jl")
end
