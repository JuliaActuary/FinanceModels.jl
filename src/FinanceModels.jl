module FinanceModels

import Dates
using FinanceCore
using OptimizationOptimJL
using OptimizationMetaheuristics
using StaticArrays
using IntervalSets
using AccessibleOptimization
using Accessors
using Transducers
using Transducers: @next, complete, __foldl__, asfoldable
import Distributions

export Cashflow, Bond, Quote, Forward, CommonEquity, Option
export NullModel, Yield, discount, accumulation, zero, forward, Equity, Volatility
export Projection, CashflowProjection
export pv
export Fit, fit

include("utils.jl")
include("Contract.jl")
include("model/Model.jl")
include("Projection.jl")
include("Fit.jl")

end
