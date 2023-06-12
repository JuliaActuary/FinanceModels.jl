module FinanceModels

import Dates
using FinanceCore
using Optimization, OptimizationOptimJL
using StaticArrays
using Accessors
using Transducers
using Transducers: @next, complete, __foldl__, asfoldable
import Distributions

export Cashflow, Bond, Quote, Forward, Equity, Option
export NullModel, Yield, discount, BlackScholesMerton
export Projection, CashflowProjection
export value
export Fit, fit

include("utils.jl")
include("Observables.jl")
include("Model.jl")
include("Projection.jl")
include("Fit.jl")

end
