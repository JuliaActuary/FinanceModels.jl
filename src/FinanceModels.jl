module FinanceModels

import Dates
using FinanceCore
using Optimization, OptimizationOptimJL
using StaticArrays
using Accessors
using Transducers
using Transducers: @next, complete, __foldl__, asfoldable

export Cashflow, Bond, Quote
export NullModel, Yield, discount
export Projection, CashflowProjection
export Fit, fit

include("utils.jl")
include("Observables.jl")
include("Model.jl")
include("Projection.jl")
include("Fit.jl")

end
