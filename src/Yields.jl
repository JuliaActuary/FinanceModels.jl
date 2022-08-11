module Yields

using Reexport
using FinanceCore
using FinanceCore: Rate, rate, discount, accumulation, Periodic, Continuous, forward
@reexport using FinanceCore: Rate, rate, discount, accumulation, Periodic, Continuous, forward 
import BSplineKit
import ForwardDiff
using LinearAlgebra
using UnicodePlots
import LsqFit
import Optim

@show FinanceCore.Continuous(0.01)
# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCurve
export LinearSpline, QuadraticSpline, 
Bootstrap,NelsonSiegel,NelsonSiegelSvensson,SmithWilson

const DEFAULT_COMPOUNDING = Yields.Continuous()

include("AbstractYieldCurve.jl")

include("utils.jl")
include("bootstrap.jl")
include("SmithWilson.jl")
include("generics.jl")
include("RateCombination.jl")
include("NelsonSiegelSvensson.jl")

include("precompiles.jl")
end
