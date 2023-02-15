module Yields

using Reexport
using SnoopPrecompile
using FinanceCore
using FinanceCore: Rate, rate, discount, accumulation, Periodic, Continuous, forward, zero
@reexport using FinanceCore: Rate, rate, discount, accumulation, Periodic, Continuous, forward, zero
import BSplineKit
import ForwardDiff # used in Bootstrap solving
using LinearAlgebra
using UnicodePlots
import LsqFit
import Optim

# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCurve
export LinearSpline, QuadraticSpline,
    Bootstrap, NelsonSiegel, NelsonSiegelSvensson, SmithWilson
export Cashflow, ZCBPrice, ZCBYield, ParYield, CMTYield, OISYield, Forward, ForwardYield,Bond, Quote,
curve

const DEFAULT_COMPOUNDING = Yields.Continuous()

include("AbstractYieldCurve.jl")

include("utils.jl")
include("Observables.jl")
include("methods.jl")
include("Constant.jl")
include("Bootstrap.jl")
# include("SmithWilson.jl")
# include("generics.jl")
# include("NelsonSiegelSvensson.jl")
include("curve.jl")
include("RateCombination.jl")

# include("precompiles.jl")


# function MethodError_hint(io::IO, ex::InexactError)
#     hint = "\nA Periodic rate requires also passing a compounding frequency." *
#     "\nFor example, call Periodic($(ex.val), 2) for a rate compounded twice per period."
#     print(io, hint)
# end

# function __init__()
#     Base.Experimental.register_error_hint(MethodError_hint, MethodError)
#     nothing
# end

end
