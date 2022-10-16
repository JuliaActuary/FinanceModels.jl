module Yields

using Reexport
using SnoopPrecompile
using FinanceCore
using FinanceCore: Rate, rate, discount, accumulation, Periodic, Continuous, forward
@reexport using FinanceCore: Rate, rate, discount, accumulation, Periodic, Continuous, forward 
import BSplineKit
import ForwardDiff
using LinearAlgebra
using UnicodePlots
import LsqFit
import Optim

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


function InexactError_hint(io::IO, ex::InexactError)
    print(io, "\nSuggestion: try calling Periodic($(ex.val), 1) instead.")
end

function __init__()
    Base.Experimental.register_error_hint(InexactError_hint, InexactError)
    nothing
end

end
