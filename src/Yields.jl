module Yields

import BSplineKit
import ForwardDiff
using LinearAlgebra
using UnicodePlots
using Roots
import LsqFit


# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCurve
export rate, discount, accumulation, forward,
LinearSpline, QuadraticSpline, Periodic, Continuous

include("AbstractYield.jl")
include("Rate.jl")
const DEFAULT_COMPOUNDING = Yields.Continuous()

include("utils.jl")
include("bootstrap.jl")
include("SmithWilson.jl")
include("generics.jl")
include("RateCombination.jl")
include("NelsonSiegelSvensson.jl")


end
