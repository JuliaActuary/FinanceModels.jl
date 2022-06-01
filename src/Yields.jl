module Yields

import BSplineKit
import ForwardDiff
using LinearAlgebra
using UnicodePlots
using Roots


# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCurve
export rate, discount, accumulation, forward,
LinearSpline, QuadraticSpline

include("Rate.jl")
const DEFAULT_COMPOUNDING = Continuous()

include("AbstractYield.jl")
include("utils.jl")
include("bootstrap.jl")
include("SmithWilson.jl")
include("generics.jl")
include("RateCombination.jl")


end
