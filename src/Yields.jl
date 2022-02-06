module Yields

import BSplineKit
import ForwardDiff
using LinearAlgebra
using UnicodePlots

# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCurve
export rate, discount, accumulation, forward, spot,
    LinearSpline, CubicSpline

include("Rate.jl")
include("AbstractYield.jl")
include("utils.jl")
include("bootstrap.jl")
include("SmithWilson.jl")
include("generics.jl")
include("RateCombination.jl")


end
