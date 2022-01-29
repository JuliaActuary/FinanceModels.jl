module Yields

import Interpolations
import ForwardDiff
using LinearAlgebra
using UnicodePlots

# don't export type, as the API of Yields.Zero is nicer and 
# less polluting than Zero and less/equally verbose as ZeroYieldCurve or ZeroCurve
export rate, discount, accumulation, forward, spot

include("Rate.jl")
include("AbstractYield.jl")
include("boostrap.jl")
include("SmithWilson.jl")
include("generics.jl")
include("RateCombination.jl")
include("utils.jl")


end
