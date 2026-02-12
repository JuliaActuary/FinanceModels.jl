using Test
using FinanceCore
using Accessors

# older Test tests below
# eventually covert these into TestItemRunner

using FinanceModels
using Test
using Transducers

include("generic.jl")
include("sp.jl")

include("Equity.jl")
include("Yield.jl")
include("CompositeYield.jl")
include("SmithWilson.jl")

# include("ActuaryUtilities.jl")
include("misc.jl")
include("NelsonSiegelSvensson.jl")
include("MonotoneConvex.jl")
include("ZeroRateCurve.jl")

include("extensions.jl")
include("Stochastic.jl")
#TODO EconomicScenarioGenerators.jl integration tests
