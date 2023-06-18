using Test
using FinanceCore


# older Test tests below
# eventually covert these into TestItemRunner

using FinanceModels
using Test

include("generic.jl")
include("sp.jl")

include("Yield.jl")
include("RateCombination.jl")
include("SmithWilson.jl")

include("ActuaryUtilities.jl")
include("misc.jl")
include("NelsonSiegelSvensson.jl")
#TODO EconomicScenarioGenerators.jl integration tests

