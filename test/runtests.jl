using Test
using FinanceCore
using Accessors

# older Test tests below
# eventually covert these into TestItemRunner

using FinanceModels
using Test

include("generic.jl")
include("sp.jl")

include("Yield.jl")
include("CompositeYield.jl")
include("SmithWilson.jl")

# TODO fix AU integration
#  include("ActuaryUtilities.jl")
include("misc.jl")
include("NelsonSiegelSvensson.jl")
#TODO EconomicScenarioGenerators.jl integration tests

