# newer TestItem tests here

using TestItemRunner

@run_package_tests

# older Test tests below
# eventually covert these into TestItemRunner

using FinanceModels
using Test

include("generic.jl")
include("bootstrap.jl")
include("RateCombination.jl")
include("SmithWilson.jl")

include("ActuaryUtilities.jl")
include("misc.jl")
include("NelsonSiegelSvensson.jl")
#TODO EconomicScenarioGenerators.jl integration tests

