using Test
using FinanceCore
using Accessors

# older Test tests below
# eventually covert these into TestItemRunner

using FinanceModels
using Test
using Transducers

include("generic.jl")
include("projection_models.jl")
include("sp.jl")

include("Equity.jl")
include("FX.jl")
include("Yield.jl")
include("CompositeYield.jl")
include("ZeroRatePrimitive.jl")
include("TransformedYield.jl")
include("YieldShiftInvariants.jl")
include("SmithWilson.jl")

# The ActuaryUtilities integration tests (downstream/ActuaryUtilities.jl) run in the
# Downstream workflow: ActuaryUtilities depends on FinanceModels, so it cannot be a test
# dependency of a FinanceModels version it does not support yet.
include("misc.jl")
include("NelsonSiegelSvensson.jl")
include("CairnsPritchard.jl")
include("MonotoneConvex.jl")
include("ZeroRateCurve.jl")
include("Extrapolation.jl")
include("regressions.jl")
include("bootstrap.jl")
include("implied_quote.jl")
include("implicit_fit.jl")

include("extensions.jl")
include("Stochastic.jl")
#TODO EconomicScenarioGenerators.jl integration tests

using Aqua
@testset "Aqua.jl" begin
    Aqua.test_all(
        FinanceModels;
        # The persistent_tasks probe spawns a subprocess that precompiles the
        # package; for a heavy dep tree (Optimization, DataInterpolations) it
        # flakily fails to precompile within the CI runner's limits on
        # macOS/Windows ("done.log was not created"). FinanceModels spawns no
        # background tasks, so the check is disabled rather than left flaky.
        persistent_tasks = false,
        # FinanceModels deliberately extends these FinanceCore functions/types for
        # contract valuation and projection (same-org packages); the projection
        # machinery (Transducers) lives here rather than in FinanceCore.
        piracies = (
            treat_as_own = [
                FinanceCore.present_value,
                FinanceCore.internal_rate_of_return,
                FinanceCore.AbstractContract,
            ],
        ),
    )
end
