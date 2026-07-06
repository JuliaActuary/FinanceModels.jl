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
include("ZeroRatePrimitive.jl")
include("TransformedYield.jl")
include("YieldShiftInvariants.jl")
include("SmithWilson.jl")

# ActuaryUtilities integration: AU < 5.9 caps FinanceCore at 2.x, so when this
# checkout is developed against an unreleased FinanceCore major the resolver can
# only satisfy the sandbox with an ancient Yields.jl-era AU — skip rather than
# exercise stale code. CI against registered versions runs these.
import ActuaryUtilities
if pkgversion(ActuaryUtilities) >= v"5"
    include("ActuaryUtilities.jl")
else
    @warn "Skipping ActuaryUtilities integration tests (resolved version $(pkgversion(ActuaryUtilities)) predates the FinanceModels-based API)"
end
include("misc.jl")
include("NelsonSiegelSvensson.jl")
include("CairnsPritchard.jl")
include("MonotoneConvex.jl")
include("ZeroRateCurve.jl")
include("instantaneous_forward.jl")
include("regressions.jl")

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
