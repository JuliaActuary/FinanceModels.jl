using Aqua

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
