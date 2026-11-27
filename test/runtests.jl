using FinanceModels
using ParallelTestRunner

# Every `.jl` file in this directory is a test file. Each one runs in its own module on one of
# several worker processes, longest first once a run has recorded the durations, so a file loads
# what it uses and does not depend on another test file.
#
# Run a subset with, for example, `Pkg.test("FinanceModels"; test_args = ["least_squares", "Blend"])`
# (a name matches the files that start with it). `--jobs=N` sets the number of workers, `--verbose`
# prints when each file starts, and `--help` lists the other options.

# Loaded into the module of every test file.
const init_code = quote
    using Test
    using FinanceCore
    using FinanceModels
    using Accessors
    using Transducers
end

testsuite = find_tests(@__DIR__)
# Not test files of this suite:
# - helpers/ holds definitions that more than one test file includes.
# - downstream/ActuaryUtilities.jl, the ActuaryUtilities integration tests, runs in the
#   Downstream workflow: ActuaryUtilities depends on FinanceModels, so it cannot be a test
#   dependency of a FinanceModels version it does not support yet.
filter!(test -> !startswith(first(test), r"helpers/|downstream/"), testsuite)

# A worker that has run the largest files holds close to 3 GiB, so the default number of workers is
# one per 3 GiB of free memory, and at most one per CPU thread.
runtests(FinanceModels, ARGS; testsuite, init_code, memory_per_worker = 3 * 2^30)
