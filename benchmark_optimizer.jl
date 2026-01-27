# Benchmark script for comparing optimizer performance
# Run this before and after changes to compare load time and runtime

using Pkg
Pkg.activate(".")

println("=" ^ 60)
println("FinanceModels Optimizer Benchmark")
println("=" ^ 60)
println()

# Measure load time
println("Measuring package load time...")
load_time = @elapsed begin
    using FinanceModels
end
println("Package load time: $(round(load_time, digits=3)) seconds")
println()

# Common benchmark cases
using BenchmarkTools

println("Setting up benchmark cases...")
println()

# Case 1: Fit Yield.Constant to ZCB prices
quotes_zcb = ZCBPrice([0.9, 0.8, 0.7, 0.6])

# Case 2: Fit Yield.NelsonSiegel to a set of quotes
# Using typical market-like quotes
quotes_ns = [
    ZCBPrice(0.99, 0.5),   # 6-month
    ZCBPrice(0.97, 1.0),   # 1-year
    ZCBPrice(0.94, 2.0),   # 2-year
    ZCBPrice(0.88, 5.0),   # 5-year
    ZCBPrice(0.78, 10.0),  # 10-year
]

# Case 3: Fit Yield.NelsonSiegelSvensson to quotes
quotes_nss = [
    ZCBPrice(0.995, 0.25),  # 3-month
    ZCBPrice(0.99, 0.5),    # 6-month
    ZCBPrice(0.97, 1.0),    # 1-year
    ZCBPrice(0.94, 2.0),    # 2-year
    ZCBPrice(0.88, 5.0),    # 5-year
    ZCBPrice(0.78, 10.0),   # 10-year
    ZCBPrice(0.65, 20.0),   # 20-year
    ZCBPrice(0.50, 30.0),   # 30-year
]

println("-" ^ 60)
println("Benchmark 1: Fit Yield.Constant to 4 ZCB prices")
println("-" ^ 60)

# Warmup
fit(Yield.Constant(), quotes_zcb)

# Benchmark
b1 = @benchmark fit(Yield.Constant(), $quotes_zcb) samples=10 evals=1
display(b1)
println()

println("-" ^ 60)
println("Benchmark 2: Fit Yield.NelsonSiegel to 5 ZCB prices")
println("-" ^ 60)

# Warmup
fit(Yield.NelsonSiegel(), quotes_ns)

# Benchmark
b2 = @benchmark fit(Yield.NelsonSiegel(), $quotes_ns) samples=10 evals=1
display(b2)
println()

println("-" ^ 60)
println("Benchmark 3: Fit Yield.NelsonSiegelSvensson to 8 ZCB prices")
println("-" ^ 60)

# Warmup
fit(Yield.NelsonSiegelSvensson(), quotes_nss)

# Benchmark
b3 = @benchmark fit(Yield.NelsonSiegelSvensson(), $quotes_nss) samples=10 evals=1
display(b3)
println()

println("=" ^ 60)
println("Summary")
println("=" ^ 60)
println("Load time:           $(round(load_time, digits=3)) s")
println("Yield.Constant:      $(round(median(b1.times) / 1e6, digits=2)) ms (median)")
println("NelsonSiegel:        $(round(median(b2.times) / 1e6, digits=2)) ms (median)")
println("NelsonSiegelSvensson:$(round(median(b3.times) / 1e6, digits=2)) ms (median)")
println()
