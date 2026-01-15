# Simple benchmark using @time for comparison
using FinanceModels
using Statistics

# Get the original implementation for comparison
function par_original(curve, time; frequency = 2)
    mat_disc = FinanceModels.discount(curve, time)
    coup_times = FinanceModels.Bond.coupon_times(time, frequency)
    coupon_pv = sum(FinanceModels.discount(curve, t) for t in coup_times)
    Δt = step(coup_times)
    r = (1 - mat_disc) / coupon_pv
    cfs = [t == last(coup_times) ? 1 + r : r for t in coup_times]
    cfs = [-1; cfs]
    r = FinanceModels.FinanceCore.internal_rate_of_return(cfs, [0; coup_times])
    frequency_inner = min(1 / Δt, max(1 / Δt, frequency))
    r = convert(FinanceModels.FinanceCore.Periodic(frequency_inner), r)
    return r
end

# Warmup
constant_curve = FinanceModels.Yield.Constant(0.04)
par_original(constant_curve, 5.0, frequency=2)
FinanceModels.par(constant_curve, 5.0, frequency=2)

println("="^80)
println("Par Function Benchmark Comparison")
println("="^80)
println()

# Benchmark function
function benchmark_fn(fn, curve, time, freq, iterations=10000)
    # Warmup
    fn(curve, time, frequency=freq)
    
    # Benchmark
    GC.gc()
    times = Float64[]
    allocs = Int[]
    
    for i in 1:iterations
        stats = @timed fn(curve, time, frequency=freq)
        push!(times, stats.time * 1e9)  # Convert to nanoseconds
        push!(allocs, stats.bytes)
    end
    
    return (
        median=median(times),
        mean=mean(times),
        min=minimum(times),
        allocs=round(Int, median(allocs))
    )
end

# Test cases
test_cases = [
    ("5-year semi-annual", 5.0, 2),
    ("10-year semi-annual", 10.0, 2),
    ("30-year semi-annual", 30.0, 2),
    ("10-year quarterly", 10.0, 4)
]

results = []
for (name, time, freq) in test_cases
    println("Test Case: $name")
    println("-"^80)
    
    println("Running original implementation...")
    orig = benchmark_fn(par_original, constant_curve, time, freq, 5000)
    
    println("Running optimized implementation...")
    opt = benchmark_fn(FinanceModels.par, constant_curve, time, freq, 5000)
    
    improvement = (orig.median - opt.median) / orig.median * 100
    alloc_reduction = (orig.allocs - opt.allocs) / orig.allocs * 100
    
    println("Original:  median=$(round(orig.median, digits=1))ns, allocs=$(orig.allocs) bytes")
    println("Optimized: median=$(round(opt.median, digits=1))ns, allocs=$(opt.allocs) bytes")
    println("⚡ Speed improvement: $(round(improvement, digits=1))%")
    println("💾 Allocation reduction: $(round(alloc_reduction, digits=1))%")
    println()
    
    push!(results, (name=name, speed=improvement, alloc=alloc_reduction))
end

# Summary
println("="^80)
println("Summary")
println("="^80)
println("Test Case                          Speed Improvement    Allocation Reduction")
println("-"^80)
for r in results
    println("$(rpad(r.name, 35))$(lpad(round(r.speed, digits=1), 7))%           $(lpad(round(r.alloc, digits=1), 7))%")
end
println()

avg_improvement = mean([r.speed for r in results])
avg_alloc_reduction = mean([r.alloc for r in results])
println("Average improvement:               $(lpad(round(avg_improvement, digits=1), 7))%           $(lpad(round(avg_alloc_reduction, digits=1), 7))%")
println("="^80)
