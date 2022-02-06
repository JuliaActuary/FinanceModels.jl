using BenchmarkTools
cmt = [5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 6.8, 7.0, 7.1, 7.15, 7.2, 7.3, 7.35, 7.4, 7.5, 7.6, 7.6, 7.7, 7.8] ./ 100
mats = collect(0.5:0.5:10.0)
curve = Yields.CMT(cmt, mats)
@benchmark Yields.CMT(cmt, mats)

@benchmark discount(curve, 10)