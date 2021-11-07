var documenterSearchIndex = {"docs":
[{"location":"api/#Yields-API-Reference","page":"API Reference","title":"Yields API Reference","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Please open an issue if you encounter any issues or confusion with the package.","category":"page"},{"location":"api/","page":"API Reference","title":"API Reference","text":"","category":"page"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [Yields]","category":"page"},{"location":"api/#Yields.AbstractYield","page":"API Reference","title":"Yields.AbstractYield","text":"An AbstractYield is an object which can be used as an argument to:\n\nzero-coupon spot rates via zero\ndiscount factor via discount\naccumulation factor via accumulation\n\nIt can be be constructed via:\n\nzero rate curve with Zero\nforward rate curve with Forward\npar rate curve with Par\ntypical OIS curve with OIS\ntypical constant maturity treasury (CMT) curve with CMT\n\n\n\n\n\n","category":"type"},{"location":"api/#Yields.BulletBondQuote","page":"API Reference","title":"Yields.BulletBondQuote","text":"BulletBondQuote(yield, price, maturity, frequency)\n\nQuote for a set of fixed interest bullet bonds with given yield, price, maturity and a given payment frequency frequency.\n\nConstruct a vector of quotes for use with SmithWilson methods, e.g. by broadcasting over an array of inputs.\n\nExamples\n\njulia> maturities = [1.2, 2.5, 3.6]\njulia> interests = [-0.02, 0.3, 0.04]\njulia> prices = [1.3, 0.1, 4.5]\njulia> frequencies = [2,1,2]\njulia> bbq = Yields.BulletBondQuote.(interests, maturities, prices, frequencies)\n\n\n\n\n\n","category":"type"},{"location":"api/#Yields.Constant","page":"API Reference","title":"Yields.Constant","text":"Constant(rate)\n\nConstruct a yield object where the spot rate is constant for all maturities. If rate is not a Rate type, will assume Periodic(1) for the compounding frequency\n\nExamples\n\njulia> y = Yields.Constant(0.05)\njulia> discount(y,2)\n0.9070294784580498     # 1 / (1.05) ^ 2\n\n\n\n\n\n","category":"type"},{"location":"api/#Yields.Continuous","page":"API Reference","title":"Yields.Continuous","text":"Continuous()\n\nA type representing continuous interest compounding frequency.\n\nExamples\n\njulia> Rate(0.01,Continuous())\nRate(0.01, Continuous())\n\nSee also: Periodic\n\n\n\n\n\n","category":"type"},{"location":"api/#Yields.Continuous-Tuple{Any}","page":"API Reference","title":"Yields.Continuous","text":"julia> Continuous(0.01)\nRate(0.01, Continuous())\n\nSee also: Periodic\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.Periodic","page":"API Reference","title":"Yields.Periodic","text":"Periodic(frequency)\n\nA type representing periodic interest compounding with the given frequency\n\nExamples\n\nCreating a semi-annual bond equivalent yield:\n\njulia> Rate(0.01,Periodic(2))\nRate(0.01, Periodic(2))\n\nSee also: Continuous\n\n\n\n\n\n","category":"type"},{"location":"api/#Yields.Periodic-Tuple{Any,Any}","page":"API Reference","title":"Yields.Periodic","text":"Periodic(rate,frequency)\n\nA convenience constructor for Rate(rate,Periodic(frequency)).\n\nExamples\n\nCreating a semi-annual bond equivalent yield:\n\njulia> Periodic(0.01,2)\nRate(0.01, Periodic(2))\n\nSee also: Continuous\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.Rate-Tuple{Any}","page":"API Reference","title":"Yields.Rate","text":"Rate(rate[,frequency=1])\nRate(rate,frequency::CompoundingFrequency)\n\nRate is a type that encapsulates an interest rate along with its compounding frequency.\n\nPeriodic rates can be constructed via Rate(rate,frequency) or Rate(rate,Periodic(frequency)).\n\nContinuous rates can be constructed via Rate(rate, Inf) or Rate(rate,Continuous()).\n\nExamples\n\njulia> Rate(0.01,Continuous())\nRate(0.01, Continuous())\n\njulia> Rate(0.01,Periodic(2))\nRate(0.01, Periodic(2))\n\njulia> Rate(0.01)\nRate(0.01, Periodic(1))\n\njulia> Rate(0.01,2)\nRate(0.01, Periodic(2))\n\njulia> Rate(0.01,Periodic(4))\nRate(0.01, Periodic(4))\n\njulia> Rate(0.01,Inf)\nRate(0.01, Continuous())\n\njulia> Rate(0.01,Continuous())\nRate(0.01, Continuous())\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.SmithWilson","page":"API Reference","title":"Yields.SmithWilson","text":"SmithWilson(zcq::Vector{ZeroCouponQuote}; ufr, α)\nSmithWilson(swq::Vector{SwapQuote}; ufr, α)\nSmithWilson(bbq::Vector{BulletBondQuote}; ufr, α)\nSmithWilson(times<:AbstractVector, cashflows<:AbstractMatrix, prices<:AbstractVector; ufr, α)\nSmithWilson(u, qb; ufr, α)\n\nCreate a yield curve object that implements the Smith-Wilson interpolation/extrapolation scheme.\n\nPositional arguments to construct a curve:\n\nQuoted instrument as the first argument: either a Vector of ZeroCouponQuotes, SwapQuotes, or BulletBondQuotes, or \nA set of times, cashflows, and prices, or\nA curve can be with u is the timepoints coming from the calibration, and qb is the internal parameterization of the curve that ensures that the calibration is correct. Users may prefer the other constructors but this mathematical constructor is also available.\n\nRequired keyword arguments:\n\nufr is the Ultimate Forward Rate, the forward interest rate to which the yield curve tends, in continuous compounding convention. \nα is the parameter that governs the speed of convergence towards the Ultimate Forward Rate. It can be typed with \\alpha[TAB]\n\n\n\n\n\n","category":"type"},{"location":"api/#Yields.Step","page":"API Reference","title":"Yields.Step","text":"Step(rates,times)\n\nCreate a yield curve object where the applicable rate is the effective rate of interest applicable until corresponding time. If rates is not a Vector{Rate}, will assume Periodic(1) type.\n\nExamples\n\njulia>y = Yields.Step([0.02,0.05], [1,2])\n\njulia>rate(y,0.5)\n0.02\n\njulia>rate(y,1.5)\n0.05\n\njulia>rate(y,2.5)\n0.05\n\n\n\n\n\n","category":"type"},{"location":"api/#Yields.SwapQuote","page":"API Reference","title":"Yields.SwapQuote","text":"SwapQuote(yield, maturity, frequency)\n\nQuote for a set of interest rate swaps with the given yield and maturity and a given payment frequency.\n\nExamples\n\njulia> maturities = [1.2, 2.5, 3.6]\njulia> interests = [-0.02, 0.3, 0.04]\njulia> prices = [1.3, 0.1, 4.5]\njulia> frequencies = [2,1,2]\njulia> swq = Yields.SwapQuote.(interests, maturities, frequencies)\n\n\n\n\n\n","category":"type"},{"location":"api/#Yields.ZeroCouponQuote","page":"API Reference","title":"Yields.ZeroCouponQuote","text":"ZeroCouponQuote(price, maturity)\n\nQuote for a set of zero coupon bonds with given price and maturity. \n\nExamples\n\njulia> prices = [1.3, 0.1, 4.5]\njulia> maturities = [1.2, 2.5, 3.6]\njulia> swq = Yields.ZeroCouponQuote.(prices, maturities)\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.:+-Tuple{Yields.AbstractYield,Yields.AbstractYield}","page":"API Reference","title":"Base.:+","text":"Yields.AbstractYield + Yields.AbstractYield\n\nThe addition of two yields will create a RateCombination. For rate, discount, and accumulation purposes the spot rates of the two curves will be added together.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.:--Tuple{Yields.AbstractYield,Yields.AbstractYield}","page":"API Reference","title":"Base.:-","text":"Yields.AbstractYield - Yields.AbstractYield\n\nThe subtraction of two yields will create a RateCombination. For rate, discount, and accumulation purposes the spot rates of the second curves will be subtracted from the first.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.convert-Tuple{Yields.CompoundingFrequency,Rate{var\"#s18\",var\"#s19\"} where var\"#s19\"<:Yields.CompoundingFrequency where var\"#s18\"<:Real}","page":"API Reference","title":"Base.convert","text":"convert(T::CompoundingFrequency,r::Rate)\n\nReturns a Rate with an equivalent discount but represented with a different compounding frequency.\n\nExamples\n\njulia> r = Rate(Periodic(12),0.01)\nRate(0.01, Periodic(12))\n\njulia> convert(Periodic(1),r)\nRate(0.010045960887181016, Periodic(1))\n\njulia> convert(Continuous(),r)\nRate(0.009995835646701251, Continuous())\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.zero-Tuple{Yields.YieldCurve,Any}","page":"API Reference","title":"Base.zero","text":"zero(curve,time)\nzero(curve,time,CompoundingFrequency)\n\nReturn the zero rate for the curve at the given time. If not specified, will use Periodic(1) compounding.\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.CMT-Union{Tuple{T}, Tuple{Array{T,1},Any}} where T<:Real","page":"API Reference","title":"Yields.CMT","text":"Takes CMT yields (bond equivalent), and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual.\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.Forward-Tuple{Any,Any}","page":"API Reference","title":"Yields.Forward","text":"Forward(rate_vector,maturities)\n\nTakes a vector of 1-period forward rates and constructs a discount curve.\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.H-Union{Tuple{T}, Tuple{Any,T,T}} where T","page":"API Reference","title":"Yields.H","text":"H(α, t1, t2)\n\nThe Smith-Wilson H function implemented in a faster way.\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.H_ordered-Tuple{Any,Any,Any}","page":"API Reference","title":"Yields.H_ordered","text":"H_ordered(α, t_min, t_max)\n\nThe Smith-Wilson H function with ordered arguments (for better performance than using min and max).\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.OIS-Union{Tuple{T}, Tuple{Array{T,1},Any}} where T<:Real","page":"API Reference","title":"Yields.OIS","text":"OIS(rates,maturities)\n\nTakes Overnight Index Swap rates, and assumes that instruments <= one year maturity are settled once and other agreements are settled quarterly with a corresponding CompoundingFrequency\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.Par-Tuple{Array{var\"#s13\",1} where var\"#s13\"<:Rate,Any}","page":"API Reference","title":"Yields.Par","text":"Par(rate, maturity)\n\nConstruct a curve given a set of bond equivalent yields and the corresponding maturities. Assumes that maturities <= 1 year do not pay coupons and that after one year, pays coupons with frequency equal to the CompoundingFrequency of the corresponding rate.\n\nExamples\n\n\njulia> par = [6.,8.,9.5,10.5,11.0,11.25,11.38,11.44,11.48,11.5] ./ 100\njulia> maturities = [t for t in 1:10]\njulia> curve = Par(par,maturities);\njulia> zero(curve,1)\nRate(0.06000000000000005, Periodic(1))\n\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.Zero-Tuple{Any,Any}","page":"API Reference","title":"Yields.Zero","text":"Zero(rates,maturities)\n\nConstruct a yield curve with given zero-coupon spot rates at the given maturities. If rates is not a Vector{Rate}, will assume Periodic(1) type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.accumulation-Union{Tuple{T}, Tuple{T,Any}} where T<:Yields.AbstractYield","page":"API Reference","title":"Yields.accumulation","text":"accumulation(rate,from,to)\n\nThe accumulation factor for the rate for times from through to. If rate is a Real number, will assume a Constant interest rate.\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.cashflows-Tuple{Any,Any,Any}","page":"API Reference","title":"Yields.cashflows","text":"cashflows(interests, maturities, frequency)\ntimepoints(zcq::Vector{ZeroCouponQuote})\ntimepoints(bbq::Vector{BulletBondQuote})\n\nProduce a cash flow matrix for a set of instruments with given interests and maturities and a given payment frequency frequency. All instruments are assumed to have their first payment at time 1/frequency and have their last payment at the largest multiple of 1/frequency less than or equal to the input maturity.\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.discount-Tuple{Any,Any}","page":"API Reference","title":"Yields.discount","text":"discount(rate,to)\ndiscount(rate,from,to)\n\nThe discount factor for the rate for times from through to. If rate is a Real number, will assume a Constant interest rate.\n\n\n\n\n\n","category":"method"},{"location":"api/#Yields.timepoints-Union{Tuple{Array{Q,1}}, Tuple{Q}} where Q<:Yields.ObservableQuote","page":"API Reference","title":"Yields.timepoints","text":"timepoints(zcq::Vector{ZeroCouponQuote})\ntimepoints(bbq::Vector{BulletBondQuote})\n\nReturn the times associated with the cashflows of the instruments.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Yields","category":"page"},{"location":"#Yields.jl","page":"Home","title":"Yields.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Stable) (Image: Dev) (Image: Build Status) (Image: Coverage) (Image: lifecycle)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Yields provides a simple interface for constructing, manipulating, and using yield curves for modeling purposes.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It's intended to provide common functionality around modeling interest rates, spreads, and miscellaneous yields across the JuliaActuary ecosystem (though not limited to use in JuliaActuary packages).","category":"page"},{"location":"#QuickStart","page":"Home","title":"QuickStart","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using Yields\n\nriskfree_maturities = [0.5, 1.0, 1.5, 2.0]\nriskfree    = [5.0, 5.8, 6.4, 6.8] ./ 100     #spot rates, annual effective if unspecified\n\nspread_maturities = [0.5, 1.0, 1.5, 3.0]      # different maturities\nspread    = [1.0, 1.8, 1.4, 1.8] ./ 100       # spot spreads\n\nrf_curve = Yields.Zero(riskfree,riskfree_maturities)\nspread_curve = Yields.Zero(spread,spread_maturities)\n\n\nyield = rf_curve + spread_curve               # additive combination of the two curves\n\ndiscount(yield,1.5) # 1 / (1 + 0.064 + 0.014) ^ 1.5","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"#Rates","page":"Home","title":"Rates","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Rates are types that wrap scalar values to provide information about how to determine discount and accumulation factors.","category":"page"},{"location":"","page":"Home","title":"Home","text":"There are two CompoundingFrequency types:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Periodic(m) for rates that compound m times per period (e.g. m times per year if working with annual rates).\nContinuous() for continuously compounding rates.","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Rate(0.05,Continuous())       # 5% continuously compounded\nContinuous(0.05)              # alternate constructor\n\nRate(0.05, Periodic(2))       # 5% compounded twice per period\nPeriodic(0.05, 2)             # alternate constructor\n\n# construct a vector of rates with the given compounding\nRate.(0.02,0.03,0.04,Periodic(2)) ","category":"page"},{"location":"#Yields","page":"Home","title":"Yields","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are a several ways to construct a yield curve object.","category":"page"},{"location":"#Bootstrapping-Methods","page":"Home","title":"Bootstrapping Methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"rates can be a vector of Rates described above, or will assume Yields.Periodic(1) if the functions are given Real number values","category":"page"},{"location":"","page":"Home","title":"Home","text":"Yields.Zero(rates,maturities)  using a vector of zero, or spot, rates\nYields.Forward(rates,maturities) using a vector of one-period (or periods-long) forward rates\nYields.Constant(rate) takes a single constant rate for all times\nYields.Step(rates,maturities) doesn't interpolate - the rate is flat up to the corresponding time in times\nYields.Par(rates,maturities) takes a series of yields for securities priced at par.Assumes that maturities <= 1 year do not pay coupons and that after one year, pays coupons with frequency equal to the CompoundingFrequency of the corresponding rate.\nYields.CMT(rates,maturities) takes the most commonly presented rate data (e.g. Treasury.gov) and bootstraps the curve given the combination of bills and bonds.\nYields.OIS(rates,maturities) takes the most commonly presented rate data for overnight swaps and bootstraps the curve.","category":"page"},{"location":"#Kernel-Methods","page":"Home","title":"Kernel Methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Yields.SmithWilson curve (used for discounting in the EU Solvency II framework) can be constructed either directly by specifying its inner representation or by calibrating to a set of cashflows with known prices.\nThese cashflows can conveniently be constructed with a Vector of Yields.ZeroCouponQuotes, Yields.SwapQuotes, or Yields.BulletBondQuotes.","category":"page"},{"location":"#Functions","page":"Home","title":"Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Most of the above yields have the following defined (goal is to have them all):","category":"page"},{"location":"","page":"Home","title":"Home","text":"discount(curve,from,to) or discount(curve,to) gives the discount factor\naccumulation(curve,from,to) or accumulation(curve,to) gives the accumulation factor\nforward(curve,from,to) gives the average rate between the two given times\nzero(curve,time) or zero(curve,time,CompoundingFrequency) gives the zero-coupon spot rate for the given time.","category":"page"},{"location":"#Combinations","page":"Home","title":"Combinations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Different yield objects can be combined with addition or subtraction. See the Quickstart for an example.","category":"page"},{"location":"","page":"Home","title":"Home","text":"When adding a Yields.AbstractYield with a scalar or vector, that scalar or vector will be promoted to a yield type via Yield(). For example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"y1 = Yields.Constant(0.05)\ny2 = y1 + 0.01                # y2 is a yield of 0.06","category":"page"},{"location":"#Internals","page":"Home","title":"Internals","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For time-variant yields (ie yield curves), the inputs are converted to spot rates and linearly interpolated (using Interpolations.jl).","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you want more precise curvature (e.g. cubic spline interpolation) you can pre-process your rates into a greater number of input points before creating the Yields representation. Yields.jl uses Interpolations.jl as it is a pure-Julia interpolations package and enables auto-differentiation (AD) in Yields.jl usage. For example, ActuaryUtilities.jl uses AD for duration and convexity.","category":"page"},{"location":"#Combination-Implementation","page":"Home","title":"Combination Implementation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Combinations track two different curve objects and are not combined into a single underlying data structure. This means that you may achieve better performance if you combine the rates before constructing a Yields representation. The exception to this is Constant curves, which do get combined into a single structure that is as performant as pre-combined rate structure.","category":"page"},{"location":"#Related-Packages","page":"Home","title":"Related Packages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"InterestRates.jl specializes in fast rate calculations aimed at valuing fixed income contracts, with business-day-level accuracy.\nComparative comments: Yields.jl does not try to provide as precise controls over the timing, structure, and interpolation of the curve. Instead, Yields.jl provides a minimal interface for common modeling needs.","category":"page"}]
}
