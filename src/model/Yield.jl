module Yield
import ..AbstractModel
import ..FinanceCore
import ..Spline as Sp
import ..ReadOnlyVector
import ..DataInterpolations
import ..Bond: coupon_times, __regular_schedule, __par_coupon
import ..__implicit_root, ..__primal, ..__ad_depth
import ..ForwardDiff

using ..FinanceCore: Continuous, Periodic, discount, accumulation, forward, pv, AbstractContract

export discount, zero, forward, par, implied_quote, pv, instantaneous_forward, knot_rates, knot_tenors, reconstruct

abstract type AbstractYieldModel <: AbstractModel end

include("Yield/Interface.jl")
include("Yield/LogTail.jl")
include("Yield/Constant.jl")
include("Yield/KnotGrid.jl")
include("Yield/InterpolatedCurves.jl")
include("Yield/Extrapolation.jl")
include("Yield/Kinks.jl")
include("Yield/DataInterpolationsCurve.jl")
include("Yield/SmithWilson.jl")
include("Yield/NelsonSiegelSvensson.jl")
include("Yield/CairnsPritchard.jl")
include("Yield/MonotoneConvex.jl")
include("Yield/ZeroRateCurve.jl")
include("Yield/ImpliedQuote.jl")
include("Yield/Generic.jl")
include("Yield/Composition.jl")
include("Yield/YieldShifts.jl")

end
