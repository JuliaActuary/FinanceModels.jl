"""
EIOPA's 2027 risk-free term structure calculations: constant-forward bootstrapping,
LLFR/UFR extrapolation, undertaking-specific VA/CSSR, and interest-rate shocks.
`calculate` builds precise analytical curves; `publication` reproduces the May 2026
demonstration workbook's staged rounding. Inputs and assumptions are supplied by
the caller; this module does not acquire market data or calculate published RCS/UFR.
See the EIOPA guide for source versions and conventions.
"""
module EIOPA
using ..FinanceModels: Yield, Spline, Fit, fit, ZeroRateCurve, ParSwapYield, ReadOnlyVector,
    Continuous, Periodic, Rate, rate, zero, forward
import ..FinanceModels: ForwardDiff, discount

include("Inputs.jl")
include("Rounding.jl")
include("Basic.jl")
include("VA.jl")
include("ShockData.jl")
include("Shocks.jl")
include("FlooredCurve.jl")
include("Calculation.jl")
end
