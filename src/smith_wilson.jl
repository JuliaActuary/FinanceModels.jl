using LinearAlgebra

export SmithWilsonYield, CalibrationInstruments
export InstrumentQuotes, ZeroCouponQuotes, SwapQuotes, BulletBondQuotes

"""
    SmithWilsonYield(ufr,alpha,u,qb)

Create a yield curve object that implements the Smith-Wilson interpolation/extrapolation scheme.
"""
struct SmithWilsonYield <: AbstractYield
    ufr    # Ultimate Forward Rate, continuous compounding
    alpha  # Speed of approach to UFR
    u      # Vector of maturities
    qb     # Q*b vector, same length as u
end

"""
    H_ordered(alpha,t_min,t_max)

The Smith-Wilson H function with ordered arguments (for better performance than using min and max).
"""
function H_ordered(alpha, t_min, t_max)
    return alpha * t_min + 0.5 * (exp(-alpha * (t_max + t_min)) - exp(-alpha * (t_max - t_min))) 
end

"""
    H(alpha,t1,t2)

The Smith-Wilson H function implemented in a faster way. Type constraints ensure that the two calls to
H_ordered return the same type.
"""
function H(alpha, t1::T, t2::T) where {T}
    return t1 < t2 ? H_ordered(alpha, t1, t2) : H_ordered(alpha, t2, t1)
end

"""
    discount(swy::SmithWilsonYield,t)

Discount factor for a Smith-Wilson yield curve.
"""
function discount(swy::SmithWilsonYield, t)
    return exp(-swy.ufr * t) * (1.0 + sum([H(swy.alpha, swy.u[midx], t) * swy.qb[midx] for midx in 1:length(swy.u)]))
end

# Utility method - this could be extended to all YieldCurves
discount(swy::SmithWilsonYield, av::AbstractVector) = [discount(swy, t) for t in av]

"""
    CalibrationInstruments(t, cf, p)

Cash flows for calibrating a yield curve, along with their prices and payment times.
"""
struct CalibrationInstruments
    t    # Column vector of maturities
    cf   # Matrix of cash flow for each [maturity, instrument]
    p    # Row vector of instrument prices
end

"""
    SmithWilsonYield(ufr, alpha, aci::CalibrationInstruments)

Calibrate a SmithWilsonYield from CalibrationInstruments
"""
function SmithWilsonYield(ufr, alpha, aci::CalibrationInstruments)
    Q = [aci.cf[tIdx, pIdx] * exp(-ufr * aci.t[tIdx]) for tIdx in 1:length(aci.t), pIdx in 1:length(aci.p)]
    Hx = [H(alpha, t1, t2) for t1 in aci.t, t2 in aci.t]
    q = transpose(sum(Q, dims=1))
    QHQ = Q' * Hx * Q
    b = QHQ \ (aci.p - q)
    Qb = Q * b
    return SmithWilsonYield(ufr, alpha, aci.t, Qb)
end

"""
Abstract type for quotes for different cash flow instruments
"""
abstract type InstrumentQuotes end

"""
    ZeroCouponQuotes(prices, maturities)

Quotes for a set of zero coupon bonds.
"""
struct ZeroCouponQuotes <: InstrumentQuotes
    prices
    maturities
end

"""
    SwapQuotes(rates, maturities, freq)

Quotes for a set of interest rate swaps with the given maturites and a given payment frequency.
"""
struct SwapQuotes <: InstrumentQuotes
    rates
    maturities
    freq
end

"""

Quotes for a set of fixed interest bullet bonds with given interests and maturites and a given payment frequency.
"""
struct BulletBondQuotes <: InstrumentQuotes
    interests
    maturities
    prices
    freq
end

# Convert ZeroCouponQuotes
function CalibrationInstruments(zcq::ZeroCouponQuotes)
    n = length(zcq.maturities)
    return CalibrationInstruments(zcq.maturities, Matrix{Float64}(I, n, n), zcq.prices)
end

# Convert SwapQuotes
function CalibrationInstruments(swq::SwapQuotes)
    n_mat = swq.freq * maximum(swq.maturities)
    n_instr = length(swq.rates)
    cf = [(mIdx <= swq.freq * swq.maturities[iIdx] ? swq.rates[iIdx] / swq.freq : 0.0) + (mIdx == swq.freq * swq.maturities[iIdx] ? 1.0 : 0.0) for mIdx in 1:n_mat, iIdx in 1:n_instr]
    return CalibrationInstruments((1:n_mat) ./ swq.freq, cf, ones(n_instr))
end

# Convert BulletBondQuotes
function CalibrationInstruments(bbq::BulletBondQuotes)
    n_mat = bbq.freq * maximum(bbq.maturities)
    n_instr = length(bbq.interests)
    cf = [(mIdx <= bbq.freq * bbq.maturities[iIdx] ? bbq.interests[iIdx] / bbq.freq : 0.0) + (mIdx == bbq.freq * bbq.maturities[iIdx] ? 1.0 : 0.0) for mIdx in 1:n_mat, iIdx in 1:n_instr]
    return CalibrationInstruments(1:n_mat ./ bbq.freq, cf, bbq.prices)
end

# Utility method for calibrating Smith-Wilson directly from quotes
SmithWilsonYield(ufr, alpha, iq::InstrumentQuotes) = SmithWilsonYield(ufr, alpha, CalibrationInstruments(iq))