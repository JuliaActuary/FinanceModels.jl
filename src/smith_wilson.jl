export SmithWilsonYield

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