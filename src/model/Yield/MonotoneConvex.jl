"""
    MonotoneConvex(rates, times)

A Monotone Convex yield curve model implementing the Hagan-West interpolation method.

This interpolation method guarantees:
- Continuous forward rates
- Positive forward rates (when input rates imply positive forwards)
- Monotone convex forward curves that match discrete forward rates at knot points

The implementation follows the Hagan-West method as described in WILMOTT magazine.

# Examples
```julia
prices = [0.98, 0.955, 0.92, 0.88, 0.830]
times = [1, 2, 3, 4, 5]
rates = @. -log(prices) / times

c = Yield.MonotoneConvex(rates, times)
zero(c, 2.5)  # Get the zero rate at t=2.5
discount(c, 2.5)  # Get the discount factor at t=2.5
```

# References
- Hagan & West, "Interpolation Methods for Curve Construction", WILMOTT magazine
- Dehlbom, "Interpolation of the yield curve" (http://uu.diva-portal.org/smash/get/diva2:1477828/FULLTEXT01.pdf)
"""
struct MonotoneConvex{T, U} <: AbstractYieldModel
    f::Vector{T}
    fᵈ::Vector{T}
    rates::Vector{T}
    times::Vector{U}
    # inner constructor ensures f consistency with rates at construction
    function MonotoneConvex(rates::Vector{T}, times::Vector{U}) where {T, U}
        f, fᵈ = __monotone_convex_fs(rates, times)
        return new{T, U}(f, fᵈ, rates, times)
    end
end


struct MonotoneConvexUnInit
end

MonotoneConvex() = MonotoneConvexUnInit()
function (m::MonotoneConvexUnInit)(times)
    rates = zeros(length(times))
    return MonotoneConvex(rates, times)
end


function __issector1(g0, g1)
    a = (g0 > 0) && (g1 >= -2 * g0) && (-0.5 * g0 >= g1)
    b = (g0 < 0) && (g1 <= -2 * g0) && (-0.5 * g0 <= g1)
    return a || b
end

function __issector2(g0, g1)
    a = (g0 > 0) && (g1 < -2 * g0)
    b = (g0 < 0) && (g1 > -2 * g0)
    return a || b
end


# Hagan West - WILMOTT magazine pgs 75-81
# Returns the g function value at x, where g represents the deviation from discrete forward
function g(x, f⁻, f, fᵈ)
    g0 = f⁻ - fᵈ
    g1 = f - fᵈ
    if sign(g0) == sign(g1)
        # sector (iv)
        η = g1 / (g1 + g0)
        α = -g0 * g1 / (g1 + g0)

        if x < η
            return α + (g0 - α) * ((η - x) / η)^2
        else
            return α + (g1 - α) * ((x - η) / (1 - η))^2
        end

    elseif __issector1(g0, g1)
        # sector (i)
        return g0 * (1 - 4 * x + 3 * x^2) + g1 * (-2 * x + 3 * x^2)
    elseif __issector2(g0, g1)
        # sector (ii)
        η = (g1 + 2 * g0) / (g1 - g0)
        if x < η
            return g0
        else
            return g0 + (g1 - g0) * ((x - η) / (1 - η))^2
        end
    else
        # sector (iii)
        η = 3 * g1 / (g1 - g0)
        if x > η
            return g1
        else
            return g1 + (g0 - g1) * ((η - x) / η)^2
        end
    end
end

function g_rate(x, f⁻, f, fᵈ)
    g0 = f⁻ - fᵈ
    g1 = f - fᵈ
    return if sign(g0) == sign(g1)
        # sector (iv)
        η = g1 / (g1 + g0)
        α = -g0 * g1 / (g1 + g0)

        if x < η
            return α * x - (g0 - α) * ((η - x)^3 / η^2 - η) / 3
        else
            return (2 * α + g0) / 3 * η + α * (x - η) + (g1 - α) / 3 * (x - η)^3 / (1 - η)^2
        end

    elseif __issector1(g0, g1)
        # sector (i)
        g0 * (x - 2 * x^2 + x^3) + g1 * (-x^2 + x^3)
    elseif __issector2(g0, g1)
        # sector (ii)
        η = (g1 + 2 * g0) / (g1 - g0)
        if x < η
            return g0 * x
        else
            return g0 * x + ((g1 - g0) * (x - η)^3 / (1 - η)^2) / 3
        end
    else
        # sector (iii)
        η = 3 * g1 / (g1 - g0)
        if x > η
            return (2 * g1 + g0) / 3 * η + g1 * (x - η)
        else
            return g1 * x - (g0 - g1) * ((η - x)^3 / η^2 - η) / 3
        end

    end
end

# Forward rate at time t for a MonotoneConvex model
function forward(mc::MonotoneConvex, t)
    f, fᵈ, times = mc.f, mc.fᵈ, mc.times
    lt = last(times)

    # Extrapolation: constant forward beyond last time
    if t >= lt
        return fᵈ[end]
    end

    i_time = __i_time(t, times)

    if i_time == 1
        # First interval: from 0 to times[1]
        x = t / times[1]
        return fᵈ[1] + g(x, f[1], f[2], fᵈ[1])
    else
        # Interval from times[i_time-1] to times[i_time]
        t_prev = times[i_time - 1]
        t_curr = times[i_time]
        x = (t - t_prev) / (t_curr - t_prev)
        return fᵈ[i_time] + g(x, f[i_time], f[i_time + 1], fᵈ[i_time])
    end
end

"""
    returns the index associated with the time t, an initial rate vector, and a time vector
"""
function __i_time(t, times)
    i_time = findfirst(x -> x > t, times)
    if isnothing(i_time)
        i_time = lastindex(times)
    end
    return i_time
end

"""
    returns a pair of vectors (f and fᵈ) used in Monotone Convex Yield Curve fitting
"""
function __monotone_convex_fs(rates, times)
    # step 1
    fᵈ = deepcopy(rates)
    for i in 2:length(times)
        fᵈ[i] = (times[i] * rates[i] - times[i - 1] * rates[i - 1]) / (times[i] - times[i - 1])
    end
    # step 2
    f = similar(rates, length(rates) + 1)
    # fill in middle elements first, then do 1st and last
    for i in 1:(length(rates) - 1)
        t_prior = if i == 1
            0
        else
            times[i - 1]
        end

        weight1 = (times[i] - t_prior) / (times[i + 1] - t_prior)
        weight2 = (times[i + 1] - times[i]) / (times[i + 1] - t_prior)
        f[i + 1] = weight1 * fᵈ[i + 1] + weight2 * fᵈ[i]
    end
    # step 3
    # collar(a,b,c) = clamp(b, a, c)
    f[1] = fᵈ[1] - 0.5 * (f[2] - fᵈ[1])

    f[end] = fᵈ[end] - 0.5 * (f[end - 1] - fᵈ[end])
    f[1] = clamp(f[1], 0, 2 * fᵈ[2])
    f[end] = clamp(f[end], 0, 2 * fᵈ[end])

    for j in 2:(length(times) - 1)
        f[j] = clamp(f[j], 0, 2 * min(fᵈ[j], fᵈ[j + 1]))
    end

    return f, fᵈ
end


function Base.zero(mc::MonotoneConvex, t)
    f, fᵈ, rates, times = mc.f, mc.fᵈ, mc.rates, mc.times
    lt = last(times)

    # Handle t=0 case (limit is instantaneous forward at t=0)
    if t == 0
        return Continuous(f[1])
    end

    # Extrapolation beyond last time point using constant forward
    if t > lt
        r_lt = rate(Base.zero(mc, lt))  # extract scalar rate for calculation
        f_lt = fᵈ[end]  # forward at last point
        return Continuous(r_lt * lt / t + f_lt * (1 - lt / t))
    end

    i_time = __i_time(t, times)

    # Calculate normalized position x in interval and interval bounds
    if i_time == 1
        # First interval: from 0 to times[1]
        x = t / times[1]
        G = g_rate(x, f[1], f[2], fᵈ[1])
        # r(t) = (1/t) * [t * fᵈ[1] + times[1] * G(x)]
        return Continuous(fᵈ[1] + times[1] * G / t)
    else
        # Interval from times[i_time-1] to times[i_time]
        t_prev = times[i_time - 1]
        t_curr = times[i_time]
        x = (t - t_prev) / (t_curr - t_prev)
        G = g_rate(x, f[i_time], f[i_time + 1], fᵈ[i_time])
        # r(t) = (1/t) * [t_prev * r_prev + (t - t_prev) * fᵈ[i] + (t_curr - t_prev) * G(x)]
        return Continuous((t_prev * rates[i_time - 1] + (t - t_prev) * fᵈ[i_time] + (t_curr - t_prev) * G) / t)
    end
end

function FinanceCore.discount(mc::MonotoneConvex, t)
    r = zero(mc, t)
    return discount(r, t)
end
