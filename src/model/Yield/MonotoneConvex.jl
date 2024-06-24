"""
An unfit Monotone Convex Yield Curve Model will simply have 


"""
struct MonotoneConvex{T,U} <: AbstractYieldModel
    f::Vector{T}
    fᵈ::Vector{T}
    rates::Vector{T}
    times::Vector{U}
    # inner constructor ensures f consistency with rates at construction
    function MonotoneConvex(rates::Vector{T}, times::Vector{U}) where {T,U}
        f, fᵈ = __monotone_convex_fs(rates, times)
        new{T,U}(f, fᵈ, rates, times)
    end
end


struct MonotoneConvexUnInit
end

struct MonotoneConvexUnoptimized{T,U}
    rates::Vector{T}
    times::Vector{U}
end

MonotoneConvex() = MonotoneConvexUnInit()
function (m::MonotoneConvexUnInit)(times)
    rates = zeros(length(times))
    MonotoneConvex(rates, times)
end



function __issector1(g0, g1)
    a = (g0 > 0) && (g1 >= -2 * g0) && (-0.5 * g0 >= g1)
    b = (g0 < 0) && (g1 <= -2 * g0) && (-0.5 * g0 <= g1)
    a || b
end

function __issector2(g0, g1)
    a = (g0 > 0) && (g1 < -2 * g0)
    b = (g0 < 0) && (g1 > -2 * g0)
    a || b
end


# Hagan West - WILMOTT magazine pgs 75-81
function g(f⁻, f, fᵈ)
    @show f⁻, f, fᵈ
    @show g0 = f⁻ - fᵈ
    @show g1 = f - fᵈ
    if sign(g0) == sign(g1)
        # sector (iv)
        η = g1 / (g1 + g0)
        α = -g0 * g1 / (g1 + g0)

        if x < η
            return x -> α + (g0 - α) * ((η - x) / η)^2
        else
            return x -> α + (g1 - α) * ((x - η) / (1 - η))^2
        end


    elseif __issector1(g0, g1)
        # sector (i)
        x -> g0 * (1 - 4 * x + 3 * x^2) + g1 * (-2 * x + 3 * x^2)
    elseif __issector2(g0, g1)
        # sector (ii)
        η = (g1 + 2 * g0) / (g1 - g0)
        if x < η
            return x -> g0
        else
            return x -> g0 + (g1 - g0) * ((x - η) / (1 - η))^2
        end
    else
        # sector (iii)
        η = 3 * g1 / (g1 - g0)
        if x > η
            return x -> g1
        else
            return x -> g1 + (g0 - g1) * ((η - x) / η)^2
        end

    end
end

function g_rate(x, f⁻, f, fᵈ)
    @show x, f⁻, f, fᵈ
    @show g0 = f⁻ - fᵈ
    @show g1 = f - fᵈ
    if sign(g0) == sign(g1)
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
        @show "(i)"
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

function forward(t, rates, times)

    t, i_time, rates, times = __monotone_convex_init(t, rates, times)
    f, fᵈ = __monotone_convex_fs(rates, times)

    x = (t - times[i_time]) / (times[i_time+1] - times[i_time])
    return fᵈ[i_time+1] + g(x, f[i_time], f[i_time+1], fᵈ[i_time+1])
end

"""
    returns the index associated with the time t, an initial rate vector, and a time vector
"""
function __i_time(t, times)
    i_time = findfirst(x -> x > t, times)
    if i_time == nothing
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
        fᵈ[i] = (times[i] * rates[i] - times[i-1] * rates[i-1]) / (times[i] - times[i-1])
    end
    # step 2
    f = similar(rates, length(rates) + 1)
    # fill in middle elements first, then do 1st and last
    for i in 1:length(rates)-1
        t_prior = if i == 1
            0
        else
            times[i-1]
        end

        weight1 = (times[i] - t_prior) / (times[i+1] - t_prior)
        weight2 = (times[i+1] - times[i]) / (times[i+1] - t_prior)
        f[i+1] = weight1 * fᵈ[i+1] + weight2 * fᵈ[i]
    end
    # step 3
    # collar(a,b,c) = clamp(b, a, c)
    f[1] = fᵈ[1] - 0.5 * (f[2] - fᵈ[1])

    f[end] = fᵈ[end] - 0.5 * (f[end-1] - fᵈ[end])
    f[1] = clamp(f[1], 0, 2 * fᵈ[2])
    f[end] = clamp(f[end], 0, 2 * fᵈ[end])

    for j in 2:(length(times)-1)
        f[j] = clamp(f[j], 0, 2 * min(fᵈ[j], fᵈ[j+1]))
    end

    return f, fᵈ
end


function Base.zero(mc::MonotoneConvex, t)
    lt = last(mc.times)
    f, fᵈ = mc.f, mc.fᵈ
    if t > lt
        r = Base.zero(mc, lt)
        i_time = __i_time(t, mc.times)
        return r * lt / t + forward(lt, mc.rates, mc.times) * (1 - lt / t)
    end
    @show i_time = __i_time(t, mc.times)
    # if the time is greater than the last input time then extrapolate using the forwards

    x = if i_time == 1
        x = t / times[i_time]
    else
        x = (t - times[i_time-1]) / (times[i_time] - times[i_time-1])
    end
    G = g(f[i_time], f[i_time+1], fᵈ[i_time])
    return Continuous(1 / t * (times[i_time] * rates[i_time] + (t - times[i_time]) * fᵈ[i_time] + (times[i_time] - times[i_time-1]) * G))

    # STATUS:
    # intermediate G/other results OK
    # need to get the right rate. Instead of last formula, trying the approach on pg 39 of 
    # http://uu.diva-portal.org/smash/get/diva2:1477828/FULLTEXT01.pdf 
    # and have added QuadGK and converted the G function to return a function of x instead of a calculated value 
    # (also need to change the signature of associated test cases)
    # QuadGK.quadgk(G,t)

end

function FinanceCore.discount(mc::MonotoneConvex, t)
    r = zero(mc, t)
    return discount(r, t)
end

function Base.zero(mc::MonotoneConvexUnoptimized, t)
    lt = last(times)
    # if the time is greater than the last input time then extrapolate using the forwards
    if t > lt
        r = myzero(lt, rates, times, f, fᵈ)
        return r * lt / t + forward(lt, rates, times) * (1 - lt / t)
    end

    t, i_time, rates, times = __monotone_convex_init(t, rates, times)
    f, fᵈ = __monotone_convex_fs(mc.rates, mc.times)
    x = (t - times[i_time]) / (times[i_time+1] - times[i_time])
    G = g_rate(x, f[i_time], f[i_time+1], fᵈ[i_time+1])
    return 1 / t * (times[i_time] * rates[i_time] + (t - times[i_time]) * fᵈ[i_time+1] + (times[i_time+1] - times[i_time]) * G)

end

function FinanceCore.discount(mc::MonotoneConvexUnoptimized, t)
    r = zero(mc, t)
    return exp(-r * t)
end

times = 1:5
rates = [0.03, 0.04, 0.047, 0.06, 0.06]
# forward(5.19, rates, times)
# myzero(1, rates, times)


# using Test
# @test forward(0.5, rates, times) ≈ 0.02875
# @test forward(1, rates, times) ≈ 0.04
# @test forward(2, rates, times) ≈ 0.0555
# @test forward(2.5, rates, times) ≈ 0.0571254591368226
# @test forward(5, rates, times) ≈ 0.05025
# @test forward(5.2, rates, times) ≈ 0.05025

# @test myzero(0.5, rates, times) ≈ 0.02625
# @test myzero(1, rates, times) ≈ 0.03
# @test myzero(2, rates, times) ≈ 0.04
# @test myzero(2.5, rates, times) ≈ 0.0431375956535047
# @test myzero(5, rates, times) ≈ 0.06
# @test myzero(5.2, rates, times) ≈ 0.059625


