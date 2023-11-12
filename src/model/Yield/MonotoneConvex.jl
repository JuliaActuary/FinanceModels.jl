
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
function g(x, f⁻, f, fᵈ)
    @show x, f⁻, f, fᵈ
    @show g0 = f⁻ - fᵈ
    @show g1 = f - fᵈ
    A = -2 * g0
    B = -2 * g1
    if sign(g0) == sign(g1)
        @show "(iv)"
        # sector (iv)
        η = g1 / (g1 + g0)
        α = -g0 * g1 / (g1 + g0)

        if x < η
            return α + (g0 - α) * ((η - x) / η)^2
        else
            return α + (g1 - α) * ((x - η) / (1 - η))^2
        end


    elseif __issector1(g0, g1)
        @show "(i)"
        # sector (i)
        g0 * (1 - 4 * x + 3 * x^2) + g1 * (-2 * x + 3 * x^2)
    elseif __issector2(g0, g1)
        @show "(ii)"
        # sector (ii)
        η = (g1 + 2 * g0) / (g1 - g0)
        if x < η
            return g0
        else
            return g0 + (g1 - g0) * ((x - η) / (1 - η))^2
        end
    else
        @show "(iii)"
        # sector (iii)
        η = 3 * g1 / (g1 - g0)
        if x > η
            return g1
        else
            return g1 + (g0 - g1) * ((η - x) / η)^2
        end

    end
end

function forward(t, rates, times)
    # the array indexing in the paper and psuedo-VBA is messy
    t = min(t, last(times))
    N = length(times)
    times = collect(times)
    rates = collect(rates)
    i_time = findfirst(x -> x > t, times)
    if i_time == nothing
        i_time = N
    end
    if !iszero(first(times))
        pushfirst!(times, zero(eltype(times)))
        pushfirst!(rates, first(rates))

    end

    f, fᵈ = __monotone_convex_fs(rates, times)

    @show x = (t - times[i_time]) / (times[i_time+1] - times[i_time])
    @show fᵈ, f
    @show fᵈ[i_time+1], g(x, f[i_time], f[i_time+1], fᵈ[i_time+1])
    return fᵈ[i_time+1] + g(x, f[i_time], f[i_time+1], fᵈ[i_time+1])




end

"""
    returns a pair of vectors (f and fᵈ) used in Monotone Convext Yield Curve fitting
"""
function __monotone_convex_fs(rates, times)
    # step 1
    fᵈ = map(2:length(times)) do i
        (times[i] * rates[i] - times[i-1] * rates[i-1]) / (times[i] - times[i-1])
    end
    pushfirst!(fᵈ, 0)
    # step 2
    f = map(2:length(times)-1) do i
        (times[i] - times[i-1]) / (times[i+1] - times[i-1]) * fᵈ[i+1] +
        (times[i+1] - times[i]) / (times[i+1] - times[i-1]) * fᵈ[i]
    end
    # step 3
    # collar(a,b,c) = clamp(b, a, c)
    pushfirst!(f, fᵈ[2] - 0.5 * (f[1] - fᵈ[2]))
    fᵈ[end], f[end-1], fᵈ[end]
    push!(f, fᵈ[end] - 0.5 * (f[end] - fᵈ[end]))
    f[1] = clamp(f[1], 0, 2 * fᵈ[2])
    f[end] = clamp(f[end], 0, 2 * fᵈ[end])

    for j in 2:(length(times)-1)
        f[j] = clamp(f[j], 0, 2 * min(fᵈ[j], fᵈ[j+1]))
    end

    return f, fᵈ
end
function myzero(t, rates, times)
    #TODO


end

times = 1:5
rates = [0.03, 0.04, 0.047, 0.06, 0.06]
forward(5.19, rates, times)


using Test
@test forward(0.5, rates, times) ≈ 0.02875
@test forward(1, rates, times) ≈ 0.04
@test forward(2, rates, times) ≈ 0.0555
@test forward(2.5, rates, times) ≈ 0.0571254591368226
@test forward(5, rates, times) ≈ 0.05025
@test forward(5.2, rates, times) ≈ 0.05025

