# The repricing residuals of an exact knot-curve fit, Rᵢ(z) = present_value(curve(z), qᵢ.instrument) - qᵢ.price,
# over trial curves built as in a loss fit. Like a loss fit's objective, every evaluation checks that the only
# dual numbers it meets are its own Jacobian's (`__own_duals`); the solve runs on Float64 knot rates.
struct __RepricingResiduals{S, I, T, Q, E}
    spline::S
    interpolation::I
    tenors::T
    quotes::Q
    extrapolation::E
    who::String
end
function __repricing_values(r, z, curve)
    own = Tuple{ForwardDiff.Tag{typeof(r), Float64}}
    return map(q -> __own_duals(c -> __quote_value(c, q) - q.price, curve, eltype(z), own, r.who, __SPLINE_FIT_HINT), r.quotes)
end
(r::__RepricingResiduals)(z) = __repricing_values(r, z, __trial_curve(r.spline, r.interpolation, z, r.tenors, r.extrapolation))

# Reuse the candidate curve for its residuals and gross quote sizes.
function __repricing_point(r, z)
    curve = __trial_curve(r.spline, r.interpolation, z, r.tenors, r.extrapolation)
    v = __repricing_values(r, z, curve)
    sizes = [__quote_size(q, curve; value = v[i] + q.price) for (i, q) in enumerate(r.quotes)]
    return v, sizes
end

# An exact fit's relative residual at which it has converged, the size of a Newton step (in knot rates) below
# which it stops, and the largest relative residual it accepts when it stops: a result that does not reprice
# every quote to `sqrt(eps)`, as `__implicit_root` requires, throws. Its iteration limit is `Fit.Exact`'s.
const __EXACT_FIT_FTOL = 1.0e-14
const __EXACT_FIT_STEPTOL = 1.0e-14
const __EXACT_FIT_RTOL = sqrt(eps(Float64))

# The exact solve's vector arithmetic, which does not depend on the quotes' types. Kept out of line, it is
# compiled once rather than inlined into the solve for each type of quotes. `__newton_step` is the Newton step
# for the scaled residuals `f`, from the Jacobian `J` of the unscaled ones, or `nothing` when the Jacobian is
# singular or not finite.
@noinline __line_search_scale(sizes) = map(__nonzero, sizes)
@noinline __scaled_residuals(v, search) = v ./ search
@noinline __newton_trial(z, α, δ) = z .- α .* δ
@noinline __relative_residuals(v, sizes) = map(__relative_residual, v, sizes)
@noinline function __newton_step(J, search, f)
    J = J ./ search
    F = lu(J; check = false)
    (all(isfinite, J) && issuccess(F)) || return nothing
    return F \ f
end
function __largest_relative_residual(v, s)
    largest = zero(float(eltype(s)))
    for i in eachindex(v, s)
        largest = max(largest, __relative_residual(v[i], s[i]))
    end
    return largest
end

# Newton's method on the repricing residuals `r` from the seed `z0`, for at most `maxiters` iterations. The
# line search on the sum of squared residuals takes each residual relative to its quote's size at the seed
# (`__quote_size`), so that it does not depend on notionals; the Newton steps do not depend on that scaling.
# It accepts only steps that decrease the sum. The stopping and acceptance rules take each residual relative
# to its quote's size on the curve reached, not at the seed: a long zero-coupon price far below its value at
# the seed would otherwise be repriced only to the seed's precision. The iteration stops at rounding:
# residuals within `__EXACT_FIT_FTOL`, a step within `__EXACT_FIT_STEPTOL`, or a step that no longer
# decreases the residuals.
function __solve_knot_rates(r::__RepricingResiduals, z0, maxiters)
    ReturnCode = Optimization.SciMLBase.ReturnCode
    v0, sizes = __repricing_point(r, z0)
    # The Jacobian's chunk follows the number of knots, so the configuration's type is known only at run time.
    # Asserting the Jacobian's type keeps the iterations inferred and compiled once for any number of knots
    # (the precompile workload's); only the Jacobian itself is compiled per chunk.
    config = ForwardDiff.JacobianConfig(r, z0)
    search = __line_search_scale(sizes)
    z, v, f = copy(z0), v0, __scaled_residuals(v0, search)
    ϕ = sum(abs2, f)
    retcode = ReturnCode.MaxIters
    for _ in 1:maxiters
        __largest_relative_residual(v, sizes) <= __EXACT_FIT_FTOL &&
            (retcode = ReturnCode.Success; break)
        δ = __newton_step(ForwardDiff.jacobian(r, z, config)::Matrix{eltype(v0)}, search, f)
        isnothing(δ) && (retcode = ReturnCode.InternalLinearSolveFailed; break)
        if maximum(abs, δ) <= __EXACT_FIT_STEPTOL
            z = __newton_trial(z, 1.0, δ)
            v, sizes = __repricing_point(r, z)
            f = __scaled_residuals(v, search)
            retcode = ReturnCode.Success
            break
        end
        α, accepted = 1.0, false
        for _ in 1:30
            zt = __newton_trial(z, α, δ)
            vt, st = __repricing_point(r, zt)
            ft = __scaled_residuals(vt, search)
            ϕt = sum(abs2, ft)
            if ϕt <= (1 - 1.0e-4 * α) * ϕ   # false for a non-finite trial
                z, v, f, ϕ, sizes, accepted = zt, vt, ft, ϕt, st, true
                break
            end
            α /= 2
        end
        accepted || (retcode = ReturnCode.Stalled; break)
    end
    final = __relative_residuals(v, sizes)
    worst = argmax(final)
    final[worst] <= __EXACT_FIT_RTOL && return z
    throw(__exact_fit_error(r.who, retcode, maxiters, z, final[worst], maturity(r.quotes[worst])))
end
function __exact_fit_error(who, retcode, maxiters, z, residual, maturity)
    ReturnCode = Optimization.SciMLBase.ReturnCode
    stopped = if retcode == ReturnCode.MaxIters
        "did not reprice its quotes within $maxiters Newton iterations"
    elseif retcode == ReturnCode.InternalLinearSolveFailed
        "stopped at knot rates $z, where the Jacobian of the quotes' pricing residuals is singular or not finite"
    else
        "cannot reduce its quotes' pricing residuals further"
    end
    return FitConvergenceError(
        retcode == ReturnCode.Success ? ReturnCode.Stalled : retcode,
        "$who $stopped: the largest relative residual is $residual, on the quote maturing at $maturity. No " *
            "curve with this interpolation may reprice the quotes (for example, a price that implies a " *
            "non-positive discount factor), or the quote prices may not determine every knot rate (a quote " *
            "whose cashflows do not depend on its own maturity's knot). fit(spline, quotes, Fit.Loss()) " *
            "minimizes the squared residuals instead."
    )
end

