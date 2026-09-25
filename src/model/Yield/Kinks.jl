# ── Kinks in the knot rates ────────────────────────────────────────────────────────────
# PCHIP, Akima, and MonotoneConvex are only piecewise smooth in their knot rates: their shape
# switches formula where a knot slope or forward crosses a threshold. This file is the one place
# that says where those switches are; other code uses two entry points:
#
# - `__build_public(spline, grid)`: the public constructors' build, which throws when the
#   caller's dual knot rates move a PCHIP or Akima kink;
# - `__near_kink(curve, z[, δ])`: whether knots `z` lie on a kink (a differentiated fit's check).
#
# `__kink_quantities(spline, z, tenors)` returns the scalar functions of the knot rates `z` whose
# zeros are the switches (none for interpolants that are linear in the knot rates).
# MonotoneConvex's sit with its formulas in `MonotoneConvex.jl`, and it has derivatives defined at
# its kinks. PCHIP's and Akima's copy DataInterpolations' branch rules, at the end of this file.
__kink_quantities(::Sp.SplineCurve, z, tenors) = eltype(z)[]
__kink_quantities(c::AbstractInterpolatedZeroCurve, z) = __kink_quantities(c.spline, z, c.tenors)

# The rounding scale of forwards and slopes derived from knots `(z, tenors)`: differences of t·z
# over the shortest interval (from 0 to the first knot, or between knots; a first knot at t = 0
# leaves no interval before it). Kink quantities within 16 times this count as zero, since
# rounding alone can move a flat or straight run of knots that far off its kink.
function __knot_noise(z, tenors)
    zmax = float(maximum(x -> abs(__primal(x)), z))
    dt = minimum(filter(>(0), diff([zero(first(tenors)); tenors])))
    return eps(zmax) * max(one(zmax), last(tenors)) / dt
end

# Whether the knots `z` of `curve` lie on a kink, up to rounding. Given a step `δ` towards an
# exact solution (a fit's Newton correction), also whether the step cannot resolve which side of
# a kink that solution lies on: its change in a kink quantity is at least the quantity itself.
function __near_kink(c::AbstractInterpolatedZeroCurve, z, δ = zero(z))
    tol = 16 * __knot_noise(z, c.tenors)
    q0 = __kink_quantities(c, z)
    q1 = __kink_quantities(c, z .- δ)
    return any(i -> abs(q1[i]) <= abs(q0[i] - q1[i]) + tol, eachindex(q0, q1))
end

# The public constructors (`ZeroRateCurve`, `Yield.Spline`, `reconstruct`) check the caller's dual
# knot rates. FinanceModels' own trial curves (optimizer candidates, bootstrap steps, the
# calibration Jacobian) do not: an optimizer only needs some slope at a kink, and a
# differentiated fit checks its fitted knots itself before its Jacobian.
function __build_public(s::Sp.SplineCurve, g::KnotGrid; extrapolation = :flat_forward)
    __check_dual_kinks(s, g)
    return __build(s, g; extrapolation)
end

__check_dual_kinks(::Sp.SplineCurve, g::KnotGrid) = nothing   # smooth, or MonotoneConvex's own semantics
function __check_dual_kinks(s::Union{Sp.PCHIP, Sp.Akima}, g::KnotGrid)
    eltype(g.rates) <: ForwardDiff.Dual || return nothing
    tol = 16 * __knot_noise(g.rates, g.tenors)
    for q in __kink_quantities(s, g.rates, g.tenors)
        (abs(ForwardDiff.value(q)) <= tol && !iszero(ForwardDiff.partials(q))) && throw(
            ArgumentError(
                "$(nameof(typeof(s))) interpolation cannot differentiate with respect to its knot rates at " *
                    "these knots: its shape switches formula here (for example, at a flat segment or a " *
                    "straight run of knots), so the derivative depends on the bump direction. Use " *
                    "Spline.Linear(), Spline.Cubic(), or Spline.MonotoneConvex(), key-rate shifts added to " *
                    "the curve (ActuaryUtilities `KeyRates`), or finite bumps."
            )
        )
    end
    return nothing
end

# ── DataInterpolations' branch rules ───────────────────────────────────────────────────
# PCHIP and Akima evaluate through DataInterpolations, which at these switches returns a
# one-sided value, NaN, or the derivative of a fallback formula (and Akima's value can jump).
# These copies of its internal branch conditions match DataInterpolations 9.2 through 10.1; the
# "DataInterpolations branch rules" tests in `test/kinks.jl` check them against the installed
# version, and must be revisited with any change to that compat bound or upstream fix.

# `du_PCHIP`: node slopes branch on the signs of the secant slopes δ, and each end slope on the
# sign of `d` and on `|d| > 3|δ₁|` when δ₁ and δ₂ differ in sign.
function __kink_quantities(::Sp.PCHIP, z, tenors)
    h = diff(tenors)
    δ = diff(z) ./ h
    q = collect(δ)
    for (h₁, h₂, δ₁, δ₂) in ((h[1], h[2], δ[1], δ[2]), (h[end], h[end - 1], δ[end], δ[end - 1]))
        d = ((2 * h₁ + h₂) * δ₁ - h₁ * δ₂) / (h₁ + h₂)
        push!(q, d, sign(δ₁) != sign(δ₂) ? abs(d) - 3 * abs(δ₁) : one(d))
    end
    return q
end

# `_akima_init!`: node slopes weight the secant slopes (padded by linear extrapolation at each
# end) by `abs` of their differences, and a node whose total weight is at most 1e-9 of the
# largest takes a fallback slope instead: the curve's value jumps there.
function __kink_quantities(::Sp.Akima, z, tenors)
    m = diff(z) ./ diff(tenors)
    n = length(m)
    m2 = 2 * m[1] - m[2]
    mn = 2 * m[n] - m[n - 1]
    Δ = diff([2 * m2 - m[1]; m2; m; mn; 2 * mn - m[n]])
    w = [abs(Δ[i]) + abs(Δ[i + 2]) for i in 1:(n + 1)]
    return [Δ; w .- 1.0e-9 * maximum(w)]
end
