# ── Differentiable knot-curve fits ─────────────────────────────────────────────────────
#
# A spline `fit` whose quotes carry ForwardDiff dual numbers (for example
# `fit(Spline.Linear(), OISYield.(dual_rates, tenors), Fit.Bootstrap())` inside a gradient)
# returns a curve whose knot rates carry the exact first-order derivatives of the
# calibration. With one knot per quote, the fitted knot rates `z` solve the square system
#
#     Rᵢ(z, p) = present_value(curve(z), qᵢ.instrument) - qᵢ.price = 0,
#
# so by the implicit function theorem `dz = -(∂R/∂z)⁻¹ (∂R/∂p) dp`. The fit itself runs on
# primal copies of the quotes, `∂R/∂z` is a ForwardDiff Jacobian taken on primal values, and
# the caller's partials enter only through one evaluation of `R` at the fitted rates with
# plain arithmetic. FinanceModels' own derivatives therefore never meet the caller's tag.

function __no_dual_time(t)
    __ad_depth(t) == 0 || throw(
        ArgumentError(
            "fit cannot differentiate with respect to quote maturities or cashflow times (got $t); " *
                "only quote prices, rates, and cashflow amounts may carry dual numbers."
        )
    )
    return nothing
end

# Primal copy of a contract: dual numbers in prices, rates, and amounts are replaced by their
# values. A contract that carries none is returned itself (`===`), so primal calibrations
# take exactly their former path. Unknown contract types are returned unchanged; a dual
# number hidden in one is caught when the fit checks its primal residuals.
__primal_contract(x) = x
function __primal_contract(c::FinanceCore.Cashflow)
    __no_dual_time(c.time)
    a = __primal(c.amount)
    return a === c.amount ? c : FinanceCore.Cashflow(a, c.time)
end
function __primal_contract(c::FinanceCore.Composite)
    a, b = __primal_contract(c.a), __primal_contract(c.b)
    return (a === c.a && b === c.b) ? c : FinanceCore.Composite(a, b)
end
function __primal_contract(b::Bond.Fixed)
    __no_dual_time(b.maturity)
    r = __primal(b.coupon_rate)
    return r === b.coupon_rate ? b : Bond.Fixed(r, b.frequency, b.maturity)
end
function __primal_contract(b::Bond.Floating)
    __no_dual_time(b.maturity)
    r = __primal(b.coupon_rate)
    return r === b.coupon_rate ? b : Bond.Floating(r, b.frequency, b.maturity, b.key)
end
function __primal_contract(c::FX.BasisSwapLeg)
    cfs = __primal_contract(c.cashflows)
    return cfs === c.cashflows ? c : FX.BasisSwapLeg(c.pair, cfs)
end
function __primal_contract(v::Union{AbstractVector, Tuple})
    p = map(__primal_contract, v)
    return all(((a, b),) -> a === b, zip(p, v)) ? v : p
end

function __primal_quote(q::FinanceCore.Quote)
    p, i = __primal(q.price), __primal_contract(q.instrument)
    return (p === q.price && i === q.instrument) ? q : FinanceCore.Quote(p, i)
end

__primal_extrapolation(e) = e
__primal_extrapolation(e::Yield.FlatForwardAt) =
    __ad_depth(e.forward) == 0 ? e : Yield.FlatForwardAt(Continuous(__primal(e.forward)))

# The quotes of a knot-curve fit, their primal copies, and whether any carried dual numbers.
# Knots sit at the quote maturities, which therefore must be primal.
function __calibration_quotes(quotes)
    qs = collect(quotes)
    foreach(q -> __no_dual_time(FinanceCore.maturity(q)), qs)
    qp = map(__primal_quote, qs)
    return qs, qp, any(i -> qp[i] !== qs[i], eachindex(qs))
end

__quotes_carry_ad(quotes) = any(q -> __primal_quote(q) !== q, quotes)

# Before solving: a dual number that survived `__primal_contract` (inside a contract type it
# does not know) would reach the solver and mix with its internal derivatives.
function __check_primal_quotes(curve, quotes)
    for q in quotes
        r = present_value(curve, q.instrument) - q.price
        __ad_depth(r) == 0 || throw(
            ArgumentError(
                "fit cannot differentiate through a quote on $(nameof(typeof(q.instrument))): dual numbers are " *
                    "supported in quote prices and in the rates and amounts of Cashflow, Composite, Bond.Fixed, " *
                    "Bond.Floating, and FX.BasisSwapLeg instruments."
            )
        )
    end
    return nothing
end

# The one dual type shared by `xs` (nothing when all are primal).
function __common_dual_type(xs)
    D = nothing
    for x in xs
        d = __ad_depth(x)
        d == 0 && continue
        d == 1 || throw(
            ArgumentError("fit supports first-order ForwardDiff derivatives only; nested dual numbers are not supported.")
        )
        T = typeof(x)
        if D === nothing
            D = T
        elseif D !== T
            throw(
                ArgumentError(
                    "the quotes carry dual numbers from different ForwardDiff calls ($D and $T); " *
                        "differentiate all calibration inputs in one call."
                )
            )
        end
    end
    return D
end

# A fit is differentiated only when it reprices its quotes: the conditions `R = 0` then hold
# at the fitted rates. Bootstrap and exact loss fits reprice to solver precision. The fit's error
# is estimated by one Newton correction of the knot rates (the largest absolute component), which
# does not depend on notionals.
const __IMPLICIT_FIT_RTOL = 1.0e-6

"""
    __implicit_knot_curve(curve, quotes, primal_quotes, extrapolation, has_ad)

Given the knot curve `curve` fitted to `primal_quotes` (one knot per quote), return it with
knot rates carrying the first-order derivatives of the calibration with respect to the dual
numbers in `quotes` and in `extrapolation`. Returns `curve` itself when there are none.
"""
function __implicit_knot_curve(curve, quotes, primal_quotes, extrapolation, has_ad)
    (has_ad || any(x -> __ad_depth(x) > 0, __extrapolation_dual(extrapolation))) || return curve
    residuals(c, qs) = [present_value(c, q.instrument) - q.price for q in qs]
    # R at the fitted rates, carrying the caller's partials: quotes and the extrapolation
    # policy are the only places they can enter.
    dual_curve = extrapolation === curve.extrapolation ? curve : Yield.reconstruct(curve; extrapolation)
    Rd = residuals(dual_curve, quotes)
    D = __common_dual_type((Rd..., __extrapolation_dual(extrapolation)...))
    D === nothing && return curve

    z0 = collect(Yield.knot_rates(curve))
    tenors = collect(Yield.knot_tenors(curve))
    length(Rd) == length(z0) || throw(
        ArgumentError("fit can differentiate only calibrations with one knot per quote (got $(length(Rd)) quotes and $(length(z0)) knots).")
    )
    # A kink of the interpolant at the fitted knots (flat quotes make adjacent forwards equal,
    # for example) leaves the refit without a derivative. Checked first: the Jacobian below
    # would otherwise meet the interpolant's own check.
    Yield.__near_kink(curve, z0) && throw(__fit_kink_error(curve))

    s, e = curve.spline, curve.extrapolation
    build(z) = Yield.__build(s, Yield.KnotGrid(Yield.Unchecked(), z, tenors); extrapolation = e)
    A = ForwardDiff.jacobian(z -> residuals(build(z), primal_quotes), z0)
    # Equilibrate each residual by its sensitivity to the knot rates, so that the conditioning
    # check does not depend on the instruments' notionals.
    w = [maximum(abs, view(A, i, :)) for i in axes(A, 1)]
    singular(κ) = ArgumentError(
        "fit cannot differentiate this calibration: the quote prices do not determine the knot rates " *
            "(the repricing Jacobian is singular or ill-conditioned, condition number $κ). " *
            "Check for quotes whose cashflows do not depend on their own maturity's knot."
    )
    (all(isfinite, A) && all(>(0), w)) || throw(singular(Inf))
    As = A ./ w
    F = lu(As; check = false)
    (issuccess(F) && cond(As) <= 1.0e10) || throw(singular(cond(As)))
    r0 = residuals(curve, primal_quotes)
    # One Newton correction of the fitted knots: a local estimate of the fit's error in knot rates.
    δ = F \ (r0 ./ w)
    worst = argmax(abs.(δ))
    abs(δ[worst]) <= __IMPLICIT_FIT_RTOL || throw(
        ArgumentError(
            "fit cannot differentiate a curve that does not reprice its quotes: a Newton correction " *
                "towards the exact fit moves its knot rate at $(tenors[worst]) by $(abs(δ[worst])) " *
                "(largest quote residual $(maximum(abs, r0))), more than 1e-6. Its derivatives assume " *
                "an exact fit; tighten the optimizer with `solve_kwargs` (for example " *
                "`solve_kwargs = (; g_tol = 1e-12)`), or use Fit.Bootstrap() with Spline.Linear()."
        )
    )
    # The fitted knots must resolve which side of each kink the exact fit lies on.
    Yield.__near_kink(curve, z0, δ) && throw(__fit_kink_error(curve))

    K = ForwardDiff.npartials(D)
    P = [ForwardDiff.partials(convert(D, Rd[i]), k) / w[i] for i in eachindex(Rd), k in 1:K]
    dZ = -(F \ P)
    W = promote_type(ForwardDiff.valtype(D), eltype(z0), eltype(dZ))
    Dz = ForwardDiff.Dual{ForwardDiff.tagtype(D), W, K}
    zd = [Dz(W(z0[i]), ForwardDiff.Partials{K, W}(ntuple(k -> W(dZ[i, k]), K))) for i in eachindex(z0)]
    return Yield.reconstruct(curve; rates = zd, extrapolation)
end

__fit_kink_error(curve) = ArgumentError(
    "fit cannot differentiate this calibration: the fitted $(nameof(typeof(curve.spline))) curve lies on, or " *
        "within the fit's precision of, a point where the interpolation switches shape (for example, flat " *
        "quotes make adjacent forwards equal), so the refitted curve has no derivative with respect to the " *
        "quotes there. Use Spline.Linear() or Spline.Cubic(), or finite bumps."
)

__extrapolation_dual(e) = ()
__extrapolation_dual(e::Yield.FlatForwardAt) = (e.forward,)
