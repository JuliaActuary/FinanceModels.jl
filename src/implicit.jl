# Implicit differentiation of solved quantities.
#
# FinanceModels differentiates roots and calibrated parameters with the implicit
# function theorem rather than through solver iterations. Internal derivatives are
# taken on primal (non-Dual) values only; a caller's ForwardDiff partials pass
# through plain arithmetic. An outer tag therefore never meets an internal one,
# whose relative order ForwardDiff cannot determine reliably.

__ad_depth(::Type) = 0
__ad_depth(::Type{<:ForwardDiff.Dual{T, V}}) where {T, V} = 1 + __ad_depth(V)
__ad_depth(x) = __ad_depth(typeof(x))

__primal(x) = x
__primal(x::ForwardDiff.Dual) = __primal(ForwardDiff.value(x))

"""
    __implicit_root(g, g_primal, x0; bracket = (-1.0, 1.0), who = "root", slope = nothing)

Solve `g_primal(x) = 0` from `x0`, falling back to a bracketed solve on `bracket`,
and return the root with first-order ForwardDiff partials of `g` propagated by the
implicit function theorem: `dx = -(∂g/∂θ) / (∂g/∂x)`.

`g_primal` must be `g` evaluated without dual numbers, or `g` with its dual numbers
stripped from the result. `slope(x)`, when given, is `∂g/∂x` on primal values in closed
form; otherwise it is `ForwardDiff.derivative(g_primal, x)`, which requires `g_primal`
to involve no dual numbers at all. The value of the result is the primal root exactly;
its partials come from one dual correction step. Nested dual numbers and a vanishing
slope throw an `ArgumentError` naming `who`.
"""
function __implicit_root(g::G, g_primal::P, x0; bracket = (-1.0, 1.0), who = "root", slope::S = nothing) where {G, P, S}
    x = try
        Roots.find_zero(g_primal, float(x0), Roots.Order1())
    catch e
        # Only a convergence failure falls back to a bracketed solve; errors raised
        # while evaluating the residual must surface.
        e isa Roots.ConvergenceFailed || rethrow()
        Roots.find_zero(g_primal, bracket, Roots.A42())
    end
    gx = g(x)
    depth = __ad_depth(gx)
    depth == 0 && return x
    depth == 1 || throw(ArgumentError("$who supports first-order ForwardDiff derivatives only; nested dual numbers are not supported"))
    s = slope === nothing ? ForwardDiff.derivative(g_primal, x) : slope(x)
    (isfinite(s) && abs(s) > sqrt(eps(typeof(s)))) ||
        throw(ArgumentError("$who has a vanishing or non-finite derivative at the solution ($s); its sensitivity is undefined"))
    # `gx` has a primal value of zero at the root, so this step keeps the root's
    # value and carries only the implicit-function partials.
    return x - (gx - __primal(gx)) / s
end
