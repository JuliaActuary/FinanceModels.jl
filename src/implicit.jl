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

# A primal solve evaluates only its own dual numbers. A caller's dual number that reaches it (inside
# a contract type a fit does not strip, or in a model field it does not optimize) would mix with the
# solver's iterations and give the derivative of the iterations, often zero. Each solve checks what
# it evaluates, every evaluation, for a ForwardDiff tag that is not one of its `own`.
# `Own` is a `Tuple` type of the solve's tags, so the check is decided at compile time.
__foreign_dual(::Type, ::Type{Own}) where {Own <: Tuple} = false
__foreign_dual(::Type{ForwardDiff.Dual{T, V, N}}, ::Type{Own}) where {T, V, N, Own <: Tuple} =
    !(T in fieldtypes(Own)) || __foreign_dual(V, Own)
__foreign_dual_error(who, hint) =
    ArgumentError(rstrip("$who: a ForwardDiff dual number reached the solve through data it does not differentiate. $hint"))

# Whether `x` holds a ForwardDiff dual number anywhere. A solve's data is checked before it starts: a
# dual number in it may reach the evaluated values only away from the start, where the check of each
# evaluation would not see it until the solver went there (or never). The types decide where they can
# (`__dual_free`, a constant per type); otherwise the values are visited, each mutable one once.
__holds_dual(x) = !__dual_free(typeof(x)) && __holds_dual(x, Base.IdSet{Any}())
function __holds_dual(x, seen)
    x isa ForwardDiff.Dual && return true
    (__dual_free(typeof(x)) || x isa Union{Type, Module, Core.SimpleVector}) && return false
    if ismutable(x)
        x in seen && return false
        push!(seen, x)
    end
    x isa AbstractArray && return any(i -> isassigned(x, i) && __holds_dual(x[i], seen), eachindex(x))
    return any(i -> isdefined(x, i) && __holds_dual(getfield(x, i), seen), 1:nfields(x))
end

# Whether no value of type `T` can hold a dual number: a concrete type whose fields (for an array, its
# elements) are of such types. Each type is visited once, so a type that contains itself ends.
function __dual_free_impl(T, visited)
    haskey(visited, T) && return visited[T]
    visited[T] = true   # a type that contains itself: its other fields decide
    free = !(T <: ForwardDiff.Dual) && isconcretetype(T) && (
        T <: AbstractArray ? __dual_free_impl(eltype(T), visited) : all(F -> __dual_free_impl(F, visited), fieldtypes(T))
    )
    visited[T] = free
    return free
end
@generated __dual_free(::Type{T}) where {T} = __dual_free_impl(T, IdDict{Any, Bool}())

# `f(x)` in a primal solve whose own tags are `Own`. `x` (element type `X`) is checked before `f`
# runs, and the value after, before any conversion. A dual number captured untyped can carry a tag of
# an own type, so the value also may not carry more dual levels than `x`.
function __own_duals(f::F, x, ::Type{X}, ::Type{Own}, who, hint) where {F, X, Own <: Tuple}
    __foreign_dual(X, Own) && throw(__foreign_dual_error(who, hint))
    v = f(x)
    (__foreign_dual(typeof(v), Own) || __ad_depth(v) > __ad_depth(X)) && throw(__foreign_dual_error(who, hint))
    return v
end

# A primal residual: a bootstrap step, or an implicit root's `g_primal`. Its only own tag is the one
# `ForwardDiff.derivative(r, x)` creates for the implicit slope, keyed by the residual's type.
struct __PrimalResidual{F}
    f::F
    who::String
    hint::String
end
(r::__PrimalResidual)(x) = __own_duals(r.f, x, typeof(x), Tuple{ForwardDiff.Tag{typeof(r), typeof(__primal(x))}}, r.who, r.hint)

# The primal root search of `__implicit_root` and of each bootstrap step: the secant method from
# `x0`, and only if it fails to converge a bracketed search over `bracket`. Errors raised while
# evaluating `g` (such as a quote that cannot be priced) surface rather than being retried. (`G`
# makes Julia specialize on `g`, which is only passed on to `find_zero`.)
function __solve_primal_root(g::G, x0, bracket) where {G}
    return try
        Roots.find_zero(g, float(x0), Roots.Order1())
    catch e
        e isa Roots.ConvergenceFailed || rethrow()
        Roots.find_zero(g, bracket, Roots.A42())
    end
end

"""
    __implicit_root(g, g_primal, x0; bracket = (-1.0, 1.0), who, hint = "", slope = nothing, scale)

Solve `g_primal(x) = 0` from `x0`, falling back to a bracketed solve on `bracket`,
and return the root with first-order ForwardDiff partials of `g` propagated by the
implicit function theorem: `dx = -(∂g/∂θ) / (∂g/∂x)`.

`g_primal` must be `g` evaluated without dual numbers, or `g` with its dual numbers
stripped from the result. `slope(x)`, when given, is `∂g/∂x` on primal values in closed
form; otherwise it is the ForwardDiff derivative of `g_primal`. `scale(x)` is the magnitude, free of units such as a
notional, that both checks at the solution use: the residual at an accepted root must be at
most `sqrt(eps)` times it, and the slope must exceed `sqrt(eps)` times it. The caller chooses
it for its residual: `implied_quote` passes the size of the quote's price and value, the
swaption critical rate the size of the terms of its slope. The solvers stop on absolute
tolerances, so `g` should already be expressed relative to such a scale. The value of the
result is the primal root exactly; its partials come from one dual correction step. Nested
dual numbers, a dual number in `g_primal`'s argument or value, a root that does not solve
`g_primal`, and a vanishing slope throw an `ArgumentError` naming `who` (and `hint`).
"""
function __implicit_root(g::G, g_primal::P, x0; bracket = (-1.0, 1.0), who, hint = "", slope::S = nothing, scale::C) where {G, P, S, C}
    residual = __PrimalResidual(g_primal, who, hint)
    x = __solve_primal_root(residual, x0, bracket)
    # A solver can accept a point whose residual is small in absolute terms only.
    tol = sqrt(eps(float(typeof(x)))) * scale(x)
    abs(residual(x)) <= tol || throw(
        ArgumentError("$who did not converge: the residual at $x is $(residual(x)), more than $tol")
    )
    gx = g(x)
    depth = __ad_depth(gx)
    depth == 0 && return x
    depth == 1 || throw(ArgumentError("$who supports first-order ForwardDiff derivatives only; nested dual numbers are not supported"))
    s = slope === nothing ? ForwardDiff.derivative(residual, x) : slope(x)
    (isfinite(s) && abs(s) > sqrt(eps(typeof(s))) * scale(x)) ||
        throw(ArgumentError("$who has a vanishing or non-finite derivative at the solution ($s); its sensitivity is undefined"))
    # `gx` has a primal value of zero at the root, so this step keeps the root's
    # value and carries only the implicit-function partials.
    return x - (gx - __primal(gx)) / s
end
