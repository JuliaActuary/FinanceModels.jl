# ── Shared knot-grid construction ──────────────────────────────────────────────────────
#
# Every curve that interpolates zero rates over a knot grid (`Yield.Spline` and
# `Yield.MonotoneConvex`, however built: `ZeroRateCurve`, `reconstruct`, `fit`) obtains its
# knot data through `KnotGrid`, so all of them copy their inputs, promote to one concrete float
# type, and raise the same `ArgumentError`s for the same invalid grids.

# Minimum knot counts. Polynomial and B-spline orders reduce to `n - 1` on short grids, so one
# knot is a flat curve; PCHIP and Akima need three (two reach an UndefRefError inside
# DataInterpolations); a one-knot MonotoneConvex curve is flat.
__min_knots(::Sp.SplineCurve) = 1
__min_knots(::Sp.PCHIP) = 3
__min_knots(::Sp.Akima) = 3


# Owned, concretely-typed float `Vector` from any iterable of reals. Always copies (even
# when handed a `Vector`), so the curve never aliases caller-owned memory; never narrows
# AD types (`float(::Dual)` is a `Dual`). Inputs that are not real numbers fail where they are
# converted (`float(String)`, `Float64(::Dual)` for mixed dual tags, …). An untyped empty input
# is `Float64[]`, so the grid's own length check reports it.
function __owned_float_vector(x)
    v = x isa AbstractVector ? x : collect(x)
    return Vector{__float_eltype(v)}(v)   # Array-from-AbstractArray always allocates a fresh copy
end

# Knot rates as continuously compounded numbers: a number is taken as continuous and a `Rate` is
# converted. A `knot_rates` view gives its stored numbers directly, so `reconstruct(c)` makes no
# extra pass. Other inputs that are not all real numbers are converted element by element.
__continuous_knot_rates(x::ReadOnlyVector{<:FinanceCore.Rate{<:Any, Continuous}}) = getfield(x, :_data)
function __continuous_knot_rates(x)
    v = x isa AbstractVector ? x : collect(x)
    return eltype(v) <: Real ? v : map(__continuous, v)
end

# Marker for the internal, unvalidated `KnotGrid` construction used by optimizer trial curves.
struct Unchecked end

"""
    KnotGrid(rates, tenors, spline::Spline.SplineCurve; who = "KnotGrid")

Internal. The owned, validated knot data behind every `AbstractInterpolatedZeroCurve`
(`Yield.Spline`, `Yield.MonotoneConvex`). `rates` and `tenors` are each
copied from any iterable (a `Vector` is copied too, so no curve aliases caller-owned
memory) and promoted **independently** to one concrete floating-point element type
(`Int` → `Float64`, `Float32` + `BigFloat` → `BigFloat`, `Float64` + `ForwardDiff.Dual` →
`Dual`, never narrowed). Each rate is a number, read as continuously compounded, or a `Rate`,
converted to its continuously compounded value; a mixture is allowed. This is the one place
where knot rates are normalized.

Throws an `ArgumentError` (prefixed with `who`, the public constructor's name) when:

- `rates` and `tenors` differ in length, or either is empty;
- any rate or tenor is not finite (`NaN`, `±Inf`);
- any tenor is negative, or the tenors are not strictly increasing (unsorted or duplicated);
- there are fewer knots than `spline` needs (`__min_knots`).

Other inputs fail where they are converted (a `MethodError` such as `float(::Type{String})` for a
tenor, or a rate that is neither a number nor a `Rate`).

A curve built from a grid takes ownership of the grid's vectors. The grid is a construction
input, not a container to keep and modify.

`KnotGrid(Yield.Unchecked(), rates, tenors)` is internal: it neither copies nor validates. It
exists only for optimizer trial curves (`fit`, bootstrap) whose grid was validated up front
and whose candidate rates may be non-finite mid-search. Fitted results are always rebuilt
through the validating form.
"""
struct KnotGrid{R, T}
    rates::Vector{R}   # continuously-compounded zero rates (finite)
    tenors::Vector{T}  # finite, ≥ 0, strictly increasing
    KnotGrid(::Unchecked, rates::Vector{R}, tenors::Vector{T}) where {R, T} = new{R, T}(rates, tenors)
end
# unchecked trial form over non-`Vector` buffers (an optimizer may hand over a view or similar)
KnotGrid(u::Unchecked, rates::AbstractVector, tenors::AbstractVector) =
    KnotGrid(u, convert(Vector, rates), convert(Vector, tenors))

function KnotGrid(rates, tenors, spline::Sp.SplineCurve; who = "KnotGrid")
    r = __owned_float_vector(__continuous_knot_rates(rates))
    t = __owned_float_vector(tenors)
    length(r) == length(t) || throw(
        ArgumentError(
            "$who: `rates` and `tenors` must have the same length (got $(length(r)) and $(length(t)))."
        )
    )
    # the knot count first, so an empty grid reports it; then the tenors, which the rates of a
    # sampled grid depend on
    k = __min_knots(spline)
    length(t) >= k || throw(ArgumentError("$who: $(spline) requires at least $k knots (got $(length(t)))."))
    all(isfinite, t) || throw(ArgumentError("$who: all tenors must be finite (got $(t))."))
    # primal values throughout: ForwardDiff orders Duals with equal values by their partials, so
    # `Dual(0, -1) >= 0` is false and `Dual(1, 1) > Dual(1, 0)` is true
    __primal(first(t)) >= 0 || throw(ArgumentError("$who: tenors must be ≥ 0 (got $(t))."))
    # finiteness is checked first so NaN cannot defeat the ordering comparison
    for i in 2:length(t)
        __primal(t[i]) > __primal(t[i - 1]) || throw(
            ArgumentError(
                "$who: tenors must be strictly increasing (sorted, no duplicates); got $(t). " *
                    "Sort the (rate, tenor) pairs before constructing."
            )
        )
    end
    all(isfinite, r) || throw(ArgumentError("$who: all rates must be finite (got $(r))."))
    return KnotGrid(Unchecked(), r, t)
end
