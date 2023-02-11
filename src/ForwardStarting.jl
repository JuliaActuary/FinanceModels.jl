"""
    ForwardStarting(curve,forwardstart)

Rebase a `curve` so that `discount`/`accumulation`/etc. are re-based so that time zero from the new curves perspective is the given `forwardstart` time.

# Examples

```julia-repl
julia> zero = [5.0, 5.8, 6.4, 6.8] ./ 100
julia> maturity = [0.5, 1.0, 1.5, 2.0]
julia> curve = Yields.Zero(zero, maturity)
julia> fwd = Yields.ForwardStarting(curve, 1.0)

julia> FinanceCore.discount(curve,1,2)
0.9275624570410582

julia> FinanceCore.discount(fwd,1) # `curve` has effectively been reindexed to `1.0`
0.9275624570410582
```

# Extended Help

While `ForwardStarting` could be nested so that, e.g. the third period's curve is the one-period forward of the second period's curve, it will be more efficient to reuse the initial curve from a runtime and compiler perspective.

`ForwardStarting` is not used to construct a curve based on forward rates. See  [`Forward`](@ref) instead.
"""
struct ForwardStarting{T,U} <: AbstractYieldCurve
    curve::U
    forwardstart::T
end
__ratetype(::Type{ForwardStarting{T,U}}) where {T,U}= __ratetype(U)

function FinanceCore.discount(c::ForwardStarting, to)
    FinanceCore.discount(c.curve, c.forwardstart, to + c.forwardstart)
end

function Base.zero(c::ForwardStarting, to,cf::C) where {C<:FinanceCore.CompoundingFrequency}
    z = forward(c.curve,c.forwardstart,to+c.forwardstart)
    return convert(cf,z)
end
