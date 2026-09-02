"""
    CairnsPritchard(c₁, c₂, b₀, b₁, b₂)
    CairnsPritchard(c₁=0.5, c₂=3.0) # used in fitting

A Cairns-Pritchard yield curve model with 2 exponential components.

The continuous zero rate at time `t` is:

``r(t) = b₀ + b₁ \\exp(-c₁ t) + b₂ \\exp(-c₂ t)``

This is a generalization of Nelson-Siegel with independent decay rates per
exponential component. Parameters and default fitting bounds:

- `c₁` decay rate for first component: `0.001 .. 10.0`
- `c₂` decay rate for second component: `0.001 .. 10.0`
- `b₀` long-term rate level: `-1.0 .. 1.0`
- `b₁` first exponential coefficient: `-10.0 .. 10.0`
- `b₂` second exponential coefficient: `-10.0 .. 10.0`

See also [`CairnsPritchardExtended`](@ref) for a 3-component variant.

# References
- Cairns, A.J.G. (1998). "Descriptive Bond-Yield and Forward-Rate Models for the British Government Securities Market". British Actuarial Journal, 4(2), 265-321.
"""
struct CairnsPritchard{T} <: AbstractYieldModel
    c₁::T
    c₂::T
    b₀::T
    b₁::T
    b₂::T

    function CairnsPritchard(c₁::T, c₂::T, b₀::T, b₁::T, b₂::T) where {T}
        (c₁ <= 0 || c₂ <= 0) && throw(DomainError("Decay parameters c must be positive"))
        return new{T}(c₁, c₂, b₀, b₁, b₂)
    end
end

# Promote mixed argument types for ForwardDiff compatibility
function CairnsPritchard(c₁, c₂, b₀, b₁, b₂)
    T = promote_type(typeof(c₁), typeof(c₂), typeof(b₀), typeof(b₁), typeof(b₂))
    return CairnsPritchard(convert(T, c₁), convert(T, c₂), convert(T, b₀), convert(T, b₁), convert(T, b₂))
end

# Default constructor with different c values to break symmetry during fitting.
# Non-zero b values provide a reasonable starting curve for the optimizer.
CairnsPritchard(c₁ = 0.5, c₂ = 3.0) = CairnsPritchard(c₁, c₂, 0.05, -0.01, -0.01)

function Base.zero(cp::CairnsPritchard, t)
    # At t=0 the formula is already well-defined: exp(0) = 1, so z(0) = b₀ + b₁ + b₂
    return Continuous(cp.b₀ + cp.b₁ * exp(-cp.c₁ * t) + cp.b₂ * exp(-cp.c₂ * t))
end
FinanceCore.discount(cp::CairnsPritchard, t) = _discount_from_zero(cp, t)

"""
    CairnsPritchardExtended(c₁, c₂, c₃, b₀, b₁, b₂, b₃)
    CairnsPritchardExtended(c₁=0.5, c₂=2.0, c₃=5.0) # used in fitting

A Cairns-Pritchard yield curve model with 3 exponential components.

The continuous zero rate at time `t` is:

``r(t) = b₀ + b₁ \\exp(-c₁ t) + b₂ \\exp(-c₂ t) + b₃ \\exp(-c₃ t)``

Parameters and default fitting bounds:

- `c₁` decay rate for first component: `0.001 .. 10.0`
- `c₂` decay rate for second component: `0.001 .. 10.0`
- `c₃` decay rate for third component: `0.001 .. 10.0`
- `b₀` long-term rate level: `-1.0 .. 1.0`
- `b₁` first exponential coefficient: `-10.0 .. 10.0`
- `b₂` second exponential coefficient: `-10.0 .. 10.0`
- `b₃` third exponential coefficient: `-10.0 .. 10.0`

See also [`CairnsPritchard`](@ref) for a 2-component variant.

# References
- Cairns, A.J.G. (1998). "Descriptive Bond-Yield and Forward-Rate Models for the British Government Securities Market". British Actuarial Journal, 4(2), 265-321.
"""
struct CairnsPritchardExtended{T} <: AbstractYieldModel
    c₁::T
    c₂::T
    c₃::T
    b₀::T
    b₁::T
    b₂::T
    b₃::T

    function CairnsPritchardExtended(c₁::T, c₂::T, c₃::T, b₀::T, b₁::T, b₂::T, b₃::T) where {T}
        (c₁ <= 0 || c₂ <= 0 || c₃ <= 0) && throw(DomainError("Decay parameters c must be positive"))
        return new{T}(c₁, c₂, c₃, b₀, b₁, b₂, b₃)
    end
end

# Promote mixed argument types for ForwardDiff compatibility
function CairnsPritchardExtended(c₁, c₂, c₃, b₀, b₁, b₂, b₃)
    T = promote_type(typeof(c₁), typeof(c₂), typeof(c₃), typeof(b₀), typeof(b₁), typeof(b₂), typeof(b₃))
    return CairnsPritchardExtended(convert(T, c₁), convert(T, c₂), convert(T, c₃), convert(T, b₀), convert(T, b₁), convert(T, b₂), convert(T, b₃))
end

# Default constructor with different c values to break symmetry during fitting.
# Non-zero b values provide a reasonable starting curve for the optimizer.
CairnsPritchardExtended(c₁ = 0.5, c₂ = 2.0, c₃ = 5.0) = CairnsPritchardExtended(c₁, c₂, c₃, 0.05, -0.01, -0.01, -0.01)

function Base.zero(cp::CairnsPritchardExtended, t)
    # At t=0 the formula is already well-defined: exp(0) = 1, so z(0) = b₀ + b₁ + b₂ + b₃
    return Continuous(cp.b₀ + cp.b₁ * exp(-cp.c₁ * t) + cp.b₂ * exp(-cp.c₂ * t) + cp.b₃ * exp(-cp.c₃ * t))
end
FinanceCore.discount(cp::CairnsPritchardExtended, t) = _discount_from_zero(cp, t)
