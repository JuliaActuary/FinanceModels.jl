# Excel publication arithmetic has 15 significant decimal digits. Normalizing
# that decimal value before rounding avoids binary products falling just below
# a half-bp boundary (e.g. 0.85 * 0.007). Keep the subsequent tie decision exact;
# converting the normalized value back to Float64 first reintroduces the bug.
function _excel_round(x::AbstractFloat, digits::Integer)
    (iszero(x) || !isfinite(x)) && return x
    value = Rational{BigInt}(x)
    ten = BigInt(10) // BigInt(1)
    exponent = floor(Int, log10(abs(x)))
    # Correct log10's approximation at powers of ten using exact comparisons.
    while abs(value) < ten^exponent
        exponent -= 1
    end
    while abs(value) >= ten^(exponent + 1)
        exponent += 1
    end
    scale = ten^(14 - exponent)
    decimal = round(BigInt, value * scale) / scale
    unit = ten^digits
    rounded = round(BigInt, decimal * unit, RoundNearestTiesAway)
    return iszero(rounded) ? copysign(zero(x), x) : oftype(x, rounded / unit)
end

_round_rate(x) = _excel_round(x, 5)
_round_va(x) = _excel_round(x, 4)
