# Excel publication arithmetic has 15 significant decimal digits. Normalizing
# that decimal value before rounding avoids binary products falling just below
# a half-bp boundary (e.g. 0.85 * 0.007). Keep the subsequent tie decision exact;
# converting the normalized value back to Float64 first reintroduces the bug.
function _excel_round_exact(x::AbstractFloat, digits::Integer)
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

_excel_round(x::AbstractFloat, digits::Integer) = _excel_round_exact(x, digits)
function _excel_round(x::Float64, digits::Integer)
    (iszero(x) || !isfinite(x)) && return x
    # 10^0 through 10^22 are exact Float64 values, so reconstructing the rounded
    # decimal by division gives the same float as the exact rational conversion.
    # At 23 digits even 1e-23 can otherwise return the adjacent float.
    if 0 <= digits <= 22
        scale = 10.0^digits
        y = x * scale
        # Below 1e8, multiplication error and 15-digit decimal normalization
        # together move y by less than 1e-6. Outside that margin around a half
        # integer they cannot change its rounded integer. Ties use the exact path.
        if abs(y) < 1.0e8 && abs(abs(y - trunc(y)) - 0.5) > 1.0e-6
            return round(y, RoundNearestTiesAway) / scale
        end
    end
    return _excel_round_exact(x, digits)
end

_round_rate(x) = _excel_round(x, 5)
_round_va(x) = _excel_round(x, 4)
