abstract type CurveMethod end


# default = Bootstrap
function curve(method::T=Bootstrap(),instruments) where {T<:CurveMethod}
    method(instruments)
end