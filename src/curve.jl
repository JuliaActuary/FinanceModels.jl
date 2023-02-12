# default = Bootstrap
function curve(method::T,instruments) where {T<:CurveMethod}
    method(instruments)
end

function curve(instruments) where {T<:CurveMethod}
    method = Bootstrap()
    method(instruments)
end