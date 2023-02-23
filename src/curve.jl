# default = Bootstrap
function curve(method::T,quotes) where {T<:CurveMethod}
    method(quotes)
end

function curve(method::V,quotes::Vector{Quote{U,Forward{N,T}}}) where {N,T<:Cashflow,U, V<:CurveMethod}
    obs = __process_forwards(quotes)
    curve(method,obs)
end

function curve(quotes) where {T<:CurveMethod}
    method = Bootstrap()
    curve(method,quotes)
end