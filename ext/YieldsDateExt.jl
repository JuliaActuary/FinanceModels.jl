module YieldsDateExt

using Dates
using DayCounts

struct DateCurve{T,D<:DayCounts.AbstractDayCount}
    curve::T
    dc::D
    reference_date::Date

    function new(curve::T,dc::D=DayCounts.ActualActual365, reference_date::Date=today() ) where {T}
        return new{T,D}(curve, dc, reference_date)
    end
end

discount(c::DateCurve, date::Date) = discount(c.curve, yearfrac(c.reference_date,date, c.dc))

end