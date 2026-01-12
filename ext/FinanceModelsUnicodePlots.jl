module FinanceModelsUnicodePlots
import FinanceModels
import UnicodePlots
import FinanceCore
# used to display simple type name in show method
# https://stackoverflow.com/questions/70043313/get-simple-name-of-type-in-julia?noredirect=1#comment123823820_70043313
name(::Type{T}) where {T} = (isempty(T.parameters) ? T : T.name.wrapper)

function Base.show(io::IO, curve::T) where {T <: FinanceModels.Yield.AbstractYieldModel}
    r = zero(curve, 1)
    ylabel = isa(r.compounding, FinanceCore.Continuous) ? "Continuous" : "Periodic($(r.compounding.frequency))"
    kind = name(typeof(curve))
    l = UnicodePlots.lineplot(
        0.0, #from
        30.0,  # to
        t -> FinanceCore.rate(FinanceModels.zero(curve, t)),
        xlabel = "time",
        ylabel = ylabel,
        compact = true,
        name = "Zero rates",
        width = 60,
        title = "Yield Curve ($kind)"
    )
    return show(io, l)
end
end
