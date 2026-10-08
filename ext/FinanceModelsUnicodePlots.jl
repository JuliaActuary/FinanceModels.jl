module FinanceModelsUnicodePlots
import FinanceModels
import UnicodePlots
import FinanceCore
# used to display simple type name in show method
# https://stackoverflow.com/questions/70043313/get-simple-name-of-type-in-julia?noredirect=1#comment123823820_70043313
name(::Type{T}) where {T} = (isempty(T.parameters) ? T : T.name.wrapper)

# 3-arg (MIME) show: rich plot display only where the display system asks for it
# (REPL results, notebooks). Overriding 2-arg `show` here would make every string
# interpolation, `print`, or log statement involving any yield model render a
# 60-character-wide plot.
function Base.show(io::IO, ::MIME"text/plain", curve::T) where {T <: FinanceModels.Yield.AbstractYieldModel}
    to = plot_end(curve)
    r = zero(curve, min(1, to))
    ylabel = isa(r.compounding, FinanceCore.Continuous) ? "Continuous" : "Periodic($(r.compounding.frequency))"
    kind = name(typeof(curve))
    l = UnicodePlots.lineplot(
        0.0, #from
        to,
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

# The plot ends at 30, or at a simulated path's last grid time, beyond which the path throws.
plot_end(curve) = 30.0
plot_end(p::FinanceModels.RatePath) = min(30.0, float(last(p.interp.t)))
end
