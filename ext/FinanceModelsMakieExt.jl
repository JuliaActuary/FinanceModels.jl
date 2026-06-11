module FinanceModelsMakieExt

using FinanceModels
import Makie

# adapted from
# https://github.com/JuliaMath/Polynomials.jl/blob/master/ext/PolynomialsMakieCoreExt.jl
# and
# https://docs.makie.org/stable/documentation/recipes/index.html#full_recipes_with_the_recipe_macro

function Makie.convert_arguments(P::Makie.PointBased, projection::Projection{C, M, K}) where {C, M, K <: CashflowProjection}
    p = collect(projection)
    amounts = [x.amount for x in p]
    times = [x.time for x in p]
    return Makie.convert_arguments(P, times, amounts)
end

function Makie.convert_arguments(P::Makie.PointBased, contract::A) where {A <: FinanceCore.AbstractContract}
    p = collect(contract)
    amounts = [x.amount for x in p]
    times = [x.time for x in p]
    return Makie.convert_arguments(P, times, amounts)
end

function Makie.convert_arguments(P::Makie.PointBased, cfs::Vector{C}) where {C <: FinanceCore.Cashflow}
    amounts = [x.amount for x in cfs]
    times = [x.time for x in cfs]
    return Makie.convert_arguments(P, times, amounts)
end


end
