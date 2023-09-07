module FinanceModelsMakieCoreExt

using FinanceModels
import MakieCore

# adapted from
# https://github.com/JuliaMath/Polynomials.jl/blob/master/ext/PolynomialsMakieCoreExt.jl
# and
# https://docs.makie.org/stable/documentation/recipes/index.html#full_recipes_with_the_recipe_macro

function MakieCore.convert_arguments(P::MakieCore.PointBased, projection::Projection{C,M,K}) where {C,M,K<:CashflowProjection}
    p = collect(projection)
    amounts = [x.amount for x in p]
    times = [x.time for x in p]
    return MakieCore.convert_arguments(P, times, amounts)
end

function MakieCore.convert_arguments(P::MakieCore.PointBased, contract::A) where {A<:FinanceCore.AbstractContract}
    p = collect(contract)
    amounts = [x.amount for x in p]
    times = [x.time for x in p]
    return MakieCore.convert_arguments(P, times, amounts)
end

function MakieCore.convert_arguments(P::MakieCore.PointBased, cfs::Vector{C}) where {C<:FinanceCore.Cashflow}
    amounts = [x.amount for x in cfs]
    times = [x.time for x in cfs]
    return MakieCore.convert_arguments(P, times, amounts)
end


end