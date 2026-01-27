using SnoopCompileCore

invs = @snoop_invalidations begin
    using FinanceModels
end

using SnoopCompile
trees = invalidation_trees(invs)
println("Found $(length(trees)) invalidation trees\n")

# Look for any FinanceModels or FinanceCore specific invalidations
fc_trees = filter(trees) do tree
    s = string(tree)
    occursin("FinanceModels", s) || occursin("FinanceCore", s)
end

if isempty(fc_trees)
    println("SUCCESS: No FinanceModels/FinanceCore invalidation causes found!")
else
    println("Found $(length(fc_trees)) FinanceModels/FinanceCore related trees:")
    for tree in fc_trees
        println("\n", tree)
    end
end
