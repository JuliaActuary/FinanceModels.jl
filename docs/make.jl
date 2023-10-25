using FinanceModels
using Documenter

makedocs(;
    modules=[FinanceModels, FinanceCore],
    authors="Alec Loudenback <alecloudenback@gmail.com> and contributors",
    repo="https://github.com/JuliaActuary/FinanceModels.jl/blob/{commit}{path}#L{line}",
    sitename="FinanceModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaActuary.github.io/FinanceModels.jl",
        assets=String[]
    ),
    pages=[
    # "Home" => "index.md",
    # "Guide" => "guide.md",
    # "Contracts" => "contracts.md",
    # "Models, Valuation, and Fitting" => "models.md",
    # "API Reference" => "api.md",
    # "Migration Guide" => "migration.md",
    # "FAQs" => "faq.md",
    ],
    warnonly=true
)

deploydocs(;
    repo="github.com/JuliaActuary/FinanceModels.jl"
)
