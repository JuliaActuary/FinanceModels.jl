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
        "Home" => "index.md",
        "Guide" => [
            "Introduction" => "introduction.md",
            "Models, Valuation, and Fitting" => "models.md",
            "Contracts" => "contracts.md",
            "Rates" => "Rates.md",
            "Migration Guide" => "migration.md",
        ],
        "Modules" => [
            "FinanceModels" => "API/FinanceModels.md",
            "FinanceCore" => "API/FinanceCore.md",
            "Spline" => "API/Spline.md",
            "Fit" => "API/Fit.md",
            "Yield" => "API/Yield.md",
            "Bond" => "API/Bond.md",
            "Equity" => "API/Equity.md",
            "Option" => "API/Option.md",
            "Volatility" => "API/Volatility.md",
        ],
        "FAQs" => "faq.md",
    ],
    warnonly=true
)

deploydocs(;
    repo="github.com/JuliaActuary/FinanceModels.jl"
)
