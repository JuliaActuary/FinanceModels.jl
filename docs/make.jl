using FinanceModels
using Documenter

makedocs(;
    modules=[FinanceModels],
    authors="Alec Loudenback <alecloudenback@gmail.com> and contributors",
    repo="https://github.com/JuliaActuary/FinanceModels.jl/blob/{commit}{path}#L{line}",
    sitename="FinanceModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaActuary.github.io/FinanceModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Developer Notes" => "developer.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaActuary/FinanceModels.jl",
)
