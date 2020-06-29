using YieldModels
using Documenter

makedocs(;
    modules=[YieldModels],
    authors="Alec Loudenback <alecloudenback@gmail.com> and contributors",
    repo="https://github.com/alecloudenback/YieldModels.jl/blob/{commit}{path}#L{line}",
    sitename="YieldModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alecloudenback.github.io/YieldModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/alecloudenback/YieldModels.jl",
)
