using InterestModels
using Documenter

makedocs(;
    modules=[InterestModels],
    authors="Alec Loudenback <alecloudenback@gmail.com> and contributors",
    repo="https://github.com/alecloudenback/InterestModels.jl/blob/{commit}{path}#L{line}",
    sitename="InterestModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alecloudenback.github.io/InterestModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/alecloudenback/InterestModels.jl",
)
