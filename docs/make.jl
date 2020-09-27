using Yields
using Documenter

makedocs(;
    modules=[Yields],
    authors="Alec Loudenback <alecloudenback@gmail.com> and contributors",
    repo="https://github.com/alecloudenback/Yields.jl/blob/{commit}{path}#L{line}",
    sitename="Yields.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alecloudenback.github.io/Yields.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/alecloudenback/Yields.jl",
)
