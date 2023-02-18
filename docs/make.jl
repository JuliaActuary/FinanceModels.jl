using Yields
using Documenter

makedocs(;
    modules=[Yields],
    authors="Alec Loudenback <alecloudenback@gmail.com> and contributors",
    repo="https://github.com/JuliaActuary/Yields.jl/blob/{commit}{path}#L{line}",
    sitename="Yields.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaActuary.github.io/Yields.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Instruments and Quotes" => "Instruments.md",
        "API Reference" => "api.md",
        "Developer Notes" => "developer.md",
        "Upgrading from prior versions" => "Upgrading.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaActuary/Yields.jl",
)
