using TKTDsimulations
using Documenter

makedocs(;
    modules=[TKTDsimulations],
    authors="Virgile Baudrot",
    repo="https://github.com/virgile-baudrot/TKTDsimulations.jl/blob/{commit}{path}#L{line}",
    sitename="TKTDsimulations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://virgile-baudrot.github.io/TKTDsimulations.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/virgile-baudrot/TKTDsimulations.jl",
)
