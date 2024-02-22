using Documenter
using KidneyExchange
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), style = :authoryear)

makedocs(;
    plugins = [bib],
    sitename = "KidneyExchange.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [KidneyExchange],
    pages = [
        "Home" => "index.md",
        "Types" => "types.md",
        "Functions" => "functions.md",
        "References" => "references.md",
    ],
)

deploydocs(;
    branch = "gh-pages",
    devbranch = "master",
    repo = "github.com/jeremyomer/KidneyExchange.jl",
)
