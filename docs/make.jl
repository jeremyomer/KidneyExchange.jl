using Documenter
using KidneyExchange
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), sorting = :nyt)

makedocs(
    bib,
    sitename = "KidneyExchange.jl",
    strict = true,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [KidneyExchange],
    pages = [
        "Home"       => "index.md",
        "Models"     => "models.md",
        "Functions"  => "functions.md",
        "Types"      => "types.md",
        "References" => "references.md"
    ]
)

deploydocs(;
    branch = "gh-pages",
    devbranch = "master",
    repo   = "github.com/jeremyomer/KidneyExchange.jl"
)
