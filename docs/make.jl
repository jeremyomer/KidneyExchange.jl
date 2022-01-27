using Documenter
using KidneyExchange

makedocs(
    sitename = "KidneyExchange",
    format = Documenter.HTML(),
    modules = [KidneyExchange]
)

deploydocs(;
    branch = "gh-pages",
    devbranch = "master",
    repo   = "github.com/jeremyomer/KidneyExchange.jl"
)
