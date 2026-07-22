using Documenter, IRKGaussLegendre

makedocs(
    sitename = "IRKGaussLegendre.jl",
    modules = [IRKGaussLegendre],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Public API" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/SciML/IRKGaussLegendre.jl.git"
)
