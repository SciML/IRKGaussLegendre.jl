using Documenter, IRKGaussLegendre

makedocs(
    sitename = "IRKGaussLegendre.jl",
    modules = [IRKGaussLegendre],
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/SciML/IRKGaussLegendre.jl.git"
)
