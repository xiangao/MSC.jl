using Pkg

ENV["GKSwstype"] = "100"

Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using MSC

DocMeta.setdocmeta!(MSC, :DocTestSetup, :(using MSC); recursive = true)

makedocs(
    sitename = "MSC.jl",
    modules = [MSC],
    format = Documenter.HTML(
        edit_link = nothing,
        repolink = "https://github.com/xiangao/MSC.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Vignettes" => [
            "Getting Started" => "vignettes/01_getting_started.md",
            "DataFrame Panels" => "vignettes/02_dataframe_panels.md",
            "Diagnostics and Placebos" => "vignettes/03_diagnostics_placebos.md",
        ],
        "Reference" => "reference.md",
    ],
    warnonly = true,
    checkdocs = :none,
    remotes = nothing,
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/xiangao/MSC.jl.git",
        devbranch = "main",
        push_preview = false,
    )
end
