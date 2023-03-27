using MOTIFs
using Documenter

DocMeta.setdocmeta!(MOTIFs, :DocTestSetup, :(using MOTIFs); recursive=true)

makedocs(;
    modules=[MOTIFs],
    authors="Shane Kuei-Hsien Chu (skchu@wustl.edu)",
    repo="https://github.com/kchu25/MOTIFs.jl/blob/{commit}{path}#{line}",
    sitename="MOTIFs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kchu25.github.io/MOTIFs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/MOTIFs.jl",
    devbranch="main",
)
