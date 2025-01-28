using PhysiCellECMCreator
using Documenter

DocMeta.setdocmeta!(PhysiCellECMCreator, :DocTestSetup, :(using PhysiCellECMCreator); recursive=true)

makedocs(;
    modules=[PhysiCellECMCreator],
    authors="Daniel Bergman <danielrbergman@gmail.com> and contributors",
    sitename="PhysiCellECMCreator.jl",
    format=Documenter.HTML(;
        canonical="https://drbergman.github.io/PhysiCellECMCreator.jl",
        edit_link=main,
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/drbergman/PhysiCellECMCreator.jl",
    devbranch="development",
)
