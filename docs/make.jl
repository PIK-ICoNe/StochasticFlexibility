using Pkg
if !(Pkg.project().path == joinpath(@__DIR__, "Project.toml"))
    Pkg.activate(@__DIR__)
end

# headless GK to fix ci
ENV["GKSwstype"] = "100"

using Documenter
using Literate



# precompile stuff now so the output won't show up in the docs
# using StochasticPrograms


# generate examples
examples = [
    joinpath(@__DIR__, "..", "examples", "flexibility_optimization.jl"),
]
OUTPUT = joinpath(@__DIR__, "src/generated")
isdir(OUTPUT) && rm(OUTPUT, recursive=true)
mkpath(OUTPUT)

for ex in examples
    Literate.markdown(ex, OUTPUT; documenter = true)
    Literate.script(ex, OUTPUT)
end

makedocs(; sitename = "Stochastic Flexibility Optimization",
    pages=[
        "Home" => "index.md",
        "Examples" => ["Simple Example" => "generated/flexibility_optimization.md"]
    ],
)

deploydocs(;
    repo="github.com/PIK-ICoNe/StochasticFlexibility",
    devbranch="main",
    push_preview=true,
)