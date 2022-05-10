
basepath = realpath(joinpath((@__DIR__, "..")))

using Pkg
Pkg.activate(basepath)
Pkg.instantiate()

##

# headless GK to fix ci
ENV["GKSwstype"] = "100"

using Documenter
using Literate



# precompile stuff now so the output won't show up in the docs
# using StochasticPrograms


# generate examples
examples = [
    joinpath(@__DIR__, "..", "manuscripts", "flexibility_by_stochastic_programing.jl"),
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
        "Flexibility analysis using Stochastic Programing" => "index.md",
        "Manuscripts" => ["Main Story" => "generated/flexibility_by_stochastic_programing.md"]
    ],
)

deploydocs(;
    repo="github.com/PIK-ICoNe/StochasticFlexibility",
    devbranch="main",
    push_preview=true,
)