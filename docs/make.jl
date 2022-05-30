
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

##

# generate examples
manuscript_path = realpath(joinpath(@__DIR__, "..", "manuscripts"))
manuscripts = [
    joinpath(manuscript_path, "flexibility_by_stochastic_programing.jl"),
]
OUTPUT = joinpath(@__DIR__, "src", "generated")
OUTPUT_SRCS = joinpath(@__DIR__, "src", "generated", "srcs")
isdir(OUTPUT) && rm(OUTPUT, recursive=true)
mkpath(OUTPUT)

##

src_path = realpath(joinpath(@__DIR__, "..", "src"))
srcs = [
    joinpath(src_path, "sp_model.jl"),
    joinpath(src_path, "evaluation_utils.jl"),
    joinpath(src_path, "plot_utils.jl"),
]

##

function set_manuscript_dir(content)
    content = replace(content, "@__DIR__" => "\"$manuscript_path\"")
    return content
end

for m in manuscripts
    Literate.markdown(m, OUTPUT; documenter = true, preprocess=set_manuscript_dir)
    Literate.script(m, OUTPUT; preprocess=set_manuscript_dir)
end

for s in srcs
    Literate.markdown(s, OUTPUT_SRCS; documenter = false)
end


##

makedocs(;
    sitename = "Stochastic Flexibility Optimization",
    pages=[
        "Flexibility analysis using Stochastic Programing" => "index.md",
        "Manuscripts" => ["Main Story" => "generated/flexibility_by_stochastic_programing.md"],
        "Sources" => [
            "Model" => "generated/srcs/sp_model.md",
            "Evaluation" => "generated/srcs/evaluation_utils.md",
            "Plotting" => "generated/srcs/plot_utils.md",
        ],
    ],

)

##

deploydocs(;
    repo="github.com/PIK-ICoNe/StochasticFlexibility",
    devbranch="main",
    push_preview=true,
)