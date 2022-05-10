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
    Literate.markdown(ex, OUTPUT)
    Literate.script(ex, OUTPUT)
end

# makedocs(;
#     modules=[BlockSystems],
#     authors="Hans Würfel <git@wuerfel.io> and contributors",
#     repo="https://github.com/hexaeder/BlockSystems.jl/blob/{commit}{path}#L{line}",
#     sitename="BlockSystems.jl",
#     format=Documenter.HTML(;
#         prettyurls=get(ENV, "CI", "false") == "true",
#         canonical="https://hexaeder.github.io/BlockSystems.jl",
#         assets=String[],
#     ),
#     pages=[
#         "Home" => "index.md",
#         "Examples" => ["Spacecraft" => "generated/spacecraft.md",
#                        # "PowerDynamics.jl Node" => "generated/pd_node.md",
#                        "Kuramoto Network" => "generated/kuramoto_network.md",
#                        "Kuramoto without ND.jl" => "generated/kuramoto_without_nd.md"]
#     ],
# )

# deploydocs(;
#     repo="github.com/hexaeder/BlockSystems.jl",
#     devbranch="main",
#     push_preview=true,
# )