#=
Here we optimize the system with different levels of guaranteed flexibility, 
including no flexibility, in absence of flexibility requests.
=#

## Everything runs in the Project environment on the basepath
basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))

costs = []
F_range = 0.:5000.:25000.
for F in F_range
    c = CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg_$F.csv"), DataFrame)[!,1]
    append!(costs, c)
end

plot(F_range, costs)
savefig(joinpath(basepath, paper_plots, "baseline_costs_F.jl"))
#-
n_samples = [collect(15:5:35); collect(40:10:100)]
run_id = "naive_flex/2023-01-18"
# select F for further exploration
for (i, F) in enumerate(F_range)
    costs_flex = []
    append!(costs_flex, costs[i])
    for ns in n_samples
        c = load_costs(run_id, "$ns_F")
        append!(costs_flex, c)
    end
end
plot([[0.];n_samples], costs_flex)

#-
invs = []
for ns in n_samples
    inv_data = CSV.read(joinpath(basepath, results, run_id, "inv$(ns)_5000.0.csv"), DataFrame)
    inv = Dict(pairs(eachcol(inv_data)))
    inv_data = nothing;
    append!(invs, [inv])
end