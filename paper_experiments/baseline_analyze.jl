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

# also plot investments(F)