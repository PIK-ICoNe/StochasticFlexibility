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
#-
F_range = 0.:2500.:25000.

costs = zeros(length(F_range))
inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
invs = Dict(([var => zeros(length(F_range)) for var in inv_vars]))
for i in eachindex(F_range)
    opt_data =  JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_$(F_range[i]).json"))
    costs[i] = opt_data["cost"]
    for var in inv_vars
        invs[var][i] = opt_data["inv"][string(var)]
    end
    opt_data = nothing
end

plot(F_range, costs, label = "cost")
savefig(joinpath(basepath, "paper_plots", "baseline_costs_F.png"))
for var in inv_vars
    plot(F_range, invs[var], label = string(var))
    savefig(joinpath(basepath, "paper_plots", "baseline_$(var)_F.png"))
end
#-

# also plot investments(F)