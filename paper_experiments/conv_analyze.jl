## Everything runs in the Project environment on the basepath

#=
In this experiment we compare flexibility potential of models with and without foreknowldge.

First we load the data for the background model. We fix the 
=#
basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using Dates
using Plots

#= Include the file containing all the dependecies and optimization functions.
=#
#-
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));

#-
base_model = JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_0.0.json"))
CB = base_model["cost"]
CIB = get_total_investment(base_model)
COB = get_operation_cost(base_model)
inv_base = base_model["inv"]
P = base_model["params"]
base_model = nothing;
#=bkg_cost = []
F = 5000.:5000.:25000.
append!(bkg_cost, CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg.csv"), DataFrame)[!,1][1])
for f in F
    append!(bkg_cost, CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg_$f.csv"), DataFrame)[!,1][1])
end
=#
#-
scen_freq = 96
run_id = "run_05_12"
n_samples = [collect(15:5:35); collect(40:10:90); collect(100:25:120)]
F_range = 5000.:5000.:25000.

data_matrix = Dict((:cost=>zeros(length(n_samples), length(F_range))))
inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
for var in inv_vars
    push!(data_matrix, var => zeros(length(n_samples), length(F_range)))
end
push!(data_matrix, :total_inv => zeros(length(n_samples), length(F_range)))
event_per_scen = 365*24 / (96 + 1)
for i in eachindex(n_samples)
    for j in eachindex(F_range)
        F = F_range[j]
        read_data = BSON.load(joinpath(basepath, "results", run_id, "conv_run_$(F)_$(scen_freq)_$(n_samples[i]).bson"))
        cost_first_stage = get_total_investment(read_data) + get_operation_cost(read_data)
        c = cost_first_stage - CB
        #c = read_data["cost"]/event_per_scen
        data_matrix[:cost][i,j] = c
        for inv_var in inv_vars
            inv = read_data[:inv][inv_var]
            data_matrix[inv_var][i,j] = inv
        end
        data_matrix[:total_inv][i,j] = get_total_investment(read_data)
        read_data = nothing;
    end
end
#-
plotpath = joinpath(basepath, "paper_plots")
savefig(heatmap(F_range, n_samples, data_matrix[:cost],
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", 
    title = "Cost of guaranteeing flexibility (CG)",
    clim = (0., Inf)), 
    joinpath(plotpath, "conv_costs.png"))
#-
savefig(heatmap(F_range, n_samples, data_matrix[:u_wind], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", 
    title = "Wind, kWp",
    clim = (0., Inf)), joinpath(plotpath, "conv_u_wind.png"))

savefig(heatmap(F_range, n_samples,data_matrix[:u_storage], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", title = "Storage, kWh",
    clim = (0., Inf)), 
    joinpath(plotpath, "conv_u_storage.png"))

savefig(heatmap(F_range, n_samples,data_matrix[:u_heat_storage], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", title = "Heat storage, kWh",
    clim = (0., Inf)),
    joinpath(plotpath, "conv_u_heatstorage.png"))

savefig(heatmap(F_range, n_samples, data_matrix[:u_pv], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", 
    title = "PV, kWp",
    clim = (0., Inf)), 
    joinpath(plotpath, "conv_u_pv.png"))

savefig(heatmap(F_range, n_samples, data_matrix[:u_heatpump], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", 
    title = "Heat pump, kWp",
    clim = (0., Inf)), 
    joinpath(plotpath, "conv_u_heatpump.png"))

savefig(heatmap(F_range, n_samples, data_matrix[:total_inv], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", 
    title = "Total investment",
    clim = (0., Inf)), 
    joinpath(plotpath, "conv_total_inv.png"))
#-