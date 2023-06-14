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
using BSON
#= Include the file containing all the dependecies and optimization functions.
=#
#-
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));

#-
#=bkg_cost = []
F = 5000.:5000.:25000.
append!(bkg_cost, CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg.csv"), DataFrame)[!,1][1])
for f in F
    append!(bkg_cost, CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg_$f.csv"), DataFrame)[!,1][1])
end
=#
base_model = JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_0.0.json"))
CB = base_model["cost"]
CIB = get_total_investment(base_model)
COB = get_operation_cost(base_model)
inv_base = base_model["inv"]
P = base_model["params"]
base_model = nothing;
#-
scen_freq = 48:48:240
F_range = 5000.:5000.:25000.
run_id = "run_05_12"
n_samples = 20

data_matrix = Dict((:cost=>zeros(length(scen_freq), length(F_range))))
inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
for var in inv_vars
    push!(data_matrix, var => zeros(length(scen_freq), length(F_range)))
end
push!(data_matrix, :total_inv => zeros(length(scen_freq), length(F_range)))

df = DataFrame()
for i in eachindex(scen_freq)
    sf = scen_freq[i]
    event_per_scen = 365*24 / (sf + 1)
    for j in eachindex(F_range)
        F = F_range[j]
        read_data = BSON.load(joinpath(basepath, "results", run_id, "run_$(F)_$(scen_freq[i]).bson"))
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
savefig(heatmap(F_range, scen_freq, data_matrix[:cost],
    xlabel = "Guaranteed flexibility",
    ylabel = "Expected offset between requests, hours", 
    clim = (0., Inf),
    title = "Cost of guaranteeing flexibility (CG)"), joinpath(plotpath, "costs.png"))
#-
savefig(heatmap(F_range, scen_freq, data_matrix[:u_wind], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Expected offset between requests, hours",
    clim = (0., Inf),
    title = "Wind, kWp"), joinpath(plotpath, "u_wind.png"))

savefig(heatmap(F_range, scen_freq,data_matrix[:u_storage], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Expected offset between requests, hours",
    clim = (0., Inf),
    title = "Storage, kWh"), joinpath(plotpath, "u_storage.png"))

savefig(heatmap(F_range, scen_freq,data_matrix[:u_heat_storage], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Expected offset between requests, hours",
    clim = (0., Inf), 
    title = "Heat storage, kWh"), 
    joinpath(plotpath, "u_heatstorage.png"))

savefig(heatmap(F_range, scen_freq, data_matrix[:u_pv], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Expected offset between requests, hours",
    clim = (0., Inf),
    title = "PV, kWp"), joinpath(plotpath, "u_pv.png"))

savefig(heatmap(F_range, scen_freq, data_matrix[:u_heatpump], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Expected offset between requests, hours",
    clim = (0., Inf),
    title = "Heat pump, kWp"), joinpath(plotpath, "u_heatpump.png"))

savefig(heatmap(F_range, scen_freq, data_matrix[:total_inv], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Expected offset between requests, hours",
    clim = (0., Inf),
    title = "Total investment"), 
    joinpath(plotpath, "total_inv.png"))
#-