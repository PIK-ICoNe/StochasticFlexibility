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
using BSON
using JSON
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
#-
scen_freq = 96
run_id = "industry_district_heating"#"no_renewables"#"run_07_28"#"run_05_12"
n_samples = [collect(15:5:35); collect(40:10:60); collect(100:25:120)]
F_range = 5000.:5000.:25000.

data_matrix = Dict((:total_cost=>zeros(length(n_samples), length(F_range)), 
:CF => zeros(length(n_samples), length(F_range)), 
:CO => zeros(length(n_samples), length(F_range)),
:CI => zeros(length(n_samples), length(F_range))))
inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
for var in inv_vars
    push!(data_matrix, var => zeros(length(n_samples), length(F_range)))
end
for i in eachindex(n_samples)
    for j in eachindex(F_range)
        F = F_range[j]
        read_data = BSON.load(joinpath(basepath, "results", run_id, "conv_run_$(F)_$(scen_freq)_$(n_samples[i]).bson"))
        #cost_first_stage = get_total_investment(read_data) + get_operation_cost(read_data)
        data_matrix[:CO][i,j] = get_operation_cost(read_data)
        #CG = cost_first_stage - CB
        data_matrix[:CF][i,j] = read_data[:cost]-CB
        data_matrix[:total_cost][i,j] = read_data[:cost]
        for inv_var in inv_vars # convert units of components into euro
            inv = read_data[:inv][inv_var]*read_data[:params][Symbol(replace(string(inv_var), "u_"=>"c_"))]
            data_matrix[inv_var][i,j] = inv
        end
        data_matrix[:CI][i,j] = get_total_investment(read_data)
        read_data = nothing;
    end
end
#-
using CairoMakie
set_theme!(Theme(fontsize=45))
hm_vars = [:CF :CI; :u_storage :u_heat_storage]
hm_labels = ["CF" "CI"; "Cost of electricity storage" "Cost of heat storage"]

fig = Figure(resolution = (2400, 1800))
for i in 1:2
    for j in 1:2
        j_m = j
        (j == 2) && (j_m+=1) # this is done to place colorbars correctly
        # having one colorbar for all plots does not work as the values are too different
        ax = Axis(fig[i,j_m], xticklabelsvisible = (i==2),
        yticklabelsvisible = (j==1), xtickformat = "{:0d}",
        xlabel = L"F_G,\, MW", ylabel = L"n_m",
        xlabelvisible = (i==2), ylabelvisible = (j==1),
        title = hm_labels[i,j])
        hm = heatmap!(ax, F_range./1000, n_samples, data_matrix[hm_vars[i,j]]', 
        colormap=:blues)
        Colorbar(fig[i,j_m+1], hm)
    end
end
display(fig)
save(joinpath(basepath,"paper_plots", "convergence", "conv_full_$(run_id).png"), fig)