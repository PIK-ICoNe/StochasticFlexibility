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
include(joinpath(basepath, "paper_experiments", "stochastic_optimization_setup.jl"));

#-
run_id = "new_parameters"
base_model = BSON.load(joinpath(basepath, "results", run_id, "baseline/baseline_0.0.bson"))
CB = base_model[:cost]
CIB = get_total_investment(base_model)
COB = get_operation_cost(base_model)
inv_base = base_model[:inv]
P = base_model[:params]
base_model = nothing;
#-
scen_freq = 96
n_samples = [collect(20:10:70); collect(80:20:120)]
F_range = [500.]
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
        ax = Axis(fig[i,j], xticklabelsvisible = (i==2),
        yticklabelsvisible = true, xtickformat = "{:0d}",
        xlabel = L"F_G,\, kW", ylabel = L"n_m",
        xlabelvisible = (i==2), ylabelvisible = (j==1),
        title = hm_labels[i,j])
        lines!(n_samples, vec(data_matrix[hm_vars[i,j]]))
    end
end
display(fig)
savepath = joinpath(basepath, "paper_plots", run_id)
if !isdir(savepath)
    mkdir(savepath)
end
save(joinpath(savepath, "conv_full.png"), fig)