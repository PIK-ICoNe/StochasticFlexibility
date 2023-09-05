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
scen_freq = 48:48:144
F_range = [250., 500., 1000., 2500., 5000.]
run_id = "new_parameters"
n_samples = 20

data_matrix = Dict((:total_cost=>zeros(length(scen_freq), length(F_range)), 
:CF => zeros(length(scen_freq), length(F_range)), 
:CF_norm => zeros(length(scen_freq), length(F_range)),
:total_flex => zeros(length(scen_freq), length(F_range)),
:CO => zeros(length(scen_freq), length(F_range)),
:CI => zeros(length(scen_freq), length(F_range))))
inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
for var in inv_vars
    push!(data_matrix, var => zeros(length(scen_freq), length(F_range)))
end
for i in eachindex(scen_freq)
    sf = scen_freq[i]
    event_per_scen = 365*24 / (sf + 1)
    for j in eachindex(F_range)
        F = F_range[j]
        read_data = BSON.load(joinpath(basepath, "results", run_id, "run_$(F)_$(scen_freq[i]).bson"))
        #cost_first_stage = get_total_investment(read_data) + get_operation_cost(read_data)
        data_matrix[:CO][i,j] = get_operation_cost(read_data)
        #CG = cost_first_stage - CB
        CF = read_data[:cost]-CB
        total_flex = (length(timesteps)-24) / (sf + 1)*F*0.4 
        norm_CF = CF/total_flex
        data_matrix[:CF][i,j] = CF
        data_matrix[:CF_norm][i,j] = norm_CF
        data_matrix[:total_flex][i,j] = total_flex

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
hm_vars = [:CF :CI; :u_storage :u_heat_storage; :CF_norm :total_flex]
hm_labels = ["CF" "CI"; "Cost of electricity storage" "Cost of heat storage"; "Normalized CF" "Total flexibility in a year"]

add_case_studies = false
fig = Figure(resolution = (2400, 1800))
for i in 1:3
    for j in 1:2
        j_m = j
        (j == 2) && (j_m+=1) # this is done to place colorbars correctly
        # having one colorbar for all plots does not work as the values are too different
        ax = Axis(fig[i,j_m], xticklabelsvisible = (i==3),
        yticklabelsvisible = (j==1), xtickformat = "{:0d}",
        xlabel = L"F_G,\, kW", ylabel = L"\Delta_t",
        xlabelvisible = (i==3), ylabelvisible = (j==1),
        title = hm_labels[i,j])
        hm = heatmap!(ax, F_range, scen_freq, data_matrix[hm_vars[i,j]]', 
        colormap=:blues)
        Colorbar(fig[i,j_m+1], hm)
        if add_case_studies # TODO update case studies
            scatter!(ax, [5,25,25], [48,48,240], color=:red, markersize=20)
            text!(ax, [5,25,25], [48,48,240], text=["small", "base", "rare"], align=(:center, :bottom), offset=(0,20), color=:gray)
        end
    end
end
display(fig)
case_studies = ""
(add_case_studies)&&(case_studies="with_markers")
savepath = joinpath(basepath, "paper_plots", run_id)
if !isdir(savepath)
    mkdir(savepath)
end
save(joinpath(savepath, "request_rate_full$(case_studies).png"), fig)