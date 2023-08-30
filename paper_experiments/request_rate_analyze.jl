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
#=bkg_cost = []
F = 5000.:5000.:25000.
append!(bkg_cost, CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg.csv"), DataFrame)[!,1][1])
for f in F
    append!(bkg_cost, CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg_$f.csv"), DataFrame)[!,1][1])
end
=#
base_model = JSON.parsefile(joinpath(basepath, "results", "baseline", "baseline_0.0.json"))
CB = base_model["cost"]
CIB = get_total_investment(base_model)
COB = get_operation_cost(base_model)
inv_base = base_model["inv"]
P = base_model["params"]
base_model = nothing;
#-
scen_freq = 48:48:144#240
F_range = 5000.:5000.:25000.
run_id = "run_07_28"#"run_05_12"
n_samples = 20

data_matrix = Dict((:total_cost=>zeros(length(scen_freq), length(F_range)), 
:CF => zeros(length(scen_freq), length(F_range)), 
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
hm_vars = [:CF :CI; :u_storage :u_heat_storage; :u_heatpump :u_heatpump]
hm_labels = ["CF" "CI"; "Cost of electricity storage" "Cost of heat storage"; "Cost of heat pump" "Cost of heat pump"]

add_case_studies = true
fig = Figure(resolution = (2400, 1800))
for i in 1:2
    for j in 1:2
        j_m = j
        (j == 2) && (j_m+=1) # this is done to place colorbars correctly
        # having one colorbar for all plots does not work as the values are too different
        ax = Axis(fig[i,j_m], xticklabelsvisible = (i==2),
        yticklabelsvisible = (j==1), xtickformat = "{:0d}",
        xlabel = L"F_G,\, MW", ylabel = L"\Delta_t",
        xlabelvisible = (i==2), ylabelvisible = (j==1),
        title = hm_labels[i,j])
        hm = heatmap!(ax, F_range./1000, scen_freq, data_matrix[hm_vars[i,j]]', 
        colormap=:blues)
        Colorbar(fig[i,j_m+1], hm)
        if add_case_studies
            scatter!(ax, [5,25,25], [48,48,240], color=:red, markersize=20)
            text!(ax, [5,25,25], [48,48,240], text=["small", "base", "rare"], align=(:center, :bottom), offset=(0,20), color=:gray)
        end
    end
end
display(fig)
case_studies = ""
(add_case_studies)&&(case_studies="with_markers")
save(joinpath(basepath,"paper_plots", "request_rate", "request_rate_full_$(run_id)_$(case_studies).png"), fig)