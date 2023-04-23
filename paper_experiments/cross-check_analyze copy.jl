## Everything runs in the Project environment on the basepath

#=
In this experiment we compare flexibility potential of models with and without foreknowldge.
=#
basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using Dates
using Plots
using JSON

#= Include the file containing all the dependecies and optimization functions.
=#
#-
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));

timesteps = 1:24*365

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;

n_samples = 20
params = [(5000., 48), (7500., 144), (25000., 240)]
run_id = "cross_check_guar"#"cross_check
n_runs = 20
#-
check_df = DataFrame()
for n in 1:5
    for i in 2:3#eachindex(params)
        F, scen_freq = params[i]
        n_i = n + n_runs*(i - 1)
        read_data = JSON.parsefile(joinpath(basepath, "results", run_id, "run_$(n_i)_$(n_samples)_$(scen_freq)_$(F).json"))
        for inv_var in inv_vars
            inv = read_data["inv"][string(inv_var)]
            append!(check_df, DataFrame(var = string(inv_var), p = i, value = inv))
        end
    end
end

#-
for n in 1:n_runs
        F, scen_freq = params[1]

        read_data = JSON.parsefile(joinpath(basepath, "results", run_id, "run_$(n)_$(n_samples)_$(scen_freq)_$(F).json"))
        for inv_var in inv_vars
            inv = read_data["inv"][string(inv_var)]
            append!(check_df, DataFrame(var = string(inv_var), p = 1, value = inv))
        end
        read_data = nothing;
end

#-
baseline_data = JSON.parsefile(joinpath(basepath, "results", "baseline", "baseline_0.0.json"))
baseline_cost = baseline_data["cost"]
baseline_total_inv = get_total_investment(baseline_data)
baseline_invs = Dict(([inv_var => baseline_data["inv"][string(inv_var)] for inv_var in inv_vars]))
opt_params = baseline_data["params"]
baseline_data = nothing

#-
F, scen_freq = params[1]


F, scen_freq = params[1]
fixed_data = JSON.parsefile(joinpath(basepath, "results", run_id, "fixed_inv", "run_6_$(n_samples)_$(scen_freq)_$(F).json"))

event_per_scen = 365*24 / (scen_freq + 1)
c_fixed = (fixed_data["cost"] - baseline_cost)/event_per_scen
inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]

inv_dict = Dict(([var => Dict(([mode => 0. for mode in [:fixed, :sto, :sto_guar]])) for var in inv_vars]))
for inv_var in inv_vars
    inv = fixed_data["inv"][string(inv_var)]
    # convert investment decisions into euro amount instead of kWp
    inv_dict[inv_var][:fixed] = inv*opt_params[replace(string(inv_var), "u_" => "c_")]
end

fixed_data = nothing;
#-
sto_data = JSON.parsefile(joinpath(basepath, "results", run_id, "run_6_$(n_samples)_$(scen_freq)_$(F).json"))

event_per_scen = 365*24 / (scen_freq + 1)
c_sto = (sto_data["cost"] - baseline_cost)/event_per_scen

for inv_var in inv_vars
    inv = fixed_data["inv"][string(inv_var)]
    # convert investment decisions into euro amount instead of kWp
    inv_dict[inv_var][:sto] = inv*opt_params[replace(string(inv_var), "u_" => "c_")]
end

sto_data = nothing;
#-
guar_data = JSON.parsefile(joinpath(basepath, "results", run_id, "guar_flex_sto", "run_6_$(n_samples)_$(scen_freq)_$(F).json"))

event_per_scen = 365*24 / (scen_freq + 1)
c_sto_guar = (guar_data["cost"] - baseline_cost)/event_per_scen

for inv_var in inv_vars
    inv = guar_data["inv"][string(inv_var)]
    # convert investment decisions into euro amount instead of kWp
    inv_dict[inv_var][:sto_guar] = inv*opt_params[replace(string(inv_var), "u_" => "c_")]
end

sto_data = nothing;