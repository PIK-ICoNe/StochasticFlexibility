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

#pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
#pars[:inv_budget] = 10^10;

n_samples = 20
params = [(5000., 48), (7500., 144), (25000., 240)]
opt_modes = Dict((:det_inv_sto_flex => "cross_check/fixed_inv",
    :sto_inv_det_flex => "cross_check/fixed_op",
    :det_inv_det_flex => "baseline",
    :sto_inv_sto_flex => "run_02_28"))
inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
#-
full_cross_check = Dict(([inv => merge(Dict(([k => 0. for k in keys(opt_modes)])), Dict((:baseline => 0.))) for inv in inv_vars]))
push!(full_cross_check, :cost => merge(Dict(([k => 0. for k in keys(opt_modes)])), Dict((:baseline => 0.))))
push!(full_cross_check, :total_inv => merge(Dict(([k => 0. for k in keys(opt_modes)])), Dict((:baseline => 0.))))
push!(full_cross_check, :flex_cost => merge(Dict(([k => 0. for k in keys(opt_modes)])), Dict((:baseline => 0.))))

full_cross_check = [merge(full_cross_check, Dict((:scen_freq => p[2], :F_guar => p[1]))) for p in params]
#-
baseline_data = JSON.parsefile(joinpath(basepath, "results", "baseline", "baseline_0.0.json"))
baseline_cost = baseline_data["cost"]
total_inv = get_total_investment(baseline_data)
for i in 1:3
    full_cross_check[i][:cost][:baseline] = baseline_cost
    full_cross_check[i][:total_inv][:baseline] = total_inv
    for inv_var in inv_vars
        full_cross_check[i][inv_var][:baseline] = baseline_data["inv"][string(inv_var)]
    end
end
baseline_data = nothing
#-
for mode in keys(opt_modes)
    path = opt_modes[mode]
    for (i, (F,s)) in enumerate(params)
        #(F, s) = params[i]
        println(i)
        println(F)
        println(s)
        if occursin("fixed", path)
            read_data = JSON.parsefile(joinpath(basepath, "results", path, "run_$(n_samples)_$(s)_nothing_nothing.json"))
            println("here")
            println("storage $(read_data["inv"]["u_storage"])")
        elseif path == "baseline" 
            read_data = JSON.parsefile(joinpath(basepath, "results", path, "baseline_$(F).json"))
        else
            read_data = JSON.parsefile(joinpath(basepath, "results", path, "run_$(n_samples)_$(s)_$(F)_$(-F).json"))
        end
        full_cross_check[i][:cost][mode] = read_data["cost"]
        if path != "baseline"
            event_per_scen = 365*24 / (s + 1)
            c = (read_data["cost"] - baseline_cost)/event_per_scen
            full_cross_check[i][:flex_cost][mode] = c
        end
        for inv_var in inv_vars
            inv = read_data["inv"][string(inv_var)]
            full_cross_check[i][inv_var][mode] = inv
            println("$mode, $inv_var = $inv")
        end
        full_cross_check[i][:total_inv][mode] = get_total_investment(read_data)
        read_data = nothing;
    end
end
#-
plt = plot()
for mode in keys(full_cross_check[1][:cost])
    scatter!(plt,[full_cross_check[i][:u_storage][mode] for i in 1:3])
end
plt
plotpath = joinpath(basepath, "results", "cross_check")

#=
    scatter!(plt, 1:3, [full_cross_check[i][var][:baseline] for i in 1:3], label = string(var))
end=#
for var in vcat(inv_vars, :total_inv)
    plt = plot();
    for mode in keys(full_cross_check[1][:cost])
        scatter!(plt, [string(p) for p in params], [full_cross_check[i][var][mode] for i in 1:3],  label = string(mode), title = string(var))
    end
    savefig(plt, joinpath(plotpath, "$var.png"))
end

for var in [:cost]
    plt = plot();
    for mode in (:det_inv_sto_flex, :sto_inv_sto_flex)
        scatter!(plt, [string(p) for p in params], [full_cross_check[i][var][mode] for i in 1:3],  label = string(mode), title = "Flexibility cost")
    end
    savefig(plt, joinpath(plotpath, "$var.png"))
end
#-
#=savefig(plot(data_matrix[:cost],
    title = "Cost of flexibility"), joinpath(plotpath, "costs.png"))
#-
savefig(plot(data_matrix[:u_wind],
    title = "Wind, kWp"), joinpath(plotpath, "u_wind.png"))

savefig(plot(data_matrix[:u_storage], title = "Storage, kWh"), joinpath(plotpath, "u_storage.png"))

savefig(plot(data_matrix[:u_heat_storage], title = "Storage, kWh"), joinpath(plotpath, "u_heatstorage.png"))
savefig(plot(data_matrix[:u_pv],  
    title = "PV, kWp"), joinpath(plotpath, "u_pv.png"))

savefig(plot(data_matrix[:u_heatpump], 
    title = "Wind, kWp"), joinpath(plotpath, "u_heatpump.png"))

savefig(plot(data_matrix[:total_inv], 
    title = "Total investment"), joinpath(plotpath, "total_inv.png"))
#-=#