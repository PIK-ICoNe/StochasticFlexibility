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

n_samples = 20
params = [(5000., 48), (7500., 144), (25000., 240)]
run_id = "debug"#"cross_check
opt_modes = Dict((:det_inv_sto_flex => run_id*"/fixed_inv",
    #:sto_inv_det_flex => run_id*"/fixed_op",
    :det_inv_det_flex => "baseline",
    :sto_inv_sto_flex => run_id))#"run_02_28"))
inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
#-
full_cross_check = Dict((:cost => Dict(([mode => zeros(3) for mode in keys(opt_modes)])),
    :flex_cost => Dict(([mode => zeros(3) for mode in keys(opt_modes)])),
    :total_inv => Dict(([mode => zeros(3) for mode in keys(opt_modes)]))))
[push!(full_cross_check, var => Dict(([mode => zeros(3) for mode in keys(opt_modes)]))) for var in inv_vars]
#-
baseline_data = JSON.parsefile(joinpath(basepath, "results", "baseline", "baseline_0.0.json"))
baseline_cost = baseline_data["cost"]
baseline_total_inv = get_total_investment(baseline_data)
baseline_invs = Dict(([inv_var => baseline_data["inv"][string(inv_var)] for inv_var in inv_vars]))
opt_params = baseline_data["params"]
baseline_data = nothing
#-
for mode in keys(opt_modes)
    path = opt_modes[mode]
    for (i, (F,s)) in enumerate(params)
        if occursin("fixed", path)
            read_data = JSON.parsefile(joinpath(basepath, "results", path, "run_$(n_samples)_$(s)_nothing_nothing.json"))
        elseif path == "baseline" 
            read_data = JSON.parsefile(joinpath(basepath, "results", path, "baseline_$(F).json"))
        else
            read_data = JSON.parsefile(joinpath(basepath, "results", path, "run_$(n_samples)_$(s)_nothing_nothing.json"))#_$(F)_$(-F).json"))
        end
        full_cross_check[:cost][mode][i] = read_data["cost"]
        if path != "baseline"
            event_per_scen = 365*24 / (s + 1)
            c = (read_data["cost"] - baseline_cost)/event_per_scen
            full_cross_check[:flex_cost][mode][i] = c
        end
        for inv_var in inv_vars
            inv = read_data["inv"][string(inv_var)]
            # convert investment decisions into euro amount instead of kWp
            full_cross_check[inv_var][mode][i] = inv*opt_params[replace(string(inv_var), "u_" => "c_")]
        end
        full_cross_check[:total_inv][mode][i] = get_total_investment(read_data)
        read_data = nothing;
    end
end
#-
plotpath = joinpath(basepath, "results", run_id)
for var in vcat(inv_vars, :total_inv, :cost)
    plt = plot();
    for mode in (:det_inv_sto_flex, :sto_inv_sto_flex, :det_inv_det_flex)
        scatter!(plt, [string(p) for p in params], full_cross_check[var][mode],  label = string(mode), title = string(var), ylims = (0., 1.1*maximum(full_cross_check[var][mode])))
    end
    savefig(plt, joinpath(plotpath, "$var.png"))
end

plt = plot();
for mode in (:det_inv_sto_flex,  :sto_inv_sto_flex)
    scatter!(plt, [string(p) for p in params], full_cross_check[:flex_cost][mode],  label = string(mode), title = "flex_cost")
end
savefig(plt, joinpath(plotpath, "flex_cost.png"))

#-
using StatsPlots
using Cairo
using Gadfly
#-
for mode in [:sto_inv_sto_flex, :det_inv_det_flex]
    inv_plt = plot()

    full_investment = zeros(3, 6)
    full_investment[:, 1] = full_cross_check[:total_inv][mode]
    for i in eachindex(inv_vars)
        full_investment[:,(i+1)] = full_cross_check[inv_vars[i]][:sto_inv_sto_flex]
    end
    groupedbar!(inv_plt, full_investment, bar_position = :stack,label=hcat("total inv", reshape(string.(inv_vars), (1,5))))
    savefig(inv_plt, joinpath(plotpath, "$mode.png"))
end

#-
full_investment = zeros(3, 5)
#full_investment[:, 1] = full_cross_check[:total_inv][:det_inv_det_flex]
for i in eachindex(inv_vars)
    full_investment[:,i] = full_cross_check[inv_vars[i]][:det_inv_det_flex]
end

m = (:sto_inv_sto_flex, :det_inv_det_flex)
m_label = Dict((:sto_inv_sto_flex => "Stochastic", :det_inv_det_flex => "Deterministic"))
df = DataFrame()
for mode in eachindex(m)
    for var in inv_vars
        for i in 1:3
            append!(df, DataFrame(mode = mode, var = string(var), p = i, value = full_cross_check[var][m[mode]][i]))
        end
    end
end

draw(PNG(joinpath(plotpath, "full_inv.png"), 3inch, 3inch), @df df Gadfly.plot(x=:mode, y=:value, color=:var, xgroup=:p,
Geom.subplot_grid(Geom.bar(position=:stack)),
#Scale.color_discrete_manual(pal...), 
Scale.x_discrete(labels = i-> m_label[m[i]]),
#Scale.y_discrete(labels = i->string(params[i])),
Guide.xlabel("Investment mode"),
Guide.ylabel("Investment decision")))
@df df Gadfly.plot(x=:mode, y=:value, color=:var, xgroup=:p,
    Geom.subplot_grid(Geom.bar(position=:stack)),
    #Scale.color_discrete_manual(pal...), 
Scale.x_discrete(labels = i-> m_label[m[i]]),
    #Scale.y_discrete(labels = i->string(inv_vars[i])),
    Guide.xlabel("Investment mode"),
    Guide.ylabel("Investment decision"))
