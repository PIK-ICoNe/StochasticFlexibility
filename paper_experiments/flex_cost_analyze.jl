#=
We focus on three particular pairs of F and scen_freq:
(5000., 48), (7500., 144), (25000., 240)
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()
using DataFrames
using CSV
using JSON
using BSON
#-
include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

#-
n_samples = 50
params = [(5000., 48), (7500., 144), (25000., 240)]
n_runs = 3
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;
savepath = joinpath(basepath, "results/flex_cost_26_04")
filetype = "bson"
#-
if isfile(joinpath(savepath, "full_cost_data.csv")) && isfile(joinpath(savepath, "full_inv_data.csv"))
    df_costs = DataFrame(CSV.File(joinpath(savepath, "full_cost_data.csv")))
    df_inv = DataFrame(CSV.File(joinpath(savepath, "full_inv_data.csv")))
else
    df_costs = DataFrame()
    df_inv = DataFrame()
    for i in eachindex(params)
        F, sf = params[i]
        for n in 1:n_runs
            for m in ("fixed_fg", "fixed_fg_inv", "unfixed_guar_flex")
                p = joinpath(savepath, "run_$(F)_$(sf)_$(n)_$(m)."*filetype)
                if isfile(p)
                    if filetype == "json"
                        opt = JSON.parsefile(p)
                        total_cost = opt["cost"]
                    elseif filetype == "bson"
                        opt = BSON.load(p)
                        total_cost = opt[:cost]
                    end
                    CI = get_total_investment(opt)
                    CO = get_operation_cost(opt)
                    CR = total_cost - CI - CO
                    append!(df_costs, DataFrame(p = i, total_cost = total_cost,
                        CI = CI, CO = CO, CR = CR, sample_num = n, mode = m))
                    if m !== "fixed_inv_fg"
                        if filetype == "json"
                            append!(df_inv, DataFrame(Dict((vcat([var => opt["inv"][var] for var in string.(inv_vars)], 
                            "mode" => m, "sample_num" => n, "p" => i)))))
                        elseif filetype == "bson"
                            append!(df_inv, DataFrame(Dict((vcat([string(var) => opt[:inv][var] for var in inv_vars], 
                            "mode" => m, "sample_num" => n, "p" => i)))))
                        end
                    end
                    opt = nothing;
                else
                    println("File $p does not exist")
                end
            end
        end
    end
    CSV.write(joinpath(savepath, "full_cost_data.csv"), df_costs)
    CSV.write(joinpath(savepath, "full_inv_data.csv"), df_inv)
end

#-
# Fist we get C^B
base_model = JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_0.0.json"))
CB = base_model["cost"]
CIB = get_total_investment(base_model)
COB = get_operation_cost(base_model) #TODO check why COB+CIB!=CB
inv_base = base_model["inv"]
P = base_model["params"]
base_model = nothing;
#-
# Plot investments for all components
df = DataFrame()
cost_vars = []#[:CO]
modes = ("fixed_fg", "fixed_fg_inv", "unfixed_guar_flex")
for p in eachindex(params)
    for (m_i, mode) in enumerate(modes)
        for (var_index, cost_var) in enumerate(cost_vars)
            append!(df, DataFrame(p = p, mode = m_i, var = var_index,
                mean_value = mean(subset(df_costs, :p => a -> a .== p, :mode => m -> m .== mode)[!,cost_vars[var_index]])))
        end
        for (var_index, inv_var) in enumerate(inv_vars)
            append!(df, DataFrame(p = p, mode = m_i, var = var_index+length(cost_vars),
                mean_value = mean(subset(df_inv, :p => a -> a .== p, 
                :mode => m -> m .== mode)[!,inv_vars[var_index]])*P[replace(string(inv_var), "u_" => "c_")]))
        end
    end
end
# Add base investments for comparison
for (var_index, inv_var) in enumerate(inv_vars)
    append!(df, DataFrame(p = 0, mode = 2, var = var_index+length(cost_vars),
        mean_value = inv_base[string(inv_var)]*P[replace(string(inv_var), "u_" => "c_")]))
end

#-

using CairoMakie
using LaTeXStrings
colors = Makie.wong_colors()
fig_inv = Figure()
ax = Axis(fig_inv[1,1], yticks = (0:3, vcat("Base model", string.(params))),
        title = "Investment decision")
labels = string.(vcat(cost_vars, inv_vars));
elements = [PolyElement(polycolor = colors[i]); for i in 1:length(labels)]
title = "Components";
Legend(fig_inv[1,2], elements, labels, title)
barlabels = ["" for i in 1:size(df)[1]];
for i in 1:size(df)[1]
    b = ["Fixed FG", "Fixed FG Inv", "OF"]#[L"OF_{|\overline{FG}}", L"OF_{|\overline{Inv}_{fg}}",L"OF"]
    if df[i, :var] == 5
        if df[i, :p] != 0
            barlabels[i] = b[df[i, :mode]]
        end
    end
end
barplot!(ax, df[!, :p], df[!, :mean_value],
dodge = df[!, :mode],
direction=:x,
flip_labels_at=0.85,
color = colors[df[!, :var]],stack = df[!, :var], bar_labels = barlabels)
fig_inv
#-
save(joinpath(basepath, "paper_plots/opt_mode_investments.png"), fig_inv);
fig_inv
#-
df_bar = DataFrame()
cost_vars = [:CR, :CO, :CI]
modes = ("fixed_fg", "fixed_fg_inv", "unfixed_guar_flex")
for p in eachindex(params)
    for (m_i, mode) in enumerate(modes)
        for (var_index, cost_var) in enumerate(cost_vars)
            append!(df_bar, DataFrame(p = p, mode = m_i*2, var = var_index, width = 1.,gap=.1,
                mean_value = mean(subset(df_costs, :p => a -> a .== p, :mode => m -> m .== mode)[!,cost_vars[var_index]])))
        end
        #total cost
        append!(df_bar, DataFrame(p=p, mode = m_i*2-1, var = 4, width = 0.9,gap=.2,
        mean_value = mean(subset(df_costs, :p => a -> a .== p, :mode => m -> m .== mode)[!,:total_cost])))
    end
end
append!(df_bar, DataFrame(p=0, mode = 4, var = 2, width = 1.,gap=.1, mean_value=COB))
append!(df_bar, DataFrame(p=0, mode = 4, var = 3, width = 1.,gap=.1, mean_value=CIB))
append!(df_bar, DataFrame(p=0, mode = 3, var = 4, width = 0.9,gap=.2, mean_value=CB))

fig_cost = Figure()
ax = Axis(fig_cost[1,1], yticks = (0:3, vcat("Base model",string.(params))),
        title = "Flexibility cost")

labels = string.(vcat(cost_vars, ["total cost"]));
elements = [PolyElement(polycolor = colors[i]); for i in 1:length(labels)]
title = "Components";
Legend(fig_cost[1,2], elements, labels, title)
barlabels = ["" for i in 1:size(df_bar)[1]];
for i in 1:size(df_bar)[1]
    b = ["Fixed FG", "Fixed FG Inv", "OF","Total cost"]#[L"OF_{|\overline{FG}}", L"OF_{|\overline{Inv}_{fg}}",L"OF"]
    if df_bar[i, :var] == length(labels)-1
        barlabels[i] = b[Int(df_bar[i, :mode]/2)]
    end
    if df_bar[i, :var] == length(labels)
        barlabels[i] == "Total cost"
    end
end
barplot!(ax, df_bar[!, :p], df_bar[!, :mean_value],
dodge = df_bar[!, :mode],
direction=:x,
#width = df_bar[!, :width],
flip_labels_at=0.95,
color = colors[df_bar[!, :var]],stack = df_bar[!, :var], bar_labels = barlabels)
fig_cost
#-
save(joinpath(basepath, "paper_plots/Flexibility_cost_with_base.png"), fig_cost)

#-
# Stack area chart for flex potential
using VegaLite
include(joinpath(basepath, "src", "evaluation_utils.jl"))

#test_model = JSON.parsefile(joinpath(savepath, "run_25000.0_240_17_unfixed_guar_flex.json"))
test_model = BSON.load(joinpath(savepath, "run_25000.0_240_17_unfixed_guar_flex.bson"))
F_pos, F_neg, F_pos_d, F_neg_d = naive_flex_potential(test_model, pv, wind, timesteps)
flex_pos = DataFrame(t = timesteps[1:100], cur = F_pos_d[1:100,1], sto = F_pos_d[1:100,2], heat = F_pos_d[1:100,3])
flex_neg = DataFrame(t = timesteps[1:100],neg_cur = F_neg[1:100,1],neg_sto=F_neg_d[1:100,2],neg_heat=F_neg_d[1:100,3])
flex_pos|>stack|>@vlplot(:area,x=:t,y={:value,stack=:zero},color="variable:n")
flex_neg|>stack|>@vlplot(:area,x=:t,y={:value,stack=:zero},color="variable:n")

#-
test_model = BSON.load(joinpath(savepath, "run_5000.0_48_17_fixed_fg_inv.bson"))
F_pos, F_neg, F_pos_d, F_neg_d = naive_flex_potential(test_model, pv, wind, timesteps)
flex_pos = DataFrame(t = timesteps[1:end-12], cur = F_pos_d[:,1], sto = F_pos_d[:,2], heat = F_pos_d[:,3])
flex_neg = DataFrame(t = timesteps[1:end-12],neg_cur = F_neg[:,1],neg_sto=F_neg_d[:,2],neg_heat=F_neg_d[:,3])
time_window = timesteps[1:end-12]#1:1003 Cost of flexibos_d[time_window,1],pos_sto=F_pos_d[time_window,2],pos_heat=F_pos_d[time_window,3],
neg_cur = F_neg[time_window,1],neg_sto=F_neg_d[time_window,2],neg_heat=F_neg_d[time_window,3]
flex_pos|>stack|>@vlplot(:area,x=:t,y={:value,stack=:zero},color="variable:n")
flex_neg|>stack|>@vlplot(:area,x=:t,y={:value,stack=:zero},color="variable:n")
flex_df|>stack|>@vlplot(:area,x=:t,y={:value,stack=:zero},color="variable:n")
