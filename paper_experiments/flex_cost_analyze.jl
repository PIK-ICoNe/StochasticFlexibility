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
# Fist we get C^B
base_model = JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_0.0.json"))
CB = base_model["cost"]
CIB = get_total_investment(base_model)
COB = get_operation_cost(base_model)
inv_base = base_model["inv"]
P = base_model["params"]
base_model = nothing;
#-
opt_params = [(5000., 48, 25), (25000., 48, 25), (25000., 240, 120)]
n_runs = 5
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;
savepath = joinpath(basepath, "results/flex_cost_06_30")#05_12")
#-
df_costs = DataFrame(CSV.File(joinpath(savepath, "val_costs.csv")))
df_inv = DataFrame(CSV.File(joinpath(savepath, "inv.csv")))
#-
# Find mean investments and costs
df_inv_mean = DataFrame()
df_cost_mean = DataFrame()
df_cost_err = DataFrame()
cost_vars = [:CG, :CR]
for p in eachindex(params)
    for m in 1:3
        sel = subset(df_costs, :p => a -> a .== p, :m => n-> n.==m)
        mean_CG = mean(sel[!,:cost])-CB
        mean_CR = mean(sel[!,:CR])
        err_CG = std(sel[!,:cost].-CB)
        err_CR = std(sel[!, :CR])  
        # Total amount of flexibility requested in a year: F*events_per_year
        total_flex = (length(timesteps)-24) / (opt_params[p][2] + 1)*opt_params[p][1]*opt_params[p][3] 
        #var = 1 -> CG
        append!(df_cost_mean, DataFrame(p=p, mode=m, var = 1,
            mean_value = mean_CG))
        #var = 2 -> CR
        append!(df_cost_mean, DataFrame(p=p, mode=m, var = 2,
            mean_value = mean_CR))
        norm_costs = (sel[!, :cost] .- CB)/total_flex
        norm_mean_CG = mean(norm_costs)
        norm_err_CG = std(norm_costs)
        norm_mean_CR = mean_CR/total_flex
        norm_err_CR = std(sel[!, :CR])/total_flex
        append!(df_cost_err, DataFrame(p=p, mode = m, 
            mean_CG = mean_CG, mean_CR = mean_CR,
            err_CG = err_CG, err_CR = err_CR,
            norm_mean_CG = norm_mean_CG, norm_mean_CR = norm_mean_CR,
            norm_err_CG = norm_err_CG, norm_err_CR = norm_err_CR,))

        for (var_index, inv_var) in enumerate(inv_vars)
            append!(df_inv_mean, DataFrame(p = p, mode = m, var = var_index,
                mean_value = mean(subset(df_inv, :p => a -> a .== p, 
                :m => mode -> mode .== m)[!,inv_vars[var_index]])*P[replace(string(inv_var), "u_" => "c_")]))
        end
    end
end
#-
using CairoMakie
using LaTeXStrings
pl_colors = Makie.wong_colors()
fig_inv = Figure()
ax = Axis(fig_inv[1,1], yticks = (0:3, vcat("Base model", string.(opt_params))),
        title = "Investment decision")
labels = string.(inv_vars);
elements = [PolyElement(polycolor = pl_colors[i]); for i in 1:length(labels)]
title = "Components";
Legend(fig_inv[1,2], elements, labels, title)

barlabels = ["" for i in 1:size(df_inv_mean)[1]];
for i in 1:size(df_inv_mean)[1]
    b = ["Fixed FG", "Fixed FG Inv", "OF"]#[L"OF_{|\overline{FG}}", L"OF_{|\overline{Inv}_{fg}}",L"OF"]
    if df_inv_mean[i, :var] == length(inv_vars)
        if df_inv_mean[i, :p] != 0
            barlabels[i] = b[df_inv_mean[i, :mode]]
        end
    end
end

barplot!(ax, df_inv_mean[!, :p], df_inv_mean[!, :mean_value],
dodge = df_inv_mean[!, :mode],
direction=:x,
flip_labels_at=0.85,
color = pl_colors[df_inv_mean[!, :var]],
stack = df_inv_mean[!, :var], bar_labels = barlabels)
fig_inv
save(joinpath(basepath, "paper_plots/opt_mode_investments.png"), fig_inv);
#-
using PlotlyJS
using LaTeXStrings
for c_var in ["CG", "CR"]
    fig = PlotlyJS.plot([
        PlotlyJS.bar(
            name=L"OF_{|\overline{FG}}",
            x=string.(opt_params), y=subset(df_cost_err,:mode => n-> n.==1)[!,Symbol("mean_"*c_var)],
            error_y=attr(type="data", array=subset(df_cost_err,:mode => n-> n.==1)[!,Symbol("err_"*c_var)])
        ),
        PlotlyJS.bar(
            name=L"OF_{|\overline{Inv}_{FG}}",
            x=string.(opt_params), y=subset(df_cost_err,:mode => n-> n.==2)[!,Symbol("mean_"*c_var)],
            error_y=attr(type="data", array=subset(df_cost_err,:mode => n-> n.==2)[!,Symbol("err_"*c_var)])
        ),
        PlotlyJS.bar(
            name=L"OF",
            x=string.(opt_params), y=subset(df_cost_err,:mode => n-> n.==3)[!,Symbol("mean_"*c_var)],
            error_y=attr(type="data", array=subset(df_cost_err,:mode => n-> n.==3)[!,Symbol("err_"*c_var)])
        )
    ])
    PlotlyJS.savefig(fig,joinpath(basepath, "paper_plots", c_var*".png"))

    fig_norm = PlotlyJS.plot([
        PlotlyJS.bar(
            name=L"OF_{|\overline{FG}}",
            x=string.(opt_params), y=subset(df_cost_err,:mode => n-> n.==1)[!,Symbol("norm_mean_"*c_var)],
            error_y=attr(type="data", array=subset(df_cost_err,:mode => n-> n.==1)[!,Symbol("norm_err_"*c_var)])
        ),
        PlotlyJS.bar(
            name=L"OF_{|\overline{Inv}_{FG}}",
            x=string.(opt_params), y=subset(df_cost_err,:mode => n-> n.==2)[!,Symbol("norm_mean_"*c_var)],
            error_y=attr(type="data", array=subset(df_cost_err,:mode => n-> n.==2)[!,Symbol("norm_err_"*c_var)])
        ),
        PlotlyJS.bar(
            name=L"OF",
            x=string.(opt_params), y=subset(df_cost_err,:mode => n-> n.==3)[!,Symbol("norm_mean_"*c_var)],
            error_y=attr(type="data", array=subset(df_cost_err,:mode => n-> n.==3)[!,Symbol("norm_err_"*c_var)])
        )
    ])

    PlotlyJS.savefig(fig_norm,joinpath(basepath, "paper_plots", c_var*"norm.png"))
end
#-

# Stack area chart for flex potential
flex_sources = [:pos_cur, :pos_sto, :pos_heat, :neg_cur, :neg_sto,:neg_heat]
plot_labels = Dict(("cur" => "Curtailment","sto" => "Electrical storage", "heat" => "Heat storage"))
include(joinpath(basepath, "src", "evaluation_utils.jl"))
fig_obj = PlotlyJS.Plot()
current_colorway = PlotlyJS.PlotlyBase._get_colorway(fig_obj)
plot_colors = Dict(("cur"=>current_colorway[1], "sto"=>current_colorway[2], "heat"=>current_colorway[3]))
plot_window = 1:100
for p in 1:3
    F, scen_freq, n_samples = opt_params[p]
    for m in ["fixed_fg", "fixed_fg_inv", "OF"]
        test_model = BSON.load(joinpath(savepath, "preval","run_$(F)_$(scen_freq)_2_$(m).bson"))
        F_pos, F_neg, F_pos_d, F_neg_d = naive_flex_potential(test_model, pv, wind, timesteps)
        flex_full = DataFrame(t=plot_window, pos_cur=F_pos_d[plot_window,1], pos_sto=F_pos_d[plot_window,2], pos_heat=F_pos_d[plot_window,3],
            neg_cur = F_neg[plot_window,1],neg_sto=F_neg_d[plot_window,2],neg_heat=F_neg_d[plot_window,3])
        pl = PlotlyJS.plot()
        for flex_var in flex_sources
            fl_sign, fl_src = split(String(flex_var),"_")
            pl_legend = false
            pl_label = ""
            if fl_sign == "pos"
                pl_label = plot_labels[fl_src]
                pl_legend = true
            end

            add_trace!(pl,PlotlyJS.scatter(
                x=plot_window, y=flex_full[plot_window,flex_var],
                stackgroup=fl_sign, mode="lines", hoverinfo="x+y",
                line=attr(width=0.5), name = pl_label, showlegend = pl_legend, line_color=plot_colors[fl_src]
            ))

        end
        PlotlyJS.savefig(pl,joinpath(basepath, "paper_plots", "flex_sources_$(p)_$(m).png"))
    end
end
#-