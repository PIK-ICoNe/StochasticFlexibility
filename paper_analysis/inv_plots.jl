#=
We focus on three particular pairs of F and scen_freq:
(5000., 48), (25000., 48), (25000., 240)
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()
using DataFrames
using CSV
using BSON

#-
include(joinpath(basepath, "src", "sp_model.jl"))
#-
opt_params = [(5000., 48, 25), (25000., 48, 25), (25000., 240, 120)]
case_names = ["small", "base", "rare"]

plain_labels = ["OFR", "OFOR", "OFIOR"]
run_id = "flex_cost_07_28"
savepath = joinpath(basepath, "results", run_id)#05_12")
#-
cost_csv = joinpath(savepath, "mean_costs.csv")
invs_csv = joinpath(savepath, "mean_invs.csv")
df_cost_err = CSV.read(cost_csv, DataFrame)
df_inv_mean = CSV.read(invs_csv, DataFrame)
#-
using CairoMakie

pl_colors = Makie.wong_colors()
fig_inv = Figure(resolution = (2400, 1800))
ax = Axis(fig_inv[1,1], yticks = (0:3, vcat("Base model", string.(case_names))),
        title = "Investment decision")
labels = ["PV", "Wind", "Electrical storage", "Heat storage", "Heat pumps"]#string.(inv_vars);
elements = [PolyElement(polycolor = pl_colors[i]); for i in 1:length(labels)]
title = "Components";
Legend(fig_inv[1,2], elements, labels, title)

barlabels = ["" for i in 1:size(df_inv_mean)[1]];
for i in 1:size(df_inv_mean)[1]
    b = plain_labels
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
stack = df_inv_mean[!, :var], bar_labels = barlabels,
label_size = 35, label_offset = 20)
fig_inv
display(fig_inv)
save(joinpath(basepath, "paper_plots/opt_mode_investments.png"), fig_inv);
#-