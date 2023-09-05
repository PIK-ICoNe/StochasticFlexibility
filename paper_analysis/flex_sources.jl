#=
We focus on three particular pairs of F and scen_freq:
(5000., 48), (25000., 48), (25000., 240)
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
#-
using BSON

include(joinpath(basepath, "src", "data_load.jl"));
include(joinpath(basepath, "src", "evaluation_utils.jl"))
include(joinpath(basepath, "src", "sp_model.jl")) # needed for default_es_pars, could be better to move it somewhere else?
#-
opt_params = [(500., 96, 20)]
case_names = ["base"]
plain_labels = ["OFR", "OFOR", "OFIOR"]
n_runs = 3
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(heat = "when2heat");
run_id = "new_parameters"
savepath = joinpath(basepath, "results", run_id)
model_names = Dict(("fixed_fg"=> "OFR", "fixed_fg_inv"=>"OFOR", "OF"=>"OFIOR"))
#-

subplot_titles=["a" "b" "c"; "d" "e" "f"; "g" "h" "i"]

using CairoMakie
# Create the figure and axis
set_theme!(Theme(fontsize=35))
fig = Figure(resolution = (2400, 1800))
fig_temp = Figure()
ax = [Axis(fig_temp[1,1]) for i in 1:9]
ax = reshape(ax, (3,3))

plot_window = 1:100
pl_colors = [(:red,0.5), (:blue,0.5), (:green,0.5)]
plot_labels = ["Curtailment","Electrical storage", "Heat storage"]
for p in eachindex(opt_params)
    F, scen_freq, n_samples = opt_params[p]
    for (m_index, m) in enumerate(plain_labels)
        Ax = Axis(fig[p,m_index])
        test_model = BSON.load(joinpath(savepath, "run_$(F)_$(scen_freq)_1_$(m).bson"))
        F_pos, F_neg, F_pos_d, F_neg_d = naive_flex_potential(test_model, pv, wind, timesteps)
        F_pos_stack = hcat(zeros(size(F_pos_d)[1]),cumsum(F_pos_d, dims=2))
        F_neg_stack = hcat(zeros(size(F_neg_d)[1]),cumsum(F_neg_d, dims =2))
        for f in 1:3 # iterate over flex sources: cur, sto, heat
            fill_between!(Ax,plot_window, 
            F_pos_stack[plot_window,f], F_pos_stack[plot_window,f+1], 
            color=pl_colors[f])
            fill_between!(Ax,plot_window, 
            F_neg_stack[plot_window,f+1], F_neg_stack[plot_window,f], 
            color=pl_colors[f])
        
        end
        if p == length(opt_params) # last row, work on x axis
            Ax.xlabel = "t, h"
        else
            Ax.xticklabelsvisible = false
        end
        if m_index == 1 # first column, set y axis
            Ax.ylabel = "F, kW"
        else
            Ax.yticklabelsvisible = false
        end
        if p == 1
            Ax.title = m
        end
        lines!(Ax,plot_window, fill(F, plot_window), 
            color=:black, linewidth=1.5, linestyle = :dash)
            lines!(Ax,plot_window, fill(-F, plot_window),
            color=:black, linewidth=1.5, linestyle = :dash)
            text!(Ax, 0,1, 
            text=subplot_titles[p,m_index], align = (:left,:top),
            fontsize = 14, font = :bold, space = :relative,
            offset = (4,-2))
            ax[p, m_index] = Ax
    end
end
linkyaxes!(ax[1,1], ax[1,2], ax[1,3])
linkyaxes!(ax[2,1], ax[2,2], ax[2,3])
linkyaxes!(ax[3,1], ax[3,2], ax[3,3])

labels = plot_labels
elements = [PolyElement(polycolor = pl_colors[i]); for i in 1:length(plot_labels)]
for p in eachindex(opt_params)
    Label(fig[p,4], case_names[p], rotation=-pi/2, tellheight=false, 
    font=:bold, width=25, justification=:left, halign= :left, valign =:center)
end
l = Legend(fig[2,5], elements, labels)
display(fig)
plot_path = joinpath(basepath, "paper_plots", run_id)
if !isdir(plot_path)
    mkdir(plot_path)
end
save(joinpath(plot_path, "flex_sources_combined.png"), fig)
#-