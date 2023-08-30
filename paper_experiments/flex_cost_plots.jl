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
using JSON
using BSON
#-
include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "data_load.jl"));
#-
opt_params = [(5000., 48, 25), (25000., 48, 25), (25000., 240, 120)]
case_names = ["small", "base", "rare"]
plain_labels = ["OFR", "OFOR", "OFIOR"]
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;
run_id = "flex_cost_07_28"#"flex_cost_05_12"
savepath = joinpath(basepath, "results", run_id)
#-
cost_csv = joinpath(savepath, "mean_costs.csv")
invs_csv = joinpath(savepath, "mean_invs.csv")
df_cost_err = CSV.read(cost_csv, DataFrame)
#-
using PlotlyJS

fig_obj = PlotlyJS.Plot();
current_colorway = PlotlyJS.PlotlyBase._get_colorway(fig_obj)

p = Dict(())
for c_var in ["CR", "CF"]
    fig_layout = Layout(width=1200, height=800,
    margin=attr(l=20,r=20,t=20,b=20, autoexpand=false),
    #paper_bgcolor="white", 
    plot_bgcolor = "white",
    legend = attr(x=0.7,y=1),
    font=attr(size = 20),
    yaxis = attr(exponentformat = "power", title_text=c_var*", Euro"))
    norm_fig_layout = Layout(width=1200, height=800,
    margin=attr(l=20,r=20,t=20,b=20, autoexpand=false),
    #paper_bgcolor="white", 
    plot_bgcolor = "white",
    legend = attr(x=0.7,y=1),
    font=attr(size = 20),
    yaxis = attr(exponentformat = "power", title_text=c_var*"/Fₛ, Euro/kWh"))
    traces = [
        PlotlyJS.bar(
            x = string.(case_names),
            y = subset(df_cost_err,:mode => n-> n.==i)[!,Symbol("mean_"*c_var)],
            error_y=attr(type="data", array=subset(df_cost_err,:mode => n-> n.==i)[!,Symbol("err_"*c_var)]),
            name = plain_labels[i], 
            #(c_var == "CG") && showlegend = false,
            marker_color = current_colorway[i])
            for i in 1:3]
    p[c_var]= PlotlyJS.Plot(traces, fig_layout)
    traces_norm = [
        PlotlyJS.bar(
            x = string.(case_names),
            y = subset(df_cost_err,:mode => n-> n.==i)[!,Symbol("norm_mean_"*c_var)],
            error_y=attr(type="data", array=subset(df_cost_err,:mode => n-> n.==i)[!,Symbol("norm_err_"*c_var)]),
            name = plain_labels[i])
            for i in 1:3]
    p["norm"*c_var]=PlotlyJS.Plot(traces_norm, norm_fig_layout)
    #relayout!(p["norm"*c_var],yaxis = attr(title_text = c_var*"/Fₛ, Euro"))
end

#Plot normalized costs:
for k in keys(p)
    PlotlyJS.savefig(p[k], joinpath(basepath, "paper_plots", "flex_cost", k*"$run_id.png"))
end
#-