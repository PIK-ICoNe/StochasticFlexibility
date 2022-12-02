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
using Plots

#= Include the file containing all the dependecies and optimization functions.
=#
#-
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));

#-
warm_up()
#-
bkg_cost = CSV.read(joinpath(basepath, "results/bkg", "cost_bkg.csv"), DataFrame)[!,1][1]

n_samples = 10:10:70 #[collect(2:2:10); collect(15:10:25)] #[10]
scen_freq = 72:48:24*30#24:24:30*24 #[24]

flex_cost = zeros(length(n_samples),length(scen_freq));
runtime = zeros(length(n_samples),length(scen_freq));
investments = zeros(length(n_samples),length(scen_freq));
#-
u_wind_data = [];
wind_invs = zeros(length(n_samples),length(scen_freq));
storage_invs = zeros(length(n_samples),length(scen_freq));

run_id = "2022-10-31"
for i in eachindex(n_samples)
     for j in eachindex(scen_freq)
        event_per_scen = 365*24 / (scen_freq[j] + 1)
        costs = load_costs(run_id, n_samples[i], scen_freq[j])
        mean_cost = mean(costs)
        invs, mean_invs = load_invs(run_id, n_samples[i], scen_freq[j])
        wind_invs[i,j] = mean_invs[:u_wind]
        storage_invs[i,j] = mean_invs[:u_storage]#
        #plot!(plt_inv, ([scen_freq[i] for n in 1:10],invs[:u_wind]), seriestype = :scatter )
        append!(u_wind_data, [invs[:u_wind]])
        runtimes = load_runtime(run_id, n_samples[i], scen_freq[j])
        mean_runtime = mean(runtimes)
        # flex_cost - cost of optimized system with mean investment and new sample TODO
        flex_cost[i,j] = (mean_cost - bkg_cost)/round(event_per_scen*n_samples[i])
        runtime[i,j] = mean_runtime
    end
end

heatmap(scen_freq, n_samples, flex_cost,
    xlabel = "Average offset between events",
    ylabel = "Scaled sample size", 
    title = "Cost of flexibility")
#-
heatmap(scen_freq, n_samples[2:end], runtime[2:end, :], 
    xlabel = "Average offset between events",
    ylabel = "Scaled sample size", 
    title = "Runtime, s")

heatmap(scen_freq, n_samples[2:end], wind_invs[2:end, :], 
    xlabel = "Average offset between events",
    ylabel = "Scaled sample size", 
    title = "Wind, kWp")

heatmap(scen_freq, n_samples[2:end],storage_invs[2:end, :], 
    xlabel = "Average offset between events",
    ylabel = "Scaled sample size", title = "Storage, kWh")

#-
#=
invs, mean_invs = load_invs(run_id, n_samples[5], scen_freq[5])
plt = plot()

plot(([scen_freq[5] for n in 1:10],invs[:u_wind]), seriestype = :scatter)

invs, mean_invs = load_invs(run_id, n_samples[6], scen_freq[5])

plot(([scen_freq[6] for n in 1:10],invs[:u_wind]), seriestype = :scatter)

plot!(plt,(5,10), seriestype = :scatter)

=#