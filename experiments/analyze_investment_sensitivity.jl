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


#n_samples = 10:10:70 #[collect(2:2:10); collect(15:10:25)] #[10]
#scen_freq = 72:48:24*30#24:24:30*24 #[24]

#n_samples = [collect(15:5:35); collect(40:10:80); 100]
#scen_freq = 48:24:120

n_samples = [collect(15:5:35);collect(40:10:90);collect(100:25:175);collect(200:50:300)]
scen_freq = [96]

flex_cost = zeros(length(n_samples),length(scen_freq));
runtime = zeros(length(n_samples),length(scen_freq));
investments = zeros(length(n_samples),length(scen_freq));
#-
u_wind_data = [];
wind_invs = zeros(length(n_samples),length(scen_freq));
storage_invs = zeros(length(n_samples),length(scen_freq));

wind_full = [];
storage_full = [];

run_id = "2022-12-23"
for i in eachindex(n_samples)
     for j in eachindex(scen_freq)
        event_per_scen = 365*24 / (scen_freq[j] + 1)
        costs = load_costs(run_id, n_samples[i], scen_freq[j])
        mean_cost = mean(costs)
        invs, mean_invs = load_invs(run_id, n_samples[i], scen_freq[j])
        wind_invs[i,j] = mean_invs[:u_wind]
        storage_invs[i,j] = mean_invs[:u_storage]#
        if i == 7
            push!(wind_full,(scen_freq[j], invs[:u_wind]))
            push!(storage_full,(scen_freq[j], invs[:u_storage]))
        end
            #plot!(plt_inv, ([scen_freq[i] for n in 1:10],invs[:u_wind]), seriestype = :scatter )
        append!(u_wind_data, [invs[:u_wind]])
        runtimes = load_runtime(run_id, n_samples[i], scen_freq[j])
        mean_runtime = mean(runtimes)
        # flex_cost - cost of optimized system with mean investment and new sample TODO
        flex_cost[i,j] = (mean_cost - bkg_cost)/round(event_per_scen)#*n_samples[i]
        runtime[i,j] = mean_runtime
    end
end
#-
#=p = plot(heatmap(z=flex_cost, y=n_samples,x=[string(s) for s in scen_freq]),
    Layout(yaxis_title="Scaled sample size", xaxis_title="Average offset between flexibility events, hours",
    xaxis_nticks = 4))
savefig(p, joinpath(basepath,"flex_cost.png"))=#
#
heatmap(scen_freq, n_samples, flex_cost,
    xlabel = "Average offset between events",
    ylabel = "Scaled sample size") 
    #title = "Cost of flexibility")
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
bkg_inv = Dict(pairs(eachcol(CSV.read(joinpath(basepath, "results/bkg", "investments_bkg.csv"), DataFrame))))
bkg_inv = Dict(([v => mean(bkg_inv[v]) for v in keys(bkg_inv)]))
bkg_op = Dict(pairs(eachcol(CSV.read(joinpath(basepath, "results/bkg", "op_bkg.csv"), DataFrame))))

timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps);

pars[:sto_ef_ch] = 0.97
pars[:sto_ef_dis] = 0.97
pars[:recovery_time] = 12

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_event_per_scen = true)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_bkg)
fix_investment!(sp_bkg, bkg_inv)
fix_operation!(sp_bkg, bkg_op, length(timesteps))
optimize!(sp_bkg)

#-
objective_value(sp_bkg)
t_max_offset = 24
t_max = minimum((length(pv), length(wind), length(demand), length(heatdemand))) - t_max_offset
recovery_time = pars[:recovery_time]
delta_t = 72 - recovery_time
pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1)
n = round(Int, 1 * pars[:event_per_scen])
println("$n total scenarios, with an average of $(pars[:event_per_scen]) events per full time period")
stime = time()
scens = poisson_events_with_offset(n, delta_t, recovery_time, 10000., t_max, F_min = 3000.)
ed = [evaluate_decision(sp_bkg, optimal_decision(sp_bkg), s) for s in scens]
#-
using PlotlyJS
wind_box = [box(y=wind_full[i][2], name = "$(scen_freq[i])") for i in eachindex(scen_freq)]
storage_box = [box(y=storage_full[i][2], name = "$(scen_freq[i])") for i in eachindex(scen_freq)]


PlotlyJS.plot(wind_box)
PlotlyJS.savefig(PlotlyJS.plot(storage_box,  Layout(yaxis_title="Investment in storage", xaxis_title = "Typical offset between ancillary service requests, hours" )), joinpath(basepath, "storage_inv.png"))

#-
od = optimal_decision(sp)

ov = objective_value(sp)
# We save objective value and result of evaluate_decision for each scenario
ovs = [objective_value(sp,2, i) for i in 1:length(scens)]
eds = [evaluate_decision(sp, od, scen) for scen in scens]