## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using DataFrames
using CSV
using Clp
using Statistics;

using Random
Random.seed!(1);

#-

#=
# Evaluating the benefit of planning for flexibility events
=#

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#

timesteps = 1:(24*365)

data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)
heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)

pv = data[timesteps, 3]
wind = data[timesteps, 4]
demand = data[timesteps, 2]
heatdemand = heatdemand_data[timesteps, 1]

data = nothing; # Free the memory
heatdemand_data = nothing;

plt = plot(timesteps, pv .* (mean(demand) / mean(pv)), label="PV (unitless)")
plot!(plt, timesteps, wind.* (mean(demand) / mean(wind)), label="Wind (unitless)")
plot!(plt, timesteps, heatdemand, label="Heat Demand")
plot!(plt, timesteps, demand, label="Electric Demand")
plt
#-

#=
Next we continue the set up. Our model comes with default parameters,
which we slightly adjust here. We use some arbitrary values to define a dummy heat demand.
=#

pars = copy(default_es_pars)

recovery_time = 12

pars[:recovery_time] = recovery_time
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 550.
pars[:penalty] = 1000000.

#=
We now do the optimization for the system without any flexibility events. The background system:
=#

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_scens_in_year = true)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_bkg)
optimize!(sp_bkg)

bkg_decision = optimal_decision(sp_bkg)

cost_bkg = objective_value(sp_bkg)

#-

plot_results(sp_bkg, pv, wind, demand)

#-

plot_heat_layer(sp_bkg, heatdemand)
#-

#=
Now we evaluate the system for different frequencies of flexibility
=#

#-

flex_interval = 18:-5:3
sps = []

Threads.@threads for i in eachindex(flex_interval)
    t_max = length(pv) - 24
    F_max = 10000.
    delta_t = flex_interval[i]*24 - recovery_time
    pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);
    n = round(Int, 10 * pars[:scens_in_year])
    scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max)
    es = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_scens_in_year = true)
    sp = instantiate(es, scens, optimizer = Clp.Optimizer)
    set_silent(sp)
    push!(sps, sp)
end

#-

#=
Optimize only the operational schedule for flexibility
=#

stime = time()

investments_bkg = get_investments(sp_bkg)

Threads.@threads for i in eachindex(sps)
    fix_investment!(sps[i], investments_bkg)
    optimize!(sps[i])
end

flex_cost_op = objective_value.(sps)

println("Optimization performed in $(time() - stime) seconds")
#-

sps_op = deepcopy(sps);

#-
#=
Optimize investments and the operational schedule for flexibility
=#

stime = time()

Threads.@threads for i in eachindex(sps)
    unfix_investment!(sps[i], investments_bkg)
    optimize!(sps[i])
end

flex_cost = objective_value.(sps)

println("Optimization performed in $(time() - stime) seconds")

#-

cost_plot = plot()
cost_plot = plot!(cost_plot, 1 ./ flex_interval, flex_cost ./ cost_bkg, label = "flex_cost")
cost_plot = plot!(cost_plot, 1 ./ flex_interval, flex_cost_op ./ cost_bkg, label = "flex_cost_op")
cost_plot
savefig(cost_plot, "costplot.png")
plot!(cost_plot, ylimits=(0.9995,1.0005))
savefig(cost_plot, "costplot2.png")

#-


println(investments_bkg)

println(get_investments(sps[end]))

# TODO print how the optimal level of investment depends on the frequency of flexibility
