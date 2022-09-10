## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using DataFrames
using CSV
using Clp
using Statistics;
using StochasticPrograms

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
pars[:c_sto_op] = 0.00001;
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
Now we set up the system for different number of samples twice, the first for optimizing,
the second as a resampled version for validating.
=#

#-
t_max = length(pv) - 24
F_max = 10000.
delta_t = 7*24 - recovery_time
pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);


n_samples = [collect(2:2:10); collect(15:5:25)]

sps = [sp_bkg for n in n_samples]
sps_resampled = [sp_bkg for n in n_samples]

t_max = length(pv) - 24
F_max = 10000.
delta_t = 7*24 - recovery_time
pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);

for i in eachindex(n_samples)
    n = round(Int, n_samples[i] * pars[:scens_in_year])
    scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max)
    scens_resampled = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max)
    es = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_scens_in_year = true)
    sp = instantiate(es, scens, optimizer = Clp.Optimizer)
    sp_resampled = instantiate(es, scens_resampled, optimizer = Clp.Optimizer)
    set_silent(sp)
    set_silent(sp_resampled)
    sps[i] = sp
    sps_resampled[i] = sp_resampled
end

#-

# Optimize the first set of models...

stime = time()

Threads.@threads for i in eachindex(sps)
    optimize!(sps[i])
end

costs = objective_value.(sps)

println("Optimization performed in $(time() - stime) seconds")
#-
# Evaluate on the second set.

costs_resampled = zeros(length(n_samples))

Threads.@threads  for i in eachindex(n_samples)
    costs_resampled[i] = evaluate_decision(sps_resampled[i], optimal_decision(sps[i]))
end

#-

cost_plot = plot()
cost_plot = plot!(n_samples, costs ./ cost_bkg)
#-

cost_plot = plot()
cost_plot = plot!(n_samples, costs_resampled ./ cost_bkg)
#-

