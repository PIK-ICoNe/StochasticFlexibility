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
# include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#

timesteps = 1:(24*365)
timesteps = 1:(6*365)

data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)
heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)

pv = data[timesteps, 3]
wind = data[timesteps, 4]
demand = data[timesteps, 2]
heatdemand = heatdemand_data[timesteps, 1]

data = nothing; # Free the memory
heatdemand_data = nothing;

# plt = plot(timesteps, pv .* (mean(demand) / mean(pv)), label="PV (unitless)")
# plot!(plt, timesteps, wind.* (mean(demand) / mean(wind)), label="Wind (unitless)")
# plot!(plt, timesteps, heatdemand, label="Heat Demand")
# plot!(plt, timesteps, demand, label="Electric Demand")
# plt
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
stime = time()

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_scens_in_year = true)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_bkg)
optimize!(sp_bkg)

bkg_decision = optimal_decision(sp_bkg)

cost_bkg = objective_value(sp_bkg)

println("Background model without flexibility set up and optimized in $(time() - stime) seconds")

#-

# plot_results(sp_bkg, pv, wind, demand)

#-

# plot_heat_layer(sp_bkg, heatdemand)
#-

#=
Now we set up the system for different number of samples twice, the first for optimizing,
the second as a resampled version for validating.
=#
#for n_samples in [40, 80]
results = []
n_samples = 160
#-

stime = time()

t_max = length(pv) - 24
F_max = 10000.
delta_t = 7*24 - recovery_time
pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);


# n_samples = 80

t_max = length(pv) - 24
F_max = 10000.
F_min = 3000.
delta_t = 7*24 - recovery_time
pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);

n = round(Int, n_samples * pars[:scens_in_year])
println("$n total scenarios, with an average of $(pars[:scens_in_year]) events per full time period")

scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
scens_resampled = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
es = define_energy_system(pv, wind, demand, heatdemand; p = pars)
# sp = instantiate(es, scens, optimizer = LShaped.Optimizer)
# set_optimizer_attribute(sp, MasterOptimizer(), Clp.Optimizer)
# set_optimizer_attribute(sp, SubProblemOptimizer(), Clp.Optimizer)
# sp_resampled = instantiate(es, scens_resampled, optimizer = LShaped.Optimizer)
# set_optimizer_attribute(sp_resampled, MasterOptimizer(), Clp.Optimizer)
# set_optimizer_attribute(sp_resampled, SubProblemOptimizer(), Clp.Optimizer)

sp = instantiate(es, scens, optimizer = Clp.Optimizer)
sp_resampled = instantiate(es, scens_resampled, optimizer = Clp.Optimizer)
set_silent(sp)
set_silent(sp_resampled)

println("Model setup and instantiation performed in $(time() - stime) seconds")

#-

# Optimize the first set of models...

stime = time()

optimize!(sp)

cost = objective_value.(sp)

println("Optimization performed in $(time() - stime) seconds")
#-
# Evaluate on the second set.

stime = time()

cost_resampled = evaluate_decision(sp_resampled, optimal_decision(sp))

println("Evaluation performed in $(time() - stime) seconds")

#-

println("Cost of optimal system with foresight relative to background:")
println("$(cost/cost_bkg)")
println("Cost of optimal system without foresight relative to background:")
println("$(cost_resampled/cost_bkg)")

#-
push!(results, (n_samples, cost/cost_bkg, cost_resampled/cost_bkg))

#end

#-
stime = time()
operation = get_operation(sp)
investments = get_investments(sp)
fix_operation!(sp_resampled, operation, length(pv))
fix_investment!(sp_resampled, investments)
optimize!(sp_resampled)

println("Fixing the first stage and optimizing the system performed in $(time() - stime) seconds")

#-
#get_penalty_array(sp_resampled)
penalized_resampled = get_penalized_scenarios(sp_resampled)
penalized = get_penalized_scenarios(sp)

println("Penalized scenarios in optimal system with foresight: $(length(penalized))")
println("Penalized scenarios in optimal system without foresight: $(length(penalized_resampled))")

#-
