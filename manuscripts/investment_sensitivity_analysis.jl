## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using DataFrames
using CSV
using Clp
using Tables
using Statistics
using StochasticPrograms;

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
pars[:penalty] = 1000000.;

#=
We now do the optimization for the system without any flexibility events. The background system:
=#

n_runs = 1

# Now we set up the system with the same set of parameters and sample size n_runs times.

n_samples = 160
#-


t_max = length(pv) - 24
F_max = 10000.
delta_t = 7*24 - recovery_time
pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1);


# n_samples = 80

t_max = length(pv) - 24
F_max = 10000.
F_min = 3000.
delta_t = 7*24 - recovery_time
pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1);

n = round(Int, n_samples * pars[:event_per_scen])
println("$n total scenarios, with an average of $(pars[:event_per_scen]) events per full time period")

#-
investments = []
results = []
all_scens = []

for i in 1:n_runs
    stime = time()
    scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
    es = define_energy_system(pv, wind, demand, heatdemand; p = pars)
    sp = instantiate(es, scens, optimizer = Clp.Optimizer)
    set_silent(sp)
    println("Sample $i instantiated")
    println("Model setup and instantiation performed in $(time() - stime) seconds")
    stime = time()
    optimize!(sp)
    println("Model optimized in $(time() - stime) secons")
    push!(investments, get_investments(sp))
    push!(results, objective_value(sp))
    push!(all_scens, scens)
end
#-
# Get background investment decision for comparison

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_event_per_scen = true)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_bkg)
optimize!(sp_bkg)
cost_bkg = objective_value(sp_bkg)
bkg_investments = get_investments(sp_bkg)

#-
mean_investments = Dict([(var => mean([inv[var] for inv in investments])) for var in keys(bkg_investments)])

#-
stime = time()
scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
es_mean = define_energy_system(pv, wind, demand, heatdemand; p = pars)
sp_mean = instantiate(es_mean, scens, optimizer = Clp.Optimizer)
fix_investment!(sp_mean, mean_investments)
optimize!(sp_mean)
println("Optimization performed in $(time() - stime) seconds")

sp_mean_cost = objective_value(sp_mean)


#-
stime = time()
scens_resampled = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
sp_mean_resampled = instantiate(es_mean, scens_resampled, optimizer = Clp.Optimizer)
set_silent(sp_mean_resampled)

cost_resampled = evaluate_decision(sp_mean_resampled, optimal_decision(sp_mean))
println("Evaluation performed in $(time() - stime) seconds")

#-
println("Cost of system with mean investment decision relative to background cost: $(sp_mean_cost/cost_bkg)")
println("Cost of system with fixed first stage and new sample relative to background cost: $(cost_resampled/cost_bkg)")

# The cost after resampling is still much higher than the background cost

#-

for inv in investments
    CSV.write(joinpath(basepath, "results", "investment_sensitivity.csv"), DataFrame(investments), append = true)
end

CSV.write(joinpath(basepath, "results", "cost_sensitivity.csv"), Tables.table(results), append = true)