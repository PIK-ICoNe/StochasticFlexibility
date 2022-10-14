## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using DataFrames
using CSV
using Tables
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
pars[:c_sto_op] = 0.00001;
pars[:penalty] = 1000000.

#=
We now do the optimization for the system without any flexibility events. The background system:
=#

# Get background investment decision for comparison

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_event_per_scen = true)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_bkg)
optimize!(sp_bkg)
cost_bkg = objective_value(sp_bkg)
bkg_investments = get_investments(sp_bkg)

CSV.write(joinpath(basepath, "results", "investments_bkg.csv"), DataFrame(bkg_investments), header = string.(keys(bkg_investments)), append = true)
CSV.write(joinpath(basepath, "results", "cost_bkg.csv"), Tables.table([cost_bkg]), header = ["Objective value"])

n_runs = 10
t_max = length(pv) - 24
F_max = 10000.
delta_t = 7*24 - recovery_time
pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1);
F_min = 3000.
#-
# Now we set up the system with the same set of parameters and sample size n_runs times.
for n_samples in 40:40:200
    n = round(Int, n_samples * pars[:event_per_scen])
    println("$n total scenarios, with an average of $(pars[:event_per_scen]) events per full time period")

    stime = time()
    scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
    es = define_energy_system(pv, wind, demand, heatdemand; p = pars)
    sp = instantiate(es, scens, optimizer = Clp.Optimizer)
    set_silent(sp)
    #println("Sample $i instantiated")
    println("Model setup and instantiation performed in $(time() - stime) seconds")
    stime = time()
    optimize!(sp)
    println("Model optimized in $(time() - stime) secons")
    investments = get_investments(sp)
    CSV.write(joinpath(basepath, "results", "investments$n.csv"), DataFrame(investments), header = string.(keys(investments)), append = true)
    CSV.write(joinpath(basepath, "results", "costs$n.csv"), Tables.table([objective_value(sp)]), header = ["Objective value"], append = true)
end
