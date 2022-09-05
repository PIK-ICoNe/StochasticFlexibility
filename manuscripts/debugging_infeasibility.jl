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

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));

#-


offset = 0
timesteps = 1:(24*365)

data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)
heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)

pv = data[timesteps .+ offset, 3]
wind = data[timesteps .+ offset, 4]
demand = data[timesteps .+ offset, 2]
heatdemand = heatdemand_data[timesteps .+ offset, 1]

data = nothing; # Free the memory
heatdemand_data = nothing;
#-

pars = copy(default_es_pars)

average_hourly_demand = mean(demand)

recovery_time = 12

pars[:recovery_time] = recovery_time
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 550.
pars[:c_sto_op] = 0.00001;

t_max = length(pv) - 24
F_max = 10000.
delta_t = 7*24 # Flex event every week
pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);
n = round(Int, 10 * pars[:scens_in_year])

scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max)

#-


es_1_scen = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true, override_no_scens_in_year = true)
es = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true)

#-

sp_no_flex = instantiate(es_1_scen, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)

set_silent(sp_no_flex)

optimize!(sp_no_flex)

no_flex_decision = optimal_decision(sp_no_flex)

objective_value(sp_no_flex)

#-

es_reg = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = 0.)

sp_reg_flex = instantiate(es_reg, scens, optimizer = Clp.Optimizer)
set_silent(sp_reg_flex)

#-

evaluate_decision(sp_reg_flex, no_flex_decision)

#-

res = [evaluate_decision(sp_reg_flex, no_flex_decision, scens[i]) for i in 1:5]

#-

