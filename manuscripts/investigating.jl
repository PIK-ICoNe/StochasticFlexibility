## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using DataFrames
using CSV
using XLSX
using Clp
using Statistics;

using Random
Random.seed!(1);

#-


#=
# Validation

We validate the model by comparing produced results with results of optimization using open_plan tool (https://open-plan-tool.org/)
=#

#- 

include(joinpath(basepath, "src", "sp_model_old_storage.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
The model contains a heat layer, but we specify zero heat demand for validation, as heat model in open_plan is under development.
=#

offset = 0
timesteps = 1:365*24

data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)
heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)

pv = data[timesteps .+ offset, 3]
wind = data[timesteps .+ offset, 4]
demand = data[timesteps .+ offset, 2]
heatdemand = heatdemand_data[timesteps .+ offset, 1]

plt = plot(timesteps, pv .* (mean(demand) / mean(pv)), label="PV (unitless)")
plot!(plt, timesteps, wind.* (mean(demand) / mean(wind)), label="Wind (unitless)")
plot!(plt, timesteps, heatdemand, label="Heat Demand")
plot!(plt, timesteps, demand, label="Electric Demand")
plt
#-

#=
Next we continue the set up. Our model comes with default parameters,
which we adjust to match parameters used in the open_plan optimization.
=#

# pars = copy(default_es_pars)

# average_hourly_demand = mean(demand)

# pars[:recovery_time] = 24
# pars[:c_storage] = 300.
# pars[:c_pv] = 800.
# pars[:c_wind] = 2500.
# pars[:c_in] = 0.9
# pars[:c_out] = 0.0001
# pars[:asset_lifetime] = 20.
# pars[:sto_ef_ch] = 0.95
# pars[:sto_ef_dis] = 0.95

# pars[:feedincap] = 1e7;

#=
The model itself is constructed by the function define_energy_system
=#


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
delta_t = 72 # Flex event every three days
pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);
n = round(Int, 20 * pars[:scens_in_year])

scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max)

es = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true)

#-

#=
This contains investment variables which we collectively call $I$, and the operational schedule $O^t$.
The total cost given a certain investment $I$ and schedule $O^t$ is denoted $C(I, O^t)$.

We now can optimize the system, initialy while ignoring flexibility:
=#

sp_no_flex = instantiate(es, scens, optimizer = Clp.Optimizer)
set_silent(sp_no_flex)

optimize!(sp_no_flex)

no_flex_decision = optimal_decision(sp_no_flex)

objective_value(sp_no_flex)

#-
sankey_results(sp_no_flex, pv, wind, demand, timesteps)
#-
plot_results(sp_no_flex, pv, wind, demand)

#-

vars = ["gci", "gco", "sto_soc", "sto_to_bus", "heat_sto_soc", "flow_energy2heat"]

plot_scenario_debug(sp_no_flex, 10; vars=vars)

#- 

plot_outcome(sp_no_flex, 1100, -1, 0.)
#-
@show pars[:scens_in_year] 
@show value.(sp_no_flex[1,:u_pv])
@show value.(sp_no_flex[1,:u_wind])
@show value.(sp_no_flex[1,:u_storage])
@show sum(value.(sp_no_flex[1,:gci]))
@show sum(value.(sp_no_flex[1,:gco]))
@show sum(value.(sp_no_flex[1,:sto_from_bus]))
@show sum(value.(sp_no_flex[1,:sto_to_bus]))
