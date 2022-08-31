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

pars[:recovery_time] = 24
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 550.
pars[:c_sto_op] = 0.00001;

pars[:scens_in_year] = 52.;

es = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true)

#-

#=
This contains investment variables which we collectively call $I$, and the operational schedule $O^t$.
The total cost given a certain investment $I$ and schedule $O^t$ is denoted $C(I, O^t)$.

We now can optimize the system, initialy while ignoring flexibility:
=#
# no_flex_pseudo_sampler
sp_no_flex = instantiate(es, simple_flex_sampler(20*pars[:scens_in_year],10000, length(pv) - 24), optimizer = Clp.Optimizer)
set_silent(sp_no_flex)

optimize!(sp_no_flex)

no_flex_decision = optimal_decision(sp_no_flex)

objective_value(sp_no_flex)

#-
#sankey_results(sp_no_flex, pv, wind, demand, timesteps)
#-
plot_results(sp_no_flex, pv, wind, demand, plot_span = 1600:1700)

#-

vars = ["gci", "gco", "sto_soc", "sto_to_bus", "heat_sto_soc", "flow_energy2heat"]

plot_scenario_debug(sp_no_flex, 100; vars=vars)

#- 

plot_outcome(sp_no_flex, 190, -1, 50., window_start=-1)
#-
@show pars[:scens_in_year] 
@show value.(sp_no_flex[1,:u_pv])
@show value.(sp_no_flex[1,:u_wind])
@show value.(sp_no_flex[1,:u_storage])
@show sum(value.(sp_no_flex[1,:gci]))
@show sum(value.(sp_no_flex[1,:gco]))
@show sum(value.(sp_no_flex[1,:sto_from_bus]))
@show sum(value.(sp_no_flex[1,:sto_to_bus]))

@show value.(sp_no_flex[2,:sto_soc],1)[1] - value.(sp_no_flex[2,:sto_soc2],1)[1]
@show value.(sp_no_flex[2,:sto_soc],1)[pars[:recovery_time]] - value.(sp_no_flex[2,:sto_soc2],1)[pars[:recovery_time]]
