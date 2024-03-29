## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
Pkg.instantiate()

using DataFrames
using CSV
using Clp
using Statistics;

using Random
Random.seed!(1);

#-

#=
# Second Validation example!

We validate the model by comparing produced results with results of optimization using open_plan tool (https://open-plan-tool.org/)
=#

#- 

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
The model contains a heat layer, but we specify zero heat demand for validation, as heat model in open_plan is under development.
=#

offset = 0
timesteps = 1:365*24

pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_pv_pik.csv"), DataFrame, header = false)
wind_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_wind_pik.csv"), DataFrame, header = false)
demand_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "demand_Industriepark.csv"), DataFrame, header = false)

heatdemand = zeros(length(timesteps))
pv = pv_data[timesteps .+ offset, 1]
wind = wind_data[timesteps .+ offset, 1]
demand = demand_data[timesteps .+ offset, 1]

pv_data = nothing;
wind_data = nothing;
demand_data = nothing; # Free the memory
heatdemand_data = nothing;

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

pars = copy(default_es_pars)

average_hourly_demand = mean(demand)

pars[:recovery_time] = 24
pars[:c_storage] = 300.
pars[:c_pv] = 800.
pars[:c_wind] = 2500.
pars[:c_in] = 0.4
pars[:c_out] = 0.04
pars[:asset_lifetime] = 20.;

#=
The model itself is constructed by the function define_energy_system
=#

es = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = false)

#-

#=
This contains investment variables which we collectively call $I$, and the operational schedule $O^t$.
The total cost given a certain investment $I$ and schedule $O^t$ is denoted $C(I, O^t)$.

We now can optimize the system, initialy while ignoring flexibility:
=#

sp_no_flex = instantiate(es, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_no_flex)

optimize!(sp_no_flex)

no_flex_decision = optimal_decision(sp_no_flex)

objective_value(sp_no_flex)

#-
sankey_results(sp_no_flex, pv, wind, demand)