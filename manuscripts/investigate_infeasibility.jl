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

pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "pv_Halle18.csv"), DataFrame, header = false)
wind_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "wind_Karholz.csv"), DataFrame, header = false)
demand_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "demand_Industriepark.csv"), DataFrame, header = false)

heatdemand = zeros(length(timesteps))
pv = pv_data[timesteps .+ offset, 1]
wind = wind_data[timesteps .+ offset, 1]
demand = demand_data[timesteps .+ offset, 1]

pv_data = nothing;
wind_data = nothing;
demand_data = nothing; # Free the memory

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
pars[:c_wind] = 1150.
pars[:c_in] = 0.165
pars[:c_out] = 0.02
pars[:asset_lifetime] = 20.
pars[:investment_budget] = 10000000.
pars[:feedincap] = 1e7;

#=
The model itself is constructed by the function define_energy_system
=#

pars[:scens_in_year] = 100.;

es = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true)

#-

#=
This contains investment variables which we collectively call $I$, and the operational schedule $O^t$.
The total cost given a certain investment $I$ and schedule $O^t$ is denoted $C(I, O^t)$.

We now can optimize the system, initialy while ignoring flexibility:
=#

sp_no_flex = instantiate(es, simple_flex_sampler(1,100.,2), optimizer = Clp.Optimizer)
set_silent(sp_no_flex)

optimize!(sp_no_flex)

no_flex_decision = optimal_decision(sp_no_flex)

objective_value(sp_no_flex)

#-
sankey_results(sp_no_flex, pv, wind, demand, [2])
#-
plot_results(sp_no_flex, pv, wind, demand)

#-
plot_outcome(sp_no_flex, 1, 1, 0., window_start = 1)

#-
@show pars[:scens_in_year] 
@show value.(sp_no_flex[1,:u_pv])
@show value.(sp_no_flex[1,:u_wind])
@show value.(sp_no_flex[1,:u_storage])
@show sum(value.(sp_no_flex[1,:gci]))
@show sum(value.(sp_no_flex[1,:gco]))
@show sum(value.(sp_no_flex[1,:sto_from_bus]))
@show sum(value.(sp_no_flex[1,:sto_to_bus]))

#-