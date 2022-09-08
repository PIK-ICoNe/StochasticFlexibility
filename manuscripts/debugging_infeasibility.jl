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
pars[:c_i] = .4;

t_max = length(pv) - 24
F_max = 10000.
delta_t = 7*24 # Flex event every week
pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);
n = round(Int, 10 * pars[:scens_in_year])

scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max)

#-


es_1_scen = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = false, override_no_scens_in_year = true)
es = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = false)

#-

sp_no_flex = instantiate(es_1_scen, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)

set_silent(sp_no_flex)

optimize!(sp_no_flex)

no_flex_decision = optimal_decision(sp_no_flex)

objective_value(sp_no_flex)

#-

es_reg = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = true)

sp_reg_flex = instantiate(es_reg, scens, optimizer = Clp.Optimizer)
set_silent(sp_reg_flex)

#-

evaluate_decision(sp_reg_flex, no_flex_decision)

#-

# Found a scenario that's infinity:

evaluate_decision(sp_reg_flex, no_flex_decision, scens[71])

#-

t_xi = scens[71].data[:t_xi]
s_xi = scens[71].data[:s_xi]
F_xi = scens[71].data[:F_xi]

plot_window = t_xi-10:t_xi+30
#-
plot_results(sp_no_flex, pv, wind, demand; plot_span = plot_window)
#-
plot_heat_layer(sp_no_flex, heatdemand; plot_span = plot_window)

# Notable: There is no heat storage use in this period...

#-

# plot_outcome(sp_no_flex, t_xi, s_xi, F_xi) # No Optimum...
#plot_results(sp_reg_flex, pv, wind, demand; plot_span = plot_window, s = 71)

infeasible_scenario_model = outcome_model(sp_reg_flex, optimal_decision(sp_no_flex), scens[71]; optimizer = subproblem_optimizer(sp_reg_flex))
set_silent(infeasible_scenario_model)
optimize!(infeasible_scenario_model) # Primal infeasible
termination_status(infeasible_scenario_model)
#-

l = list_of_constraint_types(infeasible_scenario_model)
all_constraints(infeasible_scenario_model, l[1]...) # The known values
all_constraints(infeasible_scenario_model, l[2]...) # The larger then zero bounds
all_constraints(infeasible_scenario_model, l[3]...) # The smaller than debug_cap bounds
all_constraints(infeasible_scenario_model, l[5]...) # The smaller than u bounds
all_constraints(infeasible_scenario_model, l[4]...) # The flow equations and everything else...

# l[4] are the interesting constraints...
#-
#delete(infeasible_scenario_model, feedin_constr)

# Idea: delete constraints one by one to see what works...

#-
#=
for i in 1:length(all_constraints(infeasible_scenario_model, l[4]...))
    println("Deleting $(all_constraints(infeasible_scenario_model, l[4]...)[i])")
    ism = outcome_model(sp_reg_flex, optimal_decision(sp_no_flex), scens[71]; optimizer = subproblem_optimizer(sp_reg_flex))
    set_silent(ism)
    delete(ism, all_constraints(ism, l[4]...)[i])
    optimize!(ism)
    println(termination_status(ism))
end

#-

value.(sp_no_flex[1, :heat_sto_to_bus])[4933]
value.(sp_no_flex[1, :heat_sto_from_bus])[4933]
=#