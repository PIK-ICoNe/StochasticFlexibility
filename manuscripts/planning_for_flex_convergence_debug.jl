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
using Gurobi

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
Set up system
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

#-

pars = copy(default_es_pars)

recovery_time = 12

pars[:recovery_time] = recovery_time
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 550.
pars[:c_sto_op] = 0.00001;
pars[:penalty] = 1000000.

t_max = length(pv) - 24
F_max = 10000.
delta_t = 7*24 - recovery_time
pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);

n_samples = 25

t_max = length(pv) - 24
F_max = 10000.
delta_t = 7*24 - recovery_time
pars[:scens_in_year] = t_max / (delta_t + recovery_time + 1);

es = define_energy_system(pv, wind, demand, heatdemand; p = pars)

n = round(Int, n_samples * pars[:scens_in_year])

scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max)
scens_resampled = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max)

sp = instantiate(es, scens, optimizer = Gurobi.Optimizer)
set_silent(sp)

#sp_resampled = instantiate(es, scens_resampled, optimizer = Clp.Optimizer)
#set_silent(sp_resampled)

#-

optimize!(sp; cache = true)

#-

# Search for infeasible scenarios

function find_infeasible(sp, scens) # Run this to find an infeasible scenario
    for (i, scen) in enumerate(scens)
        ed = evaluate_decision_wrapper(sp, optimal_decision(sp), scen)
        if isinf(ed)
            println("Scenario $i is infinite")
            return i, scen
        end
    end
end

#-

find_infeasible(sp, scens_resampled)

#-

scen_infeasible = scens_resampled[17]

t_i = scen_infeasible.data[:t_xi]
F_i = scen_infeasible.data[:F_xi]
s_i = scen_infeasible.data[:s_xi]

#-

plot_outcome_debug(sp, t_i, s_i, F_i)
# Whoop whoop found an infeasible one!!

#-

plot_outcome_debug(sp, t_i, s_i, 0.)
# And it's indeed some sort of bug, as F = 0. should _always_ be feasible.

#-
scen = @scenario t_xi = t_i s_xi = s_i F_xi = 0. probability = 1.
model = outcome_model(sp, optimal_decision(sp),scen, optimizer = subproblem_optimizer(sp))
optimize!(model)
compute_conflict!(model)
#-

compute_conflict!(model)
if MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
    iis_model, _ = copy_conflict(model)
    print(iis_model)
end
#-
list_of_conflicting_constraints = ConstraintRef[]
for (F, S) in list_of_constraint_types(model)[2:end]
    for con in all_constraints(model, F, S)
        if MOI.get(model, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT
            push!(list_of_conflicting_constraints, con)
        end
    end
end

#-
# A detailed look at the day in question:

all_vars=["gci", "gco", "pv_cur", "wind_cur", "sto_from_bus", "sto_to_bus", "sto_soc", "heat_sto_soc", "heat_sto_from_bus", "heat_sto_to_bus", "flow_energy2heat"]

plot_base_case_raw(sp; plot_window = t_i-2:t_i+2, vars = all_vars)

# Looks completely normal. This makes me suspect the problem is in the second stage...


#-

plot_results(sp, pv, wind, demand; plot_window = t_i-2:t_i+2)

#-

plot_heat_layer(sp, heatdemand; plot_window = t_i-2:t_i+32)


#-

# The neighbouring days:

plot_outcome_debug(sp, t_i-1, s_i, 0.)

#-

plot_outcome_debug(sp, t_i+1, s_i, F_i)

#-

plot_outcome_debug(sp, t_i+2, s_i, F_i)

#-

# So day t_i - 1 shows weird behavior where it eats the penalty by violating gco[t_xi] == gco2[1].
# even if we switch of the flex request

plot_outcome_debug(sp, t_i-1, s_i, 0.)

#-

scen_strange = @scenario t_xi = (t_i - 1) s_xi = s_i F_xi = 0. probability = 1.
sp_strange = outcome_model(sp, optimal_decision(sp), scen_strange; optimizer = subproblem_optimizer(sp))
optimize!(sp_strange)

#-

value(sp_strange[:go1])
#-

all_vars_2 = all_vars .* "2"

all_res = map(x -> value.(sp_strange[Symbol(x)]), all_vars_2)
all_res_1 = map(x -> x[1], all_res)

for (v, r) in zip(all_vars_2, all_res_1)
    println("$v : $r")
end

#-

all_res_bc = map(x -> value.(sp_strange[Symbol(x)]), all_vars)
all_res_bc_t_xi = map(x -> length(x) == 1 ? x[1] : x[t_i - 1], all_res_bc)


for (v, r) in zip(all_vars, all_res_bc_t_xi)
    println("$v : $r")
end

