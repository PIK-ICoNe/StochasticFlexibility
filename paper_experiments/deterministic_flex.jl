#=
Here we check if the baseline system optimized without guaranteed flexibility still has some level of available flexibility.
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

#-
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"))
using Gurobi
warm_up()

#-
timesteps = 1:24*365

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;

F_range = 0.:10.:100.
read_data = JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_0.0.json"))
    op = read_data["op"]
    inv = read_data["inv"]
#-
# Analyze availability of flexibility for the background system
cost = []
F = 0.
stime = time()
es_bkg = define_energy_system(pv, wind, demand, heatdemand; p=pars, override_no_event_per_scen=true, guaranteed_flex=true, F_pos=F, F_neg=F)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Gurobi.Optimizer)
set_silent(sp_bkg)
fix_investment!(sp_bkg, inv)
optimize!(sp_bkg)
cost_bkg = objective_value(sp_bkg)
#append!(cost, cost_bkg)
runtime = time()-stime
println("Model optimized in $runtime seconds")
println("Cost of model with guaranteed flex $F is $(cost_bkg)")

F_pos, F_neg = naive_flex_potential(get_investments(sp_bkg), get_operation(sp_bkg), pv, wind, pars, timesteps)

model = outcome_model(sp_bkg, optimal_decision(sp_bkg),scen, optimizer = subproblem_optimizer(sp_bkg))
scen = @scenario t_xi = 24 s_xi = 1 F_xi = 100. probability = 1.
outcome = outcome_model(sp_bkg, optimal_decision(sp_bkg), no_flex_pseudo_sampler(); optimizer = subproblem_optimizer(sp_bkg))
set_silent(outcome)
optimize!(outcome)
#-

println(MOI.get.(sp_bkg, MOI.ConstraintConflictStatus(), constraint_by_name(sp_bkg, 2, :strict_flex_in))) 

c = constraint_by_name(sp_bkg, 2, "strict_flex_in")

for t in list_of_constraint_types(sp_bkg,1)
    println(MOI.get.(sp_bkg, MOI.ConstraintConflictStatus(), t)) 
end

compute_conflict!(sp_bkg)
copy_conflict(sp_bkg)
using JuMP
list_of_conflicting_constraints = ConstraintRef[]
for (F, S) in list_of_constraint_types(sp_bkg,2)
    for con in all_constraints(sp_bkg, 2, F, S)
        #if JuMP.get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
        if MOI.get.(sp_bkg, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT
            push!(list_of_conflicting_constraints, con)
        end
    end
end