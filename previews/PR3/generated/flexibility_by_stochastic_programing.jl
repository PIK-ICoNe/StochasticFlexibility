# Everything runs in the Project environment on the basepath

basepath = realpath(joinpath("/home/runner/work/StochasticFlexibility/StochasticFlexibility/manuscripts", ".."))

using Pkg
Pkg.activate(basepath)

using DataFrames
using CSV
using Clp
using Statistics;

using Random
Random.seed!(1);

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

timesteps = 1:(24*365)
pv, wind, demand, heatdemand = load_basic_example(timesteps);

plt = plot(timesteps, pv .* (mean(demand) / mean(pv)), label="PV (unitless)")
plot!(plt, timesteps, wind.* (mean(demand) / mean(wind)), label="Wind (unitless)")
plot!(plt, timesteps, heatdemand, label="Heat Demand")
plot!(plt, timesteps, demand, label="Electric Demand")
plt

pars = copy(default_es_pars)

average_hourly_demand = mean(demand)

recovery_time = 12

pars[:recovery_time] = recovery_time
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 550.
pars[:penalty] = 1000000.;

t_max = length(pv) - 24
F_max = 10000.
F_min = 3000.
delta_t = 3*24 # Flex event every week
pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1);
n = round(Int, 10 * pars[:event_per_scen])

scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = false, override_no_event_per_scen = true)
es_no_reg = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = false)

sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)

set_silent(sp_bkg)

optimize!(sp_bkg)

bkg_decision = optimal_decision(sp_bkg)

objective_value(sp_bkg)

plot_results(sp_bkg, pv, wind, demand)

plot_heat_layer(sp_bkg, heatdemand)

prob_scen = @scenario t_xi = 195 s_xi = 1. F_xi = 20. probability = 1.
evaluate_decision(sp_bkg, bkg_decision, prob_scen)

outcome = outcome_model(sp_bkg, bkg_decision, prob_scen; optimizer = subproblem_optimizer(sp_bkg))
set_silent(outcome)
optimize!(outcome)
termination_status(outcome)

sp_no_reg = instantiate(es_no_reg, scens, optimizer = Clp.Optimizer)
set_silent(sp_no_reg)

evaluate_decision(sp_no_reg, bkg_decision)

es_reg = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = true)

sp_reg = instantiate(es_reg, scens, optimizer = Clp.Optimizer)
set_silent(sp_reg)

evaluate_decision(sp_reg, bkg_decision)

investments_nf = get_investments(sp_bkg)
fix_investment!(sp_no_reg, investments_nf)
fix_investment!(sp_reg, investments_nf)

optimize!(sp_no_reg)
optimize!(sp_reg)

termination_status(sp_no_reg) # Infeasible with F_max 10000 and 1000
termination_status(sp_reg)

reg_no_invest_decision = optimal_decision(sp_reg);

evaluate_decision(sp_reg, no_reg_no_invest_decision)

evaluate_decision(sp_reg, reg_no_invest_decision)

no_invest_relative_cost_of_flex = evaluate_decision(sp_reg, reg_no_invest_decision) / objective_value(sp_bkg) - 1.

unfix_investment!(sp_no_reg, investments_nf)
unfix_investment!(sp_reg, investments_nf)

optimize!(sp_no_reg)
optimize!(sp_reg)

termination_status(sp_no_reg)
termination_status(sp_reg)

no_reg_invest_decision = optimal_decision(sp_no_reg)
reg_invest_decision = optimal_decision(sp_reg);

evaluate_decision(sp_reg, no_reg_invest_decision)

evaluate_decision(sp_reg, reg_invest_decision)

invest_relative_cost_of_flex = evaluate_decision(sp_reg, reg_invest_decision) / objective_value(sp_bkg) - 1.

invest_relative_cost_of_flex / no_invest_relative_cost_of_flex - 1.

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

