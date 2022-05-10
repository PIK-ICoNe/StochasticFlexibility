# Everything runs in the Project environment on the basepath

basepath = realpath(joinpath("/home/runner/work/StochasticFlexibility/StochasticFlexibility/manuscripts", ".."))

using Pkg
Pkg.activate(basepath)

using DataFrames
using CSV
using Clp
using Statistics

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"))

offset = 6531
timesteps = 1:168

data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)

pv = data[timesteps .+ offset, 3]
wind = data[timesteps .+ offset, 4]
demand = data[timesteps .+ offset, 2]

data = nothing # Free the memory

#

pars = copy(default_es_pars)

average_hourly_demand = mean(demand)

pars[:recovery_time] = 24
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 800.
pars[:c_sto_op] = 0.00001

heatdemand = copy(demand)./300.
heatdemand0 = zeros(length(demand))

#

es = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true)

#

sp_no_flex = instantiate(es, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)

optimize!(sp_no_flex)

no_flex_decision = optimal_decision(sp_no_flex)

objective_value(sp_no_flex)

#

cost_pos, pot_pos, cost_neg, pot_neg = analyze_flexibility_potential(sp_no_flex, 20:50)

plot_flexibility(20:50, cost_pos, pot_pos, cost_neg, pot_neg)

#

n = 100
F_max = average_hourly_demand * 0.1 # Have unaticipated demand equal to 10% of our typical demand
t_max = length(pv) - es.parameters[2].defaults[:recovery_time]

scens = simple_flex_sampler(n, F_max, t_max)

#

sp_flex = instantiate(es, scens, optimizer = Clp.Optimizer)

evaluate_decision(sp_flex, no_flex_decision)

#

es_reg = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = false)

sp_reg_flex = instantiate(es_reg, scens, optimizer = Clp.Optimizer)


#

evaluate_decision(sp_reg_flex, no_flex_decision)

#

investments_nf = get_investments(sp_no_flex)
fix_investment!(sp_flex, investments_nf)
fix_investment!(sp_reg_flex, investments_nf)

optimize!(sp_flex)
optimize!(sp_reg_flex)

flex_no_invest_decision = optimal_decision(sp_flex)
reg_flex_no_invest_decision = optimal_decision(sp_reg_flex)

#

evaluate_decision(sp_reg_flex, flex_no_invest_decision)

#
evaluate_decision(sp_reg_flex, reg_flex_no_invest_decision)

#

relative_flex_cost = evaluate_decision(sp_reg_flex, flex_no_invest_decision) / objective_value(sp_no_flex) - 1.

#

cost_pos_flex, pot_pos_flex, cost_neg_flex, pot_neg_flex = analyze_flexibility_potential(sp_flex, 20:50)

plot_flexibility(20:50, cost_pos_flex, pot_pos_flex, cost_neg_flex, pot_neg_flex)

#

plt_av = plot();
flexibility_availability!(plt_av, pot_pos, label = "positive flexbility unaware");
flexibility_availability!(plt_av, pot_pos_flex, label = "positive flexibility aware");
flexibility_availability!(plt_av, pot_neg, label = "negative flexibility unaware");
flexibility_availability!(plt_av, pot_neg_flex, label = "negative flexibility aware");
display(plt_av)

#

unfix_investment!(sp_flex, investments_nf)
unfix_investment!(sp_reg_flex, investments_nf)

optimize!(sp_flex)
optimize!(sp_reg_flex)

flex_invest_decision = optimal_decision(sp_flex)
reg_flex_invest_decision = optimal_decision(sp_reg_flex)

evaluate_decision(sp_reg_flex, flex_invest_decision)

#
evaluate_decision(sp_reg_flex, reg_flex_invest_decision)

#

relative_flex_cost_inv = evaluate_decision(sp_reg_flex, flex_invest_decision) / objective_value(sp_no_flex) - 1.

#

relative_flex_cost_inv / relative_flex_cost - 1.

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

