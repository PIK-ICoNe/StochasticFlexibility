## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using DataFrames
using CSV
using Clp
using Statistics;

using Random
Random.seed!(1);

#-

#=
# Flexibility analysis and optimization

We take an energy system and analyze and optimize its flexibility potential.

The model is defined in the file sp_model.jl, analysis and utility functions are in the other two files.

The parametrization is not completely done yet, and running the analysis for a whole year is very time intensive.
The main optimization runs fine for a whole year, but analyzing the flexibility potential for many hours
is computationally intense.

We need to either start caching more computations, or find a way to speed up the evaluation further.

This file gives a conceptual demonstration of what we do on a two week horizon.
=#

#- 

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#

pv, wind, demand, heatdemand = load_basic_example(1:2*365);

#-

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
n = round(Int, 30 * pars[:event_per_scen])

scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)

#-
scens_plot = plot_scenario_distribution(scens)
savefig(scens_plot, "scenplot.png")
scens_plot_sign = plot_scenario_distribution_by_sign(scens)
savefig(scens_plot_sign, "scenplot_sign.png")

#-
es = define_energy_system(pv, wind, demand, heatdemand; p = pars)
#-
es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_event_per_scen = true)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_bkg)

optimize!(sp_bkg)

bkg_decision = optimal_decision(sp_bkg)
bkg_cost = objective_value(sp_bkg)
bkg_investments = get_investments(sp_bkg)
#-

sp_op = instantiate(es, scens, optimizer = Clp.Optimizer)
set_silent(sp_op)

fix_investment!(sp_op, bkg_investments)
optimize!(sp_op)

op_decision = optimal_decision(sp_op)
op_cost = objective_value(sp_op)
println("Cost of system with optimized operation relative to background: $(op_cost/bkg_cost)")
#-

sp = instantiate(es, scens, optimizer = Clp.Optimizer)
set_silent(sp)

optimize!(sp)

decision = optimal_decision(sp)
cost = objective_value(sp)
investments = get_investments(sp)
println("Cost of optimized system relative to background: $(cost/bkg_cost)")

#-
println("Additional cost of flexibility above background: $(op_cost - bkg_cost)")
println("Additional cost of flexibility above background with flex optimized investment: $(cost - bkg_cost)")
println("Relative reduction of flexibility cost due to flex optimized investment: $(1 - (cost - bkg_cost)/(op_cost - bkg_cost))")

#-
scens_resampled = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)

sp_resampled = instantiate(es, scens, optimizer = Clp.Optimizer)
fix_investment!(sp_resampled, investments)
set_silent(sp_resampled)
optimize!(sp_resampled)
cost_resampled = objective_value(sp_resampled)
println("Cost of system with optimized investment on a new sample: $(cost_resampled/bkg_cost)")
#-

sankey_results(sp_bkg, pv, wind, demand, timesteps)
