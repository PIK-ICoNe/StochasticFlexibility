#=
In this file we analyze the background model. We initialize energy system without flexibility requests.
We find investment decisions and model cost without knowledge of flexibility. 
Then we find model's potential to provide flexibility with background investment decision.
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using DataFrames
using CSV
using Tables
using Clp
using Statistics;
using StochasticPrograms

using Random

#-
include(joinpath(basepath, "src", "sp_model.jl"))
# include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"))
include(joinpath(basepath, "src", "data_load.jl"));
#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps);

pars[:sto_ef_ch] = 0.97
pars[:sto_ef_dis] = 0.97
pars[:recovery_time] = 12

#=
We now do the optimization for the system without any flexibility events. The background system:
=#

# Get background investment decision for comparison

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_event_per_scen = true)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_bkg)
optimize!(sp_bkg)
cost_bkg = objective_value(sp_bkg)
bkg_investments = get_investments(sp_bkg)
bkg_operations = get_operation(sp_bkg)
#-
CSV.write(joinpath(basepath, "results/bkg", "investments_bkg.csv"), DataFrame(bkg_investments), append = false)
CSV.write(joinpath(basepath, "results/bkg", "cost_bkg.csv"), Tables.table([cost_bkg]), header = ["Objective value"])
CSV.write(joinpath(basepath, "results/bkg", "op_bkg.csv"), DataFrame(bkg_operations), append = false)
#-

cost_pos_bkg = zeros(length(timesteps));
cost_neg_bkg = zeros(length(timesteps));

#-
# Analyze cost of flexibility for a fixed request
F = 7000. # 70% of F_max
cost_pos_bkg, cost_neg_bkg = evaluate_cost(sp_bkg, optimal_decision(sp_bkg), timesteps, F)

df_pos = DataFrame( Dict((:t => timesteps, :flex_cost => cost_pos_bkg, :flex_potential => F)))
CSV.write(joinpath(basepath, "results/bkg", "flex_pos_cost_bkg.csv"), df_pos)

df_neg = DataFrame( Dict((:t => timesteps, :flex_cost => cost_neg_bkg, :flex_potential => -F)))
CSV.write(joinpath(basepath, "results/bkg", "flex_neg_cost_bkg.csv"), df_neg)