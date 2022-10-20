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

using Plots; pyplotjs()
#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps);

#=
We now do the optimization for the system without any flexibility events. The background system:
=#

# Get background investment decision for comparison

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_event_per_scen = true)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_bkg)
optimize!(sp_bkg)[1,2]
cost_bkg = objective_value(sp_bkg)
bkg_investments = get_investments(sp_bkg)
bkg_operations = get_operation(sp_bkg)
#-
CSV.write(joinpath(basepath, "results/bkg", "investments_bkg.csv"), DataFrame(bkg_investments), header = string.(keys(bkg_investments)), append = false)
CSV.write(joinpath(basepath, "results/bkg", "cost_bkg.csv"), Tables.table([cost_bkg]), header = ["Objective value"])
CSV.write(joinpath(basepath, "results/bkg", "op_bkg.csv"), DataFrame(bkg_operations), header = string.(keys(bkg_operations)), append = false)
#-

cost_pos, pot_pos, cost_neg, pot_neg = analyze_flexibility_potential(sp_bkg, timesteps, pars[:penalty], cost_bkg, tol = 500., maxiter = 10)
#-
df_pos = DataFrame( Dict((:t => timesteps, :flex_cost => cost_pos, :flex_potential => pot_pos)))
CSV.write(joinpath(basepath, "results/bkg", "flex_pos_potential_bkg.csv"), df)

df_neg = DataFrame( Dict((:t => timesteps, :flex_cost => cost_neg, :flex_potential => pot_neg)))
CSV.write(joinpath(basepath, "results/bkg", "flex_neg_potential_bkg.csv"), df)
