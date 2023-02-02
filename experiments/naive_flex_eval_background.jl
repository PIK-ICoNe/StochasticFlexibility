#=
In this file we analyze the model's flexibility potential.
We load an optimized investment decision, find operational schedule for a new sample of determined size.
Then for the model with these fixed investment and operation we find flexibility potential at all time points.
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

#-
using Plots
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));
warm_up()
#-
timesteps = 1:24*365

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:recovery_time] = 12
pars[:sto_ef_ch] = 0.97
pars[:sto_ef_dis] = 0.97

#-
# Analyze availability of flexibility for the background system
Threads.@threads for F in 5000.:5000.:25000.
    es_bkg = define_energy_system(pv, wind, demand, heatdemand; p=pars, override_no_event_per_scen=true, guaranteed_flex=true, F_pos=F, F_neg=F)
    sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
    set_silent(sp_bkg)
    optimize!(sp_bkg)
    cost_bkg = objective_value(sp_bkg)
    bkg_investments = get_investments(sp_bkg)
    bkg_operations = get_operation(sp_bkg)
    CSV.write(joinpath(basepath, "results/naive_flex/bkg", "investments_bkg_$F.csv"), DataFrame(bkg_investments), append = false)
    CSV.write(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg.csv_$F"), Tables.table([cost_bkg]), header = ["Objective value"])
    CSV.write(joinpath(basepath, "results/naive_flex/bkg", "op_bkg.csv_$F"), DataFrame(bkg_operations), append = false)
end
#-
costs = []
append!(costs, cost_bkg)
for F in 5000.:5000.:25000.
    c = CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg.csv_$F"), DataFrame)[!,1]
    append!(costs, c)
end

plot(0.:5000.:25000., costs)

#-
costs_flex = []
append!(costs_flex, costs[2])
for n_samples in 15:5:40
    c = CSV.read(joinpath(basepath, "results/naive_flex/2023-01-17", "costs$(n_samples)_5000.0.csv"), DataFrame)[!,1]
    append!(costs_flex, c)
end

plot([[0.];collect(15:5:40)], costs_flex)

#-
invs = []
for n_samples in 15:5:40
    inv_data = CSV.read(joinpath(basepath, "results/naive_flex/2023-01-17", "inv$(n_samples)_5000.0.csv"), DataFrame)
    inv = Dict(pairs(eachcol(inv_data)))
    inv_data = nothing;
    append!(invs, [inv])
end
