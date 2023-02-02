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

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps);
pars[:recovery_time] = 12
pars[:sto_ef_ch] = 0.97
pars[:sto_ef_dis] = 0.97

#-
# Analyze availability of flexibility for the background system
inv_data = CSV.read(joinpath(basepath, "results/bkg", "investments_bkg.csv"), DataFrame)
inv = Dict(pairs(eachcol(inv_data)))
inv_data = nothing;
inv = Dict((
        [v => mean(inv[v]) for v in keys(inv)]
    ))
#-

op_data = CSV.read(joinpath(basepath, "results/bkg", "op_bkg.csv"), DataFrame)
op = Dict(pairs(eachcol(op_data)))
op_data = nothing;

#-

F_pos, F_neg = naive_flex_potential(inv, op, pv, wind, pars, timesteps)

plot(timesteps[2:end-12],F_pos[2:end])
plot(F_neg)