basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
using JSON
#-
include(joinpath(basepath,"src", "plot_utils.jl"))
#a = JSON.parsefile("results/debug/run_2_15_5000.0_-5000.0.json")
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat=true);

a = JSON.parsefile("results/run_02_16/run_100_96_25000.0_-25000.0.json")

plot_results(a, pv, wind, demand)