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

include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#
timesteps = 1:24*365
pv, wind, demand, heatdemand = load_basic_example(timesteps);

#=
Next we continue the set up. Our model comes with default parameters,
which we slightly adjust here. We use some arbitrary values to define a dummy heat demand.
=#

pars = copy(default_es_pars)

recovery_time = 12

pars[:recovery_time] = recovery_time
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 550.
pars[:penalty] = 1000000.
#-
#= 
We set up runs with variating scenario frequency and number of scenarios
=#

for n_samples in 10:40:200
    for scen_freq in 24:30*24
        optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, savefiles = true)
    end
end
