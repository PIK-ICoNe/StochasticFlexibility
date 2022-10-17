## Everything runs in the Project environment on the basepath

#=
In this experiment we analyze investment decision's sensitivity to sampling.
We variate request samples in two ways.
1. scen_freq defines the frequency of events in the sample. 
For example, scen_freq = 24 means that we expect an event roughly every day.
With scen_freq we calculate events_per_scen.
2. n_samples lets us modify the total sample size. We generate n_samples*events_per_scen requests.

For each pair of (n_samples,scen_freq) in given limits we instantiate a stochastic programming model.
The model is optimized. We save investment decisions, objective value and the sample.
We repeat this n_runs times to get a distribution.
=#
basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

#= Include the file containing all the dependecies and optimization functions.
=#
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#
timesteps = 1:24*365
pv, wind, demand, heatdemand = load_basic_example(timesteps);

#=
Next we continue the set up. Our model comes with default parameters,
which we slightly adjust here.
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

run_id = time()
mkdir("results/$run_id")
savepath = joinpath(basepath, "results/$run_id")

for n_samples in [[1];collect(2:2:10); collect(15:5:25)]
    for scen_freq in 24:24:30*24
        param_id = "$n_samples_$scen_freq"
        savefiles = Dict(([var => joinpath(savepath, string(var)*param_id*".csv")
        for var in [:runtime, :scen, :costs, :inv]]))
            Threads.@threads for n in 1:n_runs
            optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, savefiles = savefiles)
        end
    end
end
