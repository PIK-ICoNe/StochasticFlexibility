## Everything runs in the Project environment on the basepath

#=
In this experiment we analyze investment decision's sensitivity to sampling.
=#
basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

#using Dates

#= Include the file containing all the dependecies and optimization functions.
=#
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));

#-
# ARGS[1] should be "debug" or run_id
debug = false
if ARGS[1] == "debug"
    debug = true
end

!debug && (warm_up())

debug && (timesteps = 1:50)
#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat=true);

#= 
We set up runs with variating guaranteed flexibility and number of scenarios
=#
#-
#run_id = Dates.now()
run_id  = ARGS[1]
stime = time()
mkdir(joinpath(basepath, results, run_id))
savepath = joinpath(basepath, "results", run_id)

n_samples = [collect(15:5:35); collect(40:10:90); collect(100:25:175)]
scen_freq = 96
debug && (n_samples = 1:15)
sample_param = []
for n in n_samples
    for F in 5000.:5000.:25000.
        push!(sample_param, (n,F))
    end
end

n, F = [Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))]

#-
println("n_samples = $(n), F = $(F)")
param_id = "$(n)_$(F)"
savefile_lock = ReentrantLock()
sp, rt = optimize_sp(pv, wind, demand, heatdemand, pars, n, scen_freq, savefile_lock, savefiles = true, savepath = savepath, F_pos = F, F_neg = -F)
println("Runtime in seconds: $(time()-stime)")