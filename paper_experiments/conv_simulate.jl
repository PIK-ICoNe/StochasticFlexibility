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
@assert length(ARGS)>=1

if occursin("debug", ARGS[1])
    debug = true
    println("Debug run")
end

!debug && (warm_up())


#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#
#-
timesteps = 1:24*365
debug && (timesteps = 1:50)
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat=true);
pars[:inv_budget] = 10^10;
#= 
We set up runs with variating guaranteed flexibility and number of scenarios
=#
#-
run_id  = ARGS[1]
stime = time()
if !isdir(joinpath(basepath, "results", run_id))
    mkdir(joinpath(basepath, "results", run_id))
end
savepath = joinpath(basepath, "results", run_id)
scen_freq = 96
debug && (scen_freq = 3+pars[:recovery_time])
if debug
    n, F = 2, 5000.
else
    n_samples = [collect(15:5:35); collect(40:10:90); collect(100:25:150)]
    F_range = 5000.:5000.:25000.
    @assert length(n_samples)*length(F_range) == Base.parse(Int,(ENV["SLURM_ARRAY_TASK_COUNT"]))
    sample_param = []
    for n in n_samples
        for F in F_range
            push!(sample_param, (n,F))
        end
    end
    n, F = sample_param[Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))]
end
#-
println("n_samples = $(n), F = $(F)")
#param_id = "$(n)_$(F)"
filepath = joinpath(savepath, "conv_run_$(F)_$(scen_freq)_$(n).bson")
sp, rt = optimize_sp(pv, wind, demand, heatdemand, pars, n, scen_freq, 
savefiles = true, savepath = filepath, 
F_pos = F, F_neg = -F, F_max = F, F_min = F*0.6, resample = true)
println("Runtime in seconds: $(time()-stime)")