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
include(joinpath(basepath, "paper_experiments", "stochastic_optimization_setup.jl"));

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
pv, wind, demand, heatdemand, pars = load_max_boegl(heat="when2heat");

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
n_samples = 20 # base n_samples for scen_freq = 48
debug && (scen_freq = 3+pars[:recovery_time])
if debug
    F = 5000.
else
    scen_freq = 48:48:240
    F_range = 5000.:2500.:25000.
    @assert length(scen_freq)*length(F_range) == Base.parse(Int,(ENV["SLURM_ARRAY_TASK_COUNT"]))
    sample_param = []
    for sf in scen_freq
        for F in F_range
            push!(sample_param, (sf,F))
        end
    end
    sf, F = sample_param[Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))]
end
#-
println("scen_freq = $(sf), F = $(F)")

filepath = joinpath(savepath, "run_$(F)_$(sf).bson")
sp, rt = optimize_sp(pv, wind, demand, heatdemand, pars, round(Int, n_samples*sf/48), sf, 
savefiles = true, savepath = filepath, 
F_pos = F, F_neg = -F, F_max = F, F_min = F*0.6, resample = true)
println("Runtime in seconds: $(time()-stime)")
