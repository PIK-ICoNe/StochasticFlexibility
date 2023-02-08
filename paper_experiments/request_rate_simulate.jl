## Everything runs in the Project environment on the basepath

#=
In this experiment we analyze investment decision's sensitivity to sampling.
=#

debug = false
if ARGS[1] == "debug"
    debug = true
end
    

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using Dates

#= Include the file containing all the dependecies and optimization functions.
=#
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));

#-
warm_up()
#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#
timesteps = 1:24*365

debug && timetsteps = 1:10

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat=true);

#= 
We set up runs with variating guaranteed flexibility and number of scenarios
=#
#-
run_id = Dates.now()
#run_id  = "2023-01-18"
stime = time()
mkdir(joinpath(basepath, results, run_id))
savepath = joinpath(basepath, results, run_id)

n_samples = [collect(15:5:35); collect(40:10:100)]
scen_freq = 24*4:48:24*14

debug && n_samples = [2, 4]
debug && scen_freq = [3, 4]
# More work for the debug run.

sample_param = []
for n in n_samples
    for s in scen_freq
        push!(sample_param, (n,s))
    end
end

F = 15000. # choose later
n, s = [Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))]

#-
println("n_samples = $(n), scen_freq = $(s)")
param_id = "$(n)_$(s)"
№mkdir(joinpath(savepath, param_id))
savefiles = Dict(([var => joinpath(savepath, string(var)*param_id*".csv")
    for var in [:runtime, :scen, :costs, :inv]]))
instantiate_files(savefiles)
savefile_lock = ReentrantLock()
sp, rt = optimize_sp(pv, wind, demand, heatdemand, pars, n, s, savefile_lock, savefiles = savefiles, F_pos = F, F_neg = -F)

operations = get_operation(sp)
CSV.write(joinpath(savepath, "operation"*param_id*".csv"), DataFrame(operations), append = false)

#-
println("Runtime in seconds: $(time()-stime)")