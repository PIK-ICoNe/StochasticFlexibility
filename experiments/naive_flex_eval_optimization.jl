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
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat=true);
pars[:recovery_time] = 12
pars[:sto_ef_ch] = 0.97
pars[:sto_ef_dis] = 0.97
#= 
We set up runs with variating scenario frequency and number of scenarios
=#
#-
#run_id = Dates.now()
run_id  = "2023-01-18"
stime = time()
#mkdir("results/$run_id")
savepath = joinpath(basepath, "results/naive_flex/$run_id")

#n_samples = [collect(15:5:35); collect(40:10:80); collect(100:25:250); collect(300:50:500); collect(600:100:1000)]
n_samples = [collect(15:5:35); collect(40:10:90); collect(100:25:175)]
scen_freq = 96
sample_param = []
for n in n_samples
    for F in 5000.:5000.:25000.
        push!(sample_param, (n,F))
    end
end

n, F = sample_param[1]#[Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))]

#-
println("n_samples = $(n), F = $(F)")
param_id = "$(n)_$(F)"
savefiles = Dict(([var => joinpath(savepath, string(var)*param_id*".csv")
for var in [:runtime, :scen, :costs, :inv]]))
instantiate_files(savefiles)
savefile_lock = ReentrantLock()
sp, rt = optimize_sp(pv, wind, demand, heatdemand, pars, n, scen_freq, savefile_lock, savefiles = savefiles, F_pos = F, F_neg = -F)

operations = get_operation(sp)
CSV.write(joinpath(savepath, "operation"*param_id*".csv"), DataFrame(operations), append = false)

#-
println("Runtime in seconds: $(time()-stime)")

