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

#=
Here we optimize operational schedule and 
investment decision separartely, 
fixing one of them to that of the baseline scenario.

We focus on three particular pairs of F and scen_freq:
(5000., 24), (7500., 144), (25000., 240)
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

#-
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"))

#-
# ARGS[1] should be "debug" or run_id
debug = true
@assert length(ARGS)>=1

if occursin("debug", ARGS[1])
    debug = true
    println("Debug run")
end

!debug && (warm_up())

run_id  = ARGS[1]
stime = time()
if !isdir(joinpath(basepath, "results", run_id))
    mkdir(joinpath(basepath, "results", run_id))
end
savepath = joinpath(basepath, "results", run_id)
#-
timesteps = 1:24*365
debug && (timesteps = 1:50)

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;

if debug
    scen_freq, F = 3+pars[:recovery_time], 5000.
	n_samples = 5
else
    n_samples = 20
    F_range = 5000.:5000.:25000.
    params = [(5000., 24), (7500., 144), (25000., 240)]
    F, scen_freq = sample_param[Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))]
end
#-
# Fix investment, optimize operation

fixed_inv_path = joinpath(savepath, "fixed_inv")
fixed_ops_path = joinpath(savepath, "fixed_op")

for path in (fixed_inv_path, fixed_ops_path)
	if !isdir(path)
		mkdir(path)
	end
end
savefile_lock = ReentrantLock()

read_data = JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_$(F).json"))
op = read_data["op"]
inv = read_data["inv"]
    
sp, = optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, savefile_lock, fixed_invs = inv, savefiles = false)


optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, sf, savefile_lock, fixed_operation = op, savefiles = true, savepath = fixed_ops_path)

opt_params = Dict((:t_max_offset => 24, :n_samples => n_samples, :scen_freq => scen_freq, 
:F_guar_pos => F_pos, :F_guar_neg => F_neg))
all_data = get_all_data(sp_bkg)
all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp))))
open(joinpath(savepath, "run_$(n_samples)_$(scen_freq)_$(F_pos)_$(F_neg).json"), "w") do f
JSON.print(f,all_data)
end