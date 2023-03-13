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
warm_up()
#-
timesteps = 1:24*365

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;

n_samples = 20
params = [(5000., 24), (7500., 144), (25000., 240)]
# Fix investment, optimize operation

fixed_inv_path = joinpath(basepath, "results/fixed_inv")
fixed_ops_path = joinpath(basepath, "results/fixed_op")

for path in (fixed_inv_path, fixed_ops_path)
	if !isdir(path)
		mkdir(path)
	end
end
savefile_lock = ReentrantLock()
Threads.@threads for (F, scen_freq) in params
    read_data = JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_$(F).json"))
    op = read_data["op"]
    inv = read_data["inv"]
    
    optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, savefile_lock, fixed_invs = inv, savefiles = true, savepath = fixed_inv_path)
    optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, sf, savefile_lock, fixed_operation = op, savefiles = true, savepath = fixed_ops_path)
end