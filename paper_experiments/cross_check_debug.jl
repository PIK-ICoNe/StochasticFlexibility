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
debug = false

!debug && (warm_up())

stime = time();
savepath = joinpath(basepath, "results", "debug")
if !isdir(savepath)
    mkdir(savepath)
end

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
    params = [(7500., 144), (25000., 240)]
    F, scen_freq = params[1]#[Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))]
	println("F = $F")
	println("scen_freq = $(scen_freq)")
end
#-
# Fix investment, optimize operation

fixed_inv_path = joinpath(savepath, "fixed_inv")

if !isdir(fixed_inv_path)
    mkdir(fixed_inv_path)
end

#-
savefile_lock = ReentrantLock()

t_max_offset = 24
t_max = minimum((length(pv), length(wind), length(demand), length(heatdemand))) - t_max_offset
recovery_time = pars[:recovery_time]
@assert scen_freq > recovery_time
delta_t = scen_freq - recovery_time
pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1)
n = round(Int, n_samples * pars[:event_per_scen])
println("$n total scenarios, with an average of $(pars[:event_per_scen]) events per full time period")
scens = poisson_events_with_offset(n, delta_t, recovery_time, F, t_max, F_min = 3000.)

read_data = JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_$(F).json"))
invs = read_data["inv"]
@show invs
read_data = nothing;
sp_bkg, = optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, savefile_lock, fixed_invs = invs, savefiles = false, savepath = fixed_inv_path, F_max = F, scens = scens)

optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, savefile_lock, savefiles = false, savepath = savepath, F_max = F, scens = scens, sp_bck = sp_bkg)
