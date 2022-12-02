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
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps);
pars[:recovery_time] = 12
pars[:sto_ef_ch] = 0.97
pars[:sto_ef_dis] = 0.97
#= 
We set up runs with variating scenario frequency and number of scenarios
=#

run_id = "2022-11-14"
stime = time()
#mkdir("results/$run_id")
savepath = joinpath(basepath, "results/$run_id")

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_event_per_scen = true)
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_bkg)
optimize!(sp_bkg)
cost_bkg = objective_value(sp_bkg)
bkg_investments = get_investments(sp_bkg)
bkg_operations = get_operation(sp_bkg)
#-
CSV.write(joinpath(basepath, "results/$(run_id)", "investments_bkg.csv"), DataFrame(bkg_investments), append = false)
CSV.write(joinpath(basepath, "results/$(run_id)", "cost_bkg.csv"), Tables.table([cost_bkg]), header = ["Objective value"])
CSV.write(joinpath(basepath, "results/$(run_id)", "op_bkg.csv"), DataFrame(bkg_operations), append = false)

sp_bkg = nothing;

n_samples = 50
scen_freq = 72
param_id = "$(n_samples)_$(scen_freq)"
t_max = minimum((length(pv), length(wind), length(demand), length(heatdemand))) - 24
recovery_time = pars[:recovery_time]
delta_t = scen_freq - recovery_time
pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1)
n = round(Int, n_samples * pars[:event_per_scen])
es = define_energy_system(pv, wind, demand, heatdemand; p = pars)
scens = poisson_events_with_offset(n, delta_t, recovery_time, 10000., t_max, F_min = 3000.)
sp = instantiate(es, scens, optimizer = Clp.Optimizer)
set_silent(sp)
optimize!(sp)
cost = objective_value(sp)
investments = get_investments(sp)
operations = get_operation(sp)
CSV.write(joinpath(basepath, "results/$(run_id)", "investments.csv"), DataFrame(bkg_investments), append = false)
CSV.write(joinpath(basepath, "results/$(run_id)", "cost.csv"), Tables.table([cost_bkg]), header = ["Objective value"])
CSV.write(joinpath(basepath, "results/$(run_id)", "op.csv"), DataFrame(bkg_operations), append = false)
#-
fix_investments!(sp, bkg_investments)
optimize!(sp)
cost = objective_value(sp)
investments = get_investments(sp)
operations = get_operation(sp)
CSV.write(joinpath(basepath, "results/$(run_id)", "investments_fix.csv"), DataFrame(bkg_investments), append = false)
CSV.write(joinpath(basepath, "results/$(run_id)", "cost_fix.csv"), Tables.table([cost_bkg]), header = ["Objective value"])
CSV.write(joinpath(basepath, "results/$(run_id)", "op_fix.csv"), DataFrame(bkg_operations), append = false)

n_samples_array = 70:10:150

n_samples
for n_samples in [collect(2:2:10); collect(15:5:20)]
    for scen_freq in 24:24:30*24
        param_id = "$(n_samples)_$(scen_freq)"
        savefiles = Dict(([var => joinpath(savepath, string(var)*param_id*".csv")
        for var in [:runtime, :scen, :costs, :inv]]))
        instantiate_files(savefiles)
        savefile_lock = ReentrantLock()
        Threads.@threads for n in 1:n_runs
            optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, savefile_lock, savefiles = savefiles)
        end
    end
end
#-
println(time()-stime)