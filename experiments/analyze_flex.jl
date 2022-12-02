#=
In this file we analyze the model's flexibility potential.
We load an optimized investment decision, find operational schedule for a new sample of determined size.
Then for the model with these fixed investment and operation we find flexibility potential at all time points.
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

#-
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));
warm_up()
#-
timesteps = 1:24*365
n_samples = 15
scen_freq = 48

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps);
pars[:recovery_time] = 12
pars[:sto_ef_ch] = 0.97
pars[:sto_ef_dis] = 0.97

#-
# Analyze cost of flexibility for an optimized system
run_id = "2022-10-25"
invs, mean_invs = load_invs(run_id, n_samples, scen_freq)

savefile_lock = ReentrantLock()
sp, rt = optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, savefile_lock)
#-
F = 7000.

cost_pos, pen_pos = evaluate_cost_penalty(sp, timesteps[end-12],  1, F, pars[:penalty])
cost_neg, pen_neg = evaluate_cost_penalty(sp, timesteps[end-12], -1, F, pars[:penalty])

df_pos = DataFrame(Dict((:t => timesteps[end-12], :F => F, :cost => cost_pos, :penalty => pen_pos)))
df_neg = DataFrame(Dict((:t => timesteps[end-12], :F => F, :cost => cost_neg, :penalty => pen_neg)))
CSV.write(joinpath(basepath, "results", "flex_cost_pos_$(n_samples)_$(scen_freq).csv"), df_pos, writeheader = true)
CSV.write(joinpath(basepath, "results", "flex_cost_neg_$(n_samples)_$(scen_freq).csv"), df_neg, writeheader = true)
savefile_lock = ReentrantLock()
Threads.@threads for i in eachindex(timesteps[1:2])
    println("Timestep $i")
    cost_pos, pen_pos = evaluate_cost_penalty(sp, timesteps[i],  1, F, pars[:penalty])
    cost_neg, pen_neg = evaluate_cost_penalty(sp, timesteps[i], -1, F, pars[:penalty])
    df_pos = DataFrame(Dict((:t => timesteps[i], :F => F, :cost => cost_pos)))
    df_neg = DataFrame(Dict((:t => timesteps[i], :F => F, :cost => cost_neg)))
    lock(savefile_lock)
    CSV.write(joinpath(basepath, "results", "flex_cost_pos_$(n_samples)_$(scen_freq).csv"), df_pos, append = true)
    CSV.write(joinpath(basepath, "results", "flex_cost_neg_$(n_samples)_$(scen_freq).csv"), df_neg, append = true)
    unlock(savefile_lock)
end

#=
df_pos = DataFrame( Dict((:t => timesteps, :flex_cost => cost_pos, :flex_potential => F)))
CSV.write(joinpath(basepath, "results/$run_id", "flex_pos_cost$(n_samples)_$(scen_freq).csv"), df_pos)

df_neg = DataFrame( Dict((:t => timesteps, :flex_cost => cost_neg, :flex_potential => -F)))
CSV.write(joinpath(basepath, "results/$run_id", "flex_neg_cost$(n_samples)_$(scen_freq).csv"), df_neg)
#-
for filename in ["flex_neg_cost$(n_samples)_$(scen_freq).csv", "flex_pos_cost$(n_samples)_$(scen_freq).csv"]
    open(joinpath(basepath, "results/$run_id", filename), "w") do f
        CSV.write(f,[], writeheader=true, header = ["Flexibility_cost", "t", "F"])
    end
end
=#