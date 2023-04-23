#=
We focus on three particular pairs of F and scen_freq:
(5000., 48), (7500., 144), (25000., 240)
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()
using DataFrames

using JSON
using BSON
using Clp
using Statistics;
using StochasticPrograms
using Random

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "data_load.jl"))
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));
#-
# ARGS[1] should be "debug" or run_id
debug = false
@assert length(ARGS)>=1

n_runs = 20 # number of repeats with the same parameters
if length(ARGS) > 1
    n_runs = Base.parse(Int,ARGS[2])
    println(n_runs)
end

if occursin("debug", ARGS[1])
    debug = true
    println("Debug run")
end

!debug && (warm_up())

run_id  = ARGS[1]
stime = time()
savepath = joinpath(basepath, "results", run_id)
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
    params = [(5000., 48), (7500., 144), (25000., 240)]
    i = Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))
    F, scen_freq = params[((i-1)÷n_runs)+1]
	println("F = $F")
	println("scen_freq = $(scen_freq)")
end
#-
#=
C^B and C^{FG} are already optimized in results/baseline
=#
samp = JSON.parsefile(joinpath(basepath, "samples","scen_freq$scen_freq", "sample$((i-1)%n_runs+1).json"))
scens = reconstruct_sample(samp)
samp = nothing;
#-
fg = JSON.parsefile(joinpath(basepath, "results/baseline", "baseline_$(F).json"))
inv_fg = fg["inv"]
op_fg = fg["op"]
fg = nothing; #free memory
#-
#es_no_guar = define_energy_system(pv, wind, demand, heatdemand; p = pars)
F_pos = F
F_neg = -F
es_guar = define_energy_system(pv, wind, demand, heatdemand;
p = pars, guaranteed_flex=true, F_pos=F_pos, F_neg=F_neg)
#-
filename = joinpath(savepath, "run_$(F)_$(scen_freq)_$((i-1)%n_runs+1)_fixed_fg.bson")
if !isfile(filename) # Check if the optimization is already complete
    # useful for reruns if only some failed due to time or memory limits
    stime = time()
    sp = instantiate(es_guar, scens, optimizer = Clp.Optimizer)
    fix_investment!(sp, inv_fg)
    fix_operation!(sp, op_fg, length(timesteps))

    println("Model setup and instantiation performed in $(time() - stime) seconds")
    stime = time()
    set_silent(sp)
    optimize!(sp)
    runtime = time() - stime
    println("Model optimized in $runtime seconds")
    @assert objective_value(sp) != Inf
    stime = time()
    opt_params = Dict((:F_min => 3000., :F_max => F, :t_max_offset => 24, :n_samples => n_samples, :scen_freq => scen_freq, 
                #:F_guar_pos => F_pos, :F_guar_neg => F_neg
                ))
    all_data = get_all_data(sp)
    all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp))))
    #=open(filename, "w") do f
        JSON.print(f,all_data)
    end=#
    bson(filename, all_data)
    println("Data extracted and saved in $(time() - stime) seconds")
    all_data = nothing;
else
    println("Fixed FG optimization for F = $F, scen_freq = $scen_freq and sample $((i-1)%n_runs+1) skipped.")
end
#-
filename = joinpath(savepath, "run_$(F)_$(scen_freq)_$((i-1)%n_runs+1)_fixed_fg_inv.bson")
if !isfile(filename)
    stime = time()
    sp = instantiate(es_guar, scens, optimizer = Clp.Optimizer)
    fix_investment!(sp, inv_fg)
    println("Model setup and instantiation performed in $(time() - stime) seconds")
    set_silent(sp)
    stime = time()
    optimize!(sp)
    runtime = time() - stime
    println("Model optimized in $runtime seconds")
    @assert objective_value(sp) != Inf
    stime = time()
    opt_params = Dict((:F_min => 3000., :F_max => F, :t_max_offset => 24, :n_samples => n_samples, :scen_freq => scen_freq, 
                #:F_guar_pos => F_pos, :F_guar_neg => F_neg
                ))
    all_data = get_all_data(sp)
    all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp))))
    bson(filename, all_data)
    println("Data extracted and saved in $(time() - stime) seconds")
    all_data = nothing;
else 
    println("Fixed investment optimization for F = $F, scen_freq = $scen_freq and sample $((i-1)%n_runs+1) skipped.")

end
#-
filename = joinpath(savepath, "run_$(F)_$(scen_freq)_$((i-1)%n_runs+1)_unfixed_guar_flex.bson")
if !isfile(filename)
    stime = time()
    sp = instantiate(es_guar, scens, optimizer = Clp.Optimizer)
    println("Model setup and instantiation performed in $(time() - stime) seconds")
    stime = time()
    set_silent(sp)
    optimize!(sp)
    runtime = time() - stime
    println("Model optimized in $runtime seconds")
    @assert objective_value(sp) != Inf
    stime = time()
    opt_params = Dict((:F_min => 3000., :F_max => F, :t_max_offset => 24, :n_samples => n_samples, :scen_freq => scen_freq, 
                :F_guar_pos => F_pos, :F_guar_neg => F_neg
                ))
    all_data = get_all_data(sp)
    all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp))))
    bson(filename, all_data)
    println("Data extracted and saved in $(time() - stime) seconds")
else
    println("OF optimization for F = $F, scen_freq = $scen_freq and sample $((i-1)%n_runs+1) skipped.")
end