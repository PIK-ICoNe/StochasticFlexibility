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
using CSV

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "data_load.jl"))
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));
#-
# ARGS[1] should be "debug" or run_id
debug = false
@assert length(ARGS)>=1

n_runs = 3 # number of repeats with the same parameters
n_val = 3 # number of validation repeats, reevaluating for n_val new samples
if length(ARGS) > 1
    n_runs = Base.parse(Int,ARGS[2])
    println(n_runs)
    if length(ARGS) > 2
        n_val = Base.parse(Int,ARGS[2])
    end
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

csv_path = joinpath(savepath, "val_costs.csv")
if !isfile(csv_path)
#m - Opt. mode (1 - fixed inv., 2 - free OF opt.)
    open(csv_path, "w") do f
        CSV.write(f,[], writeheader = true, header = ["p", "s_or", "s_val", "m", "cost", "CR"])
    end
end
#-
timesteps = 1:24*365
debug && (timesteps = 1:50)

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;

if debug
    scen_freq, F = 3+pars[:recovery_time], 5000.
	n_samples = 5
    s_or = 1
    s_val = 2
    m = 1
else
    opt_params = [(5000., 48, 25), (7500., 144,70), (25000., 240, 120)]
    i = Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))
    #p = ((i-1)÷n_runs)+1
    #s_or = ((i-1)%n_runs+1)
    run_params = []
    for p in 1:3
        #F, scen_freq, n_samples = params[p]
        for m in 1:2
            for s_or in 1:n_runs
                for s_val in [collect(1:s_or-1);collect(s_or+1:n_val)]
                    push!(run_params, (p, m, s_or, s_val))
                end
            end
        end
    end
    p, m, s_or, s_val = run_params[Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))]
    F, scen_freq, n_samples = opt_params[p]
	println("F = $F")
	println("scen_freq = $(scen_freq)")
    println("Original sample $s_or")
    println("Validation sample $s_val")
end
#-

#es_no_guar = define_energy_system(pv, wind, demand, heatdemand; p = pars)
F_pos = F
F_neg = -F
es_guar = define_energy_system(pv, wind, demand, heatdemand;
p = pars, guaranteed_flex=true, F_pos=F_pos, F_neg=F_neg)
#-
if m == 1
    mode = "fixed_fg_inv"
elseif m == 2
    mode = "unfixed_guar_flex"
end
#modes = Dict(("fixed_fg_inv" => 1, "unfixed_guar_flex" => 2))

    filename = joinpath(savepath, "preval", "run_$(F)_$(scen_freq)_$(s_or)_$mode.bson")
    if isfile(filename) # Check if the pre-validation optimization is already complete
        opt_data = BSON.load(filename)
        inv = opt_data[:inv]
        ops = opt_data[:op]
        opt_data = nothing;
        println("Pre-validation data for sample $(s_or) and mode $mode extracted")
        savefile_lock = ReentrantLock()
        # First we open the file with validation costs
        val_costs = CSV.read(csv_path, DataFrame, header = true)
            println("Evaluating the decision on sample $s_val")
            lock(savefile_lock)
            if size(val_costs)[1] != 0
                sel = subset(val_costs, :p => a -> a .== p, :s_or => c -> c .== (s_or), :s_val => s -> s .== s_val, :m => n-> n.==m)
            else
                sel = []
            end
            if size(sel)[1] == 0
                println("Starting the evaluation")
                stime = time()
                samp = JSON.parsefile(joinpath(basepath, "samples","scen_freq$(scen_freq)_n_$(n_samples)", "sample$(s_or).json"))
                scens_validation = reconstruct_sample(samp)
                samp = nothing;
                sp = instantiate(es_guar, scens_validation, optimizer = Clp.Optimizer)
                fix_investment!(sp, inv)
                fix_operation!(sp, ops, length(timesteps))
                set_silent(sp)
                optimize!(sp)
                runtime = time() - stime
                println("Model optimized in $runtime seconds")
                @assert objective_value(sp) != Inf
                CSV.write(csv_path, DataFrame(p=p, s_or = (s_or), s_val = s_val, m = m, cost = objective_value(sp), CR = get_servicing_cost(sp)), append = true)
            end
            unlock(savefile_lock)
        
    else
        println("Base optimization for sample $(s_or) with mode $mode is missing")
    end

