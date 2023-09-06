#=
We focus on three particular pairs of F and scen_freq:
(5000., 48), (25000., 48), (25000., 240)
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
include(joinpath(basepath, "src", "stochastic_optimization_setup.jl"));
#-
# ARGS[1] should be "debug" or run_id
debug = false
@assert length(ARGS)>=1

n_runs = 20 # number of repeats with the same parameters
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
    mkpath(savepath)
end

csv_path = joinpath(savepath, "val_costs.csv")
if !isfile(csv_path)
#m - opt_mode = ["fixed_fg","fixed_fg_inv","OF"]
    open(csv_path, "w") do f
        CSV.write(f,[], writeheader = true, header = ["p", "s_or", "s_val", "m", "cost", "CR"])
    end
end

inv_csv = joinpath(savepath, "inv.csv")
if !isfile(inv_csv)
    #m - opt_mode = ["fixed_fg","fixed_fg_inv","OF"]
    open(inv_csv, "w") do f
        CSV.write(f,[], writeheader = true, 
        header = ["p", "s_or", "m", "u_pv", "u_wind", "u_storage", "u_heat_storage", "u_heatpump"])
    end
end
#-
#-
timesteps = 1:24*365
debug && (timesteps = 1:50)

pv, wind, demand, heatdemand, pars = load_max_boegl(heat = "when2heat");

if debug
    scen_freq, F = 3+pars[:recovery_time], 500.
	n_samples = 5
else
    opt_params = [(500., 96, 60), (5000., 96, 60), (5000., 48, 30)]
    opt_mode = ["OFR","OFOR","OFIOR"]
    run_params = []
    for p in eachindex(opt_params)
        for m_index in 1:3
            for s_or in 1:n_runs
                for s_val in [collect(1:s_or-1);collect(s_or+1:n_val)]
                    push!(run_params, (p, m_index, s_or, s_val))
                end
            end
        end
    end
    i = Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))
    p, m_index, s_or, s_val = run_params[i]
    F, scen_freq, n_samples = opt_params[p]
	println("F = $F")
	println("scen_freq = $(scen_freq)")
    println("Sample #$(s_or)")
    println(opt_mode[m_index])
end
#-
#=
C^B and C^{FG} are already optimized in results/baseline
=#
samp = JSON.parsefile(joinpath(basepath, "samples","scen_freq$(scen_freq)_F_$(F)_n_$(n_samples)", "sample$(s_or).json"))
scens = reconstruct_sample(samp)

samp = JSON.parsefile(joinpath(basepath, "samples","scen_freq$(scen_freq)_F_$(F)_n_$(n_samples)", "sample$(s_val).json"))
scens_val = reconstruct_sample(samp)
samp = nothing;
t_max = 24*365 - 24
delta_t = scen_freq - pars[:recovery_time]
pars[:event_per_scen] = t_max / (delta_t + pars[:recovery_time] + 1)
#-
fg = BSON.load(joinpath(basepath, "results", run_id, "baseline/baseline_$(F).bson"))
inv_fg = fg[:inv]
op_fg = fg[:op]
fg = nothing; #free memory
#-
F_pos = F
F_neg = -F
es_guar = define_energy_system(pv, wind, demand, heatdemand;
p = pars, guaranteed_flex=true, F_pos=F_pos, F_neg=F_neg)
#-
savefile_lock = ReentrantLock()
# First we open the file with validation costs
val_costs = CSV.read(csv_path, DataFrame, header = true)
invests = CSV.read(inv_csv, DataFrame, header = true)
println("Evaluating the decision on sample $s_val")

lock(savefile_lock)
if size(val_costs)[1] != 0
    sel = subset(val_costs, :p => a -> a .== p, :s_or => c -> c .== (s_or), :s_val => s -> s .== s_val, :m => n-> n.==m_index)
else
    sel = []
end
if size(invests)[1] != 0
    sel_inv = subset(val_costs, :p => a -> a .== p, :s_or => c -> c .== (s_or), :m => n-> n.==m_index)
else
    sel_inv = []
end
if size(sel)[1] == 0
    println("Starting the evaluation")
    stime = time()
    filename = "run_$(F)_$(scen_freq)_$(s_or)_$(opt_mode[m_index])"
    if !isfile(filename) # Check if the optimization is already complete
        # useful for reruns if only some failed due to time or memory limits
        savefiles = true
    else 
        savefiles = false
    end
    stime = time()
    fixed_invs = nothing
    fixed_ops = nothing
    if m_index == 1 # fixed FG
        fixed_invs = inv_fg
        fixed_ops = op_fg
    elseif m_index == 2
        fixed_invs = inv_fg
    elseif m_index == 3
        # no need to fix anything
    end
    sp, rt = optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, 
        savefiles = savefiles, savepath = savepath, filename = filename, 
        F_pos = F, F_neg = -F, F_max = F, F_min = F*0.6, resample = true,
        fixed_invs = fixed_invs, fixed_ops = fixed_ops, scens = scens, resample_scens = scens_val)
    println("Runtime in seconds: $(time()-stime)")
	if size(sel_inv)[1] == 0
        invs = get_investments(sp)
        CSV.write(inv_csv, DataFrame(p=p, s_or = s_or, m = m_index, 
        u_pv = invs[:u_pv], u_wind = invs[:u_wind], u_storage = invs[:u_storage], 
        u_heat_storage = invs[:u_heat_storage], u_heatpump = invs[:u_heatpump]), append = true)
    end
    CSV.write(csv_path, DataFrame(p=p, s_or = s_or, s_val = s_val, m = m_index, cost = objective_value(sp), CR = get_servicing_cost(sp)), append = true)
else
    println("$(opt_mode[m_index]) optimization for F = $F, scen_freq = $scen_freq and sample $(s_or) skipped.")
end

unlock(savefile_lock)