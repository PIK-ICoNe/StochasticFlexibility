#=
Here we optimize the system with different levels of guaranteed flexibility, 
including no flexibility, in absence of flexibility requests.
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

#-
include(joinpath(basepath, "paper_experiments", "stochastic_optimization_setup.jl"))
warm_up()
#-
timesteps = 1:24*365

pv, wind, demand, heatdemand, pars = load_max_boegl(heat = "when2heat");

#-
run_id  = ARGS[1]
stime = time()

savepath = joinpath(basepath,"results", run_id, "baseline")
if !isdir(savepath)
    mkdir(savepath)
end

F_range = [0., 500.]
i = Base.parse(Int,(ENV["SLURM_ARRAY_TASK_ID"]))
i = 2
F = F_range[i]
println("Optimizing for GF $F")
# Analyze availability of flexibility for the background system
#Threads.@threads for F in F_range
stime = time()
filename = "baseline_$(F)"
if !isfile(joinpath(savepath, filename*".bson"))
    es_bkg = define_energy_system(pv, wind, demand, heatdemand; p=pars, regularized =false, override_no_event_per_scen=true, guaranteed_flex=true, F_pos=F, F_neg=-F)
    sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Cbc.Optimizer)
    sp, rt = optimize_sp(pv, wind, demand, heatdemand, pars, 1, 0, 
    savefiles = true, savepath = savepath, filename = "baseline_$(F)",
    sp_bck = sp_bkg, scens = no_flex_pseudo_sampler(),
    F_pos = F, F_neg = -F, F_max = F, F_min = F*0.6, resample = false)
    println("Runtime in seconds: $(time()-stime)")
else
    println("Skipping optimization for $F")
end
#-
