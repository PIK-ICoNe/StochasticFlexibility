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
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"))
warm_up()
#-
timesteps = 1:24*365

pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;

#-
savepath = joinpath(basepath,"results","baseline")
if !isdir(savepath)
    mkdir(savepath)
end

F_range = 0.:2500.:25000.
# Analyze availability of flexibility for the background system
Threads.@threads for F in F_range
    stime = time()
    es_bkg = define_energy_system(pv, wind, demand, heatdemand; p=pars, override_no_event_per_scen=true, guaranteed_flex=true, F_pos=F, F_neg=-F)
    sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
    set_silent(sp_bkg)
    optimize!(sp_bkg)
    cost_bkg = objective_value(sp_bkg)
    bkg_investments = get_investments(sp_bkg)
    bkg_operations = get_operation(sp_bkg)
    runtime = time()-stime
    println("Model optimized in $runtime seconds")
    #CSV.write(joinpath(basepath, "results/bkg", "investments_bkg_$F.csv"), DataFrame(bkg_investments), append = false)
    #CSV.write(joinpath(basepath, "results/bkg", "cost_bkg_$F.csv"), Tables.table([cost_bkg]), header = ["Objective value"], append = false)
    #CSV.write(joinpath(basepath, "results/bkg", "op_bkg_$F.csv"), DataFrame(bkg_operations), append = false)
    opt_params = Dict((:F_guar_pos => F, :F_guar_neg => -F))
        all_data = get_all_data(sp_bkg)
        all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp_bkg))))
        open(joinpath(savepath, "baseline_$(F).json"), "w") do f
            JSON.print(f,all_data)
        end
end
#-
