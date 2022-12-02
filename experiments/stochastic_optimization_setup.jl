using DataFrames
using CSV
using Tables
using Clp
using Statistics;
using StochasticPrograms

using Random

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

function optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq, savefile_lock; F_min = 3000., F_max = 10000., t_max_offset = 24, savefiles = nothing, invs = nothing)
    t_max = minimum((length(pv), length(wind), length(demand), length(heatdemand))) - t_max_offset
    recovery_time = pars[:recovery_time]
    delta_t = scen_freq - recovery_time
    pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1)
    n = round(Int, n_samples * pars[:event_per_scen])
    println("$n total scenarios, with an average of $(pars[:event_per_scen]) events per full time period")
    stime = time()
    scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
    es = define_energy_system(pv, wind, demand, heatdemand; p = pars)
    sp = instantiate(es, scens, optimizer = Clp.Optimizer)
    if maximum(heatdemand) == 0.
        fix!(decision_by_name(sp, 1, :u_heatpump), 0.)
        fix!(decision_by_name(sp, 1, :u_heat_storage), 0.)
    end
    set_silent(sp)
    if !isnothing(invs)
        fix_investment!(sp, invs)
    end
    println("Model setup and instantiation performed in $(time() - stime) seconds")
    stime = time()
    optimize!(sp)
    runtime = time() - stime
    println("Model optimized in $runtime seconds")
    lock(savefile_lock)
    if !isnothing(savefiles)
        if :scen in keys(savefiles)
            CSV.write(savefiles[:scen], DataFrame([s.data for s in scens]), append = true)
        end
        if :inv in keys(savefiles)
            investments = get_investments(sp)
            CSV.write(savefiles[:inv], DataFrame(investments), append = true)
        end
        if :costs in keys(savefiles)
            CSV.write(savefiles[:costs], Tables.table([objective_value(sp)]), append = true)
        end
        if :runtime in keys(savefiles)
            CSV.write(savefiles[:runtime], Tables.table([runtime]), append = true)
        end
    end
    unlock(savefile_lock)
    return sp, runtime
end

function instantiate_files(savefiles)
    header = Dict((:runtime => ["Runtime, seconds"],
    :costs => ["Objective value"],
    :inv => ["u_pv","u_wind","u_heatpump","u_storage","u_heat_storage"],
    :scen => ["t_xi", "s_xi", "F_xi"]
    ))
    for var in keys(savefiles)
        open(savefiles[var], "w") do f
            CSV.write(f,[], writeheader=true, header=header[var])
        end
    end
end

function warm_up()
    stime = time()
    pv, wind, demand, heatdemand = load_basic_example(1:24);
    pars = copy(default_es_pars)
    es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_event_per_scen = true)
    sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
    set_silent(sp_bkg)
    optimize!(sp_bkg)
    println("Warm up performed in $(time() - stime) seconds")
end

function get_background_model(pv, wind, demand, heatdemand, pars)
    stime = time()
    bkg = CSV.read(joinpath(basepath, "results/bkg", "investments_bkg.csv"), DataFrame)
    bkg_inv = Dict(pairs(eachcol(bkg)))
    bkg_inv = Dict((
        [v => mean(bkg_inv[v]) for v in keys(bkg_inv)]
    ))

    bkg = CSV.read(joinpath(basepath, "results/bkg", "op_bkg.csv"), DataFrame)
    bkg_op = Dict(pairs(eachcol(bkg)))

    bkg = nothing;
    es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, override_no_event_per_scen = true)
    sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
    fix_investment!(sp_bkg, bkg_inv)
    fix_operation!(sp_bkg, bkg_op, length(timesteps))
    set_silent(sp_bkg)
    optimize!(sp_bkg)
    println("Background model instntiated and optimized in $(time()-stime) seconds")
    return sp_bkg
end