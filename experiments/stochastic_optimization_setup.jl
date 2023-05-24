using DataFrames

using JSON
using Clp
using Statistics;
using StochasticPrograms
using Random

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

function optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq;
    F_min = 3000., F_max = 10000., F_pos=nothing, F_neg=nothing, t_max_offset = 24,
    savefiles=false, savepath=nothing, fixed_invs=nothing, fixed_ops = nothing,
    scens=nothing, sp_bck=nothing, resample = false, resample_scens = nothing)
    number_of_hours = minimum((length(pv),length(wind),length(demand),length(heatdemand)))
    if savefiles
        @assert !isnothing(savepath)
        println("Ready to save files")
    end
    if isnothing(scens)
        println("Sampling")
        t_max = number_of_hours - t_max_offset
        recovery_time = pars[:recovery_time]
        @assert scen_freq > recovery_time
        delta_t = scen_freq - recovery_time
        pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1)
        n = round(Int, n_samples * pars[:event_per_scen])
        println("$n total scenarios, with an average of $(pars[:event_per_scen]) events per full time period")
        scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
    end
    stime = time()

    if isnothing(sp_bck)
        println("Don't have sp provided")
        if !isnothing(F_pos) || !isnothing(F_neg)
            println("Defining system with FG")
            es = define_energy_system(pv, wind, demand, heatdemand;
            p = pars, guaranteed_flex=true, F_pos=F_pos, F_neg=F_neg)
        else
            es = define_energy_system(pv, wind, demand, heatdemand; p = pars)
        end
    
        sp = instantiate(es, scens, optimizer = Clp.Optimizer)
    else
        sp = sp_bck
        println("Warning: ignoring everything but fixed_inv and using provided stochastic program")
    end
    # Prevent investment in the heat components if there is no heat demand
    if maximum(heatdemand) == 0.
        fix(decision_by_name(sp, 1, "u_heatpump"), 0.)
        fix(decision_by_name(sp, 1, "u_heat_storage"), 0.)
    end
    set_silent(sp)
    if !isnothing(fixed_invs)
        fix_investment!(sp, fixed_invs)
    else
        unfix_investment!(sp)
    end
    if !isnothing(fixed_ops)
        fix_operation!(sp, fixed_ops, number_of_hours)
    else
        #unfix_operation!(sp, number_of_hours)
    end

    println("Model setup and instantiation performed in $(time() - stime) seconds")
    stime = time()
    optimize!(sp)
    runtime = time() - stime
    println("Model optimized in $runtime seconds")
    @assert objective_value(sp) != Inf
    if resample
        stime = time()
        println("Starting resampling")
        inv = get_investments(sp)
        ops = get_operation(sp)
        if isnothing(resample_scens)
            scens_rs = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
        else 
            scens_rs = resample_scens
        end
        sp = instantiate(es, scens_rs, optimizer = Clp.Optimizer)
        println("Model setup and instantiation performed in $(time() - stime) seconds")
        fix_investment!(sp, inv)
        fix_operation!(sp, ops, number_of_hours)
        stime = time()
        set_silent(sp)
        optimize!(sp)
        runtime = time() - stime
        println("Model reoptimized in $runtime seconds")
        @assert objective_value(sp) != Inf
    end
    if savefiles
        opt_params = Dict((:F_min => F_min, :F_max => F_max, :t_max_offset => t_max_offset, :n_samples => n_samples, :scen_freq => scen_freq, 
            :F_guar_pos => F_pos, :F_guar_neg => F_neg))
        all_data = get_all_data(sp, rec=false, scen=false)
        all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp))))
        if occursin("json", savepath)
            open(savepath, "w") do f
                JSON.print(f,all_data)
            end# joinpath(savepath, "run_$(n_samples)_$(scen_freq)_$(F_pos)_$(F_neg).json")
        elseif occursin("bson", savepath)
            bson(savepath, all_data)
        end
    end
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

"""
Get optimized background model with fixed investment decision read from file.
"""
# TODO revise, add asserts or provide savepath
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