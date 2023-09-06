using DataFrames

using JSON
using BSON
using Clp
using Cbc
using Gurobi
using Statistics;
using StochasticPrograms
using Random
using CSV

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

function optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq;
    F_min = 300., F_max = 1000., F_pos=nothing, F_neg=nothing, t_max_offset = 24,
    savefiles=false, savepath=nothing, filename = nothing, fixed_invs=nothing, fixed_ops = nothing,
    scens=nothing, sp_bck=nothing, resample = false, resample_scens = nothing, opt = "Cbc", cap_constraint = "improved_heat")
    number_of_hours = minimum((length(pv),length(wind),length(demand),length(heatdemand)))
    if opt == "Cbc"
        optim = Cbc.Optimizer
    elseif opt == "Gurobi"
        optim = Gurobi.Optimizer
    else
        optim = Clp.Optimizer
    end
    @show optim
    println("Hourly electricity demand: mean = $(mean(demand)), max = $(maximum(demand))")
    println("Hourly heat demand: mean = $(mean(heatdemand)), max = $(maximum(heatdemand))")
    println("Type of constraint used: $cap_constraint")
    println("Hourly PV production at max. investment: mean = $(mean(pv)*pars[:max_pv]), max = $(maximum(pv)*pars[:max_pv])")
    println("Hourly wind production at max. investment: mean = $(mean(wind)*pars[:max_wind]), max = $(maximum(wind)*pars[:max_wind])")

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
    println("Expected hourly flex. request: mean = $(mean([scens[i].data[:F_xi] for i in eachindex(scens)])), max = $(maximum([scens[i].data[:F_xi] for i in eachindex(scens)]))")

    stime = time()

    if isnothing(sp_bck)
        println("Don't have sp provided")
        if !isnothing(F_pos) || !isnothing(F_neg)
            println("Defining system with FG")
            es = define_energy_system(pv, wind, demand, heatdemand;
            p = pars, guaranteed_flex=true, F_pos=F_pos, F_neg=F_neg, cap_constraint = cap_constraint)
        else
            es = define_energy_system(pv, wind, demand, heatdemand; p = pars, cap_constraint = cap_constraint)
        end
    
        sp = instantiate(es, scens, optimizer = optim)
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
    gci, gco = get_penalty_array(sp)
    println("$(length(gci[gci.>0.])) scenarios used gci to balance the request")
    println("$(length(gco[gco.>0.])) scenarios used gco to balance the request")
    #penalty_DF = DataFrame(gci = gci, gco = gco)
    #CSV.write(joinpath(savepath, "penalty_info.csv"), penalty_DF)
    println("Model optimized in $runtime seconds")
    @assert objective_value(sp) != Inf
    first_stage_cost = get_total_investment(sp)/sp.stages[1].parameters[:asset_lifetime] * 365 * 24 / number_of_hours + get_operation_cost(sp)
    println("Cost/kWh = $(first_stage_cost/(sum(demand)+sum(heatdemand)))")
    if resample
        if savefiles
            opt_params = Dict((:F_min => F_min, :F_max => F_max, :t_max_offset => t_max_offset, :n_samples => n_samples, :scen_freq => scen_freq, 
                :F_guar_pos => F_pos, :F_guar_neg => F_neg))
            all_data = get_all_data(sp, rec=false, scen=false)
            all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp),
                                                        :gci_penalty => gci, :gco_penalty => gco)))
            bson(joinpath(savepath, filename*"before_resampling.bson"), all_data)
        end
        stime = time()
        println("Starting resampling")
        inv = get_investments(sp)
        ops = get_operation(sp)
        if isnothing(resample_scens)
            scens_rs = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
        else 
            scens_rs = resample_scens
        end
        sp = instantiate(es, scens_rs, optimizer = optim)
        println("Model setup and instantiation performed in $(time() - stime) seconds")
        fix_investment!(sp, inv)
        fix_operation!(sp, ops, number_of_hours)
        stime = time()
        set_silent(sp)
        optimize!(sp)
        runtime = time() - stime
        println("Model reoptimized in $runtime seconds")
        @assert objective_value(sp) != Inf
        first_stage_cost = get_total_investment(sp)/sp.stages[1].parameters[:asset_lifetime] * 365 * 24 / number_of_hours + get_operation_cost(sp)
        println("After resampling Cost/kWh = $(first_stage_cost/(sum(demand)+sum(heatdemand)))")
        gci, gco = get_penalty_array(sp)
        println("$(length(gci[gci.>0.])) scenarios used gci to balance the request")
        println("$(length(gco[gco.>0.])) scenarios used gco to balance the request")
    end
    if savefiles
        opt_params = Dict((:F_min => F_min, :F_max => F_max, :t_max_offset => t_max_offset, :n_samples => n_samples, :scen_freq => scen_freq, 
            :F_guar_pos => F_pos, :F_guar_neg => F_neg))
        all_data = get_all_data(sp, rec=false, scen=false)
        all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp),
                                                    :gci_penalty => gci, :gco_penalty => gco)))
        bson(joinpath(savepath, filename*".bson"), all_data)
    end
    return sp, runtime
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
