using DataFrames
using CSV
using Tables
using Clp
using Statistics;
using StochasticPrograms

using Random

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

function optimize_sp(pv, wind, demand, heatdemand, pars, n_samples, scen_freq; F_min = 3000., F_max = 10000., t_max_offset = 24, savefiles = nothing)
    t_max = minimum((length(pv), length(wind), length(demand), length(heatdemand))) - t_max_offset
    delta_t = scen_freq - recovery_time
    pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1)
    n = round(Int, n_samples * pars[:event_per_scen])
    println("$n total scenarios, with an average of $(pars[:event_per_scen]) events per full time period")
    stime = time()
    scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
    es = define_energy_system(pv, wind, demand, heatdemand; p = pars)
    sp = instantiate(es, scens, optimizer = Clp.Optimizer)
    set_silent(sp)
    println("Model setup and instantiation performed in $(time() - stime) seconds")
    stime = time()
    optimize!(sp)
    runtime = time() - stime
    println("Model optimized in $runtime seconds")
    investments = get_investments(sp)
    if !isnothing(savefiles)
        if :scen in keys(savefiles)
            CSV.write(savefiles[:scen], DataFrame([s.data for s in scens]), header = string.(keys(scens[1].data)), append = true)
        end
        if :inv in keys(savefiles)
            CSV.write(savefiles[:inv], DataFrame(investments), header = string.(keys(investments)), append = true)
        end
        if :costs in keys(savefiles)
            CSV.write(savefiles[:costs], Tables.table([objective_value(sp)]), header = ["Objective value"], append = true)
        end
        if :runtime in keys(savefiles)
            CSV.write(savefiles[:runtime], Tables.table([runtime]), header = ["Runtime"], append = true)
        end
    end
    return sp
end