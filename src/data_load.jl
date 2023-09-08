using DataFrames
using CSV

## TODO make sure somehow that default_es_pars is known

# pass heat = "district" for district heating, heat = "when2heat" for germany
# "district" is default
function load_max_boegl(;
    basepath=joinpath(@__DIR__,".."),
    time_steps=nothing,
    offset=0,
    heat="district",
    scale_factor=1.0
)

    pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "pv_Halle18.csv"), DataFrame, header=false)
    wind_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "wind_Karholz.csv"), DataFrame, header=false)
    demand_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "demand_Industriepark.csv"), DataFrame, header=false)

    if isnothing(time_steps)
        timesteps = 1:length(pv_data[:, 1])
    else
        timesteps = time_steps
    end

    if heat == "when2heat"
        heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)
        heatdemand = heatdemand_data[timesteps.+offset, 1] ./ (50.0 * scale_factor)
    elseif heat == "district"
        heatdemand_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_heatdemand_apartment_block(KFW40)_300km2_12300MWhperyear.csv"), DataFrame, header=false)
        heatdemand = heatdemand_data[timesteps.+offset, 1] ./ scale_factor
    else
        heatdemand = zeros(length(timesteps))
    end
    pv = pv_data[timesteps.+offset, 1]
    wind = wind_data[timesteps.+offset, 1]
    demand = demand_data[timesteps.+offset, 1]

    pv_data = nothing
    wind_data = nothing
    demand_data = nothing # Free the memory

    pars = copy(default_es_pars)

    pars[:recovery_time] = 12
    pars[:c_storage] = 300.0
    pars[:c_pv] = 800.0
    pars[:c_wind] = 1150.0
    pars[:c_in] = 0.165
    pars[:c_out] = 0.02
    pars[:asset_lifetime] = 10 # increased
    pars[:inv_budget] = 10000000.0
    pars[:sto_ef_ch] = 0.97
    pars[:sto_ef_dis] = 0.97
    pars[:feedincap] = 1e7
    pars[:max_pv] = 2 * 10^4 # increased by factor of 2
    pars[:max_wind] = 10^4
    pars[:inv_budget] = 10^10
    println("Hourly electricity demand: min = $(minimum(demand)),  mean = $(mean(demand)), max = $(maximum(demand))")
    println("Hourly heat demand: min = $(minimum(heatdemand)),  mean = $(mean(heatdemand)), max = $(maximum(heatdemand))")
    println("Hourly PV production at max. investment: mean = $(mean(pv)*pars[:max_pv]), max = $(maximum(pv)*pars[:max_pv])")
    println("Hourly wind production at max. investment: mean = $(mean(wind)*pars[:max_wind]), max = $(maximum(wind)*pars[:max_wind])")
    return pv, wind, demand, heatdemand, pars
end

function load_max_boegl_district_heating(
    timesteps;
    basepath=joinpath(@__DIR__, ".."),
    offset=0,
    heat=false,
    heat_timeseries="district"
)
    pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "pv_Halle18.csv"), DataFrame, header=false)
    wind_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "wind_Karholz.csv"), DataFrame, header=false)
    demand_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "demand_Industriepark.csv"), DataFrame, header=false)
    if heat
        heatdemand_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_heatdemand_apartment_block(KFW40)_300km2_12300MWhperyear.csv"), DataFrame, header=false)
        heatdemand = heatdemand_data[timesteps.+offset, 1]
    else
        heatdemand = zeros(length(timesteps))
    end
    pv = pv_data[timesteps.+offset, 1]
    wind = wind_data[timesteps.+offset, 1]
    demand = demand_data[timesteps.+offset, 1]

    pv_data = nothing
    wind_data = nothing
    demand_data = nothing # Free the memory

    pars = copy(default_es_pars)

    pars[:recovery_time] = 12
    pars[:c_storage] = 300.0
    pars[:c_pv] = 800.0
    pars[:c_wind] = 1150.0
    pars[:c_in] = 0.165
    pars[:c_out] = 0.02
    pars[:asset_lifetime] = 5.0
    pars[:inv_budget] = 10000000.0
    pars[:sto_ef_ch] = 0.97
    pars[:sto_ef_dis] = 0.97
    pars[:feedincap] = 1e7
    pars[:max_pv] = 10^4
    pars[:max_wind] = 10^4
    return pv, wind, demand, heatdemand, pars
end

function load_basic_example(
    timesteps;
    basepath=joinpath(@__DIR__, ".."),
    offset=0,
)
    data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)
    heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)

    pv = data[timesteps.+offset, 3]
    wind = data[timesteps.+offset, 4]
    demand = data[timesteps.+offset, 2]
    heatdemand = heatdemand_data[timesteps.+offset, 1]

    data = nothing
    heatdemand_data = nothing # Free the memory

    pars = copy(default_es_pars)

    pars[:recovery_time] = 12
    pars[:c_storage] = 100.0
    pars[:c_pv] = 300.0
    pars[:c_wind] = 550.0
    pars[:penalty] = 1000000.0
    return pv, wind, demand, heatdemand, pars
end

function load_quartier(
    timesteps;
    basepath=joinpath(@__DIR__, ".."),
    offset=0,
)
    pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_pv_pik.csv"), DataFrame, header=false)
    wind_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_wind_pik.csv"), DataFrame, header=false)
    demand_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "demand_Industriepark.csv"), DataFrame, header=false)

    heatdemand = zeros(length(timesteps))
    pv = pv_data[timesteps.+offset, 1]
    wind = wind_data[timesteps.+offset, 1]
    demand = demand_data[timesteps.+offset, 1]

    pv_data = nothing
    wind_data = nothing
    demand_data = nothing # Free the memory

    pars = copy(default_es_pars)

    pars[:recovery_time] = 12
    pars[:c_storage] = 100.0
    pars[:c_pv] = 300.0
    pars[:c_wind] = 550.0
    pars[:penalty] = 1000000.0

    return pv, wind, demand, heatdemand, pars
end

function load_costs(
    run_id,
    n_samples, 
    scen_freq;
    basepath=joinpath(@__DIR__, ".."),
    savepath=joinpath(basepath, "results")
)
    return CSV.read(joinpath(savepath, run_id, "costs$(n_samples)_$(scen_freq).csv"), DataFrame)[!, 1]
end

function load_costs(
    run_id, 
    param::String;
    basepath=joinpath(@__DIR__, ".."),
    savepath=joinpath(basepath, "results"),
)
    return CSV.read(joinpath(savepath, run_id, "costs" * param * ".csv"), DataFrame)[!, 1]
end

function load_invs(
    run_id, 
    n_samples, 
    scen_freq;
    basepath=joinpath(@__DIR__, ".."),
    savepath=joinpath(basepath, "results"),
)
    inv_data = CSV.read(joinpath(savepath, run_id, "inv$(n_samples)_$(scen_freq).csv"), DataFrame)
    inv = Dict(pairs(eachcol(inv_data)))
    inv_data = nothing
    mean_inv = Dict((
        [v => mean(inv[v]) for v in keys(inv)]
    ))
    return inv, mean_inv
end

function load_ops(
    run_id,
    n_samples,
    scen_freq;
    basepath=joinpath(@__DIR__, ".."),
    savepath=joinpath(basepath, "results")
)
    op_data = CSV.read(joinpath(savepath, run_id, "operation$(n_samples)_$(scen_freq).csv"), DataFrame)
    op = Dict(pairs(eachcol(op_data)))
    op_data = nothing
    mean_op = Dict((
        [v => mean(op[v]) for v in keys(op)]
    ))
    return op, mean_op
end

function load_runtime(
    run_id, 
    n_samples, 
    scen_freq;
    basepath=joinpath(@__DIR__, ".."),
    savepath=joinpath(basepath, "results")
)
    return CSV.read(joinpath(savepath, run_id, "runtime$(n_samples)_$(scen_freq).csv"), DataFrame)[!, 1]
end