using DataFrames
using CSV

# pass heat = "district" for district heating, heat = "when2heat" for germany
    # "district" is default
function load_max_boegl(timesteps; offset = 0, heat = nothing, scale_factor = 1.) 
    pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "pv_Halle18.csv"), DataFrame, header = false)
    wind_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "wind_Karholz.csv"), DataFrame, header = false)
    demand_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "demand_Industriepark.csv"), DataFrame, header = false)
    if !isnothing(heat)
        if heat == "when2heat"
            heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)
            heatdemand = heatdemand_data[timesteps .+ offset, 1]./scale_factor
        else
            heatdemand_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_heatdemand_apartment_block(KFW40)_300km2_12300MWhperyear.csv"), DataFrame, header = false)
            heatdemand = heatdemand_data[timesteps .+ offset, 1]./scale_factor
        end
    else
        heatdemand = zeros(length(timesteps))
    end
    pv = pv_data[timesteps .+ offset, 1]
    wind = wind_data[timesteps .+ offset, 1]
    demand = demand_data[timesteps .+ offset, 1]

    pv_data = nothing;
    wind_data = nothing;
    demand_data = nothing; # Free the memory

    pars = copy(default_es_pars)

    pars[:recovery_time] = 12
    pars[:c_storage] = 300.
    pars[:c_pv] = 800.
    pars[:c_wind] = 1150.
    pars[:c_in] = 0.165
    pars[:c_out] = 0.02
    pars[:asset_lifetime] = 5.
    pars[:inv_budget] = 10000000.
    pars[:sto_ef_ch] = 0.97
    pars[:sto_ef_dis] = 0.97
    pars[:feedincap] = 1e7
    pars[:max_pv] = 10^4
    pars[:max_wind] = 10^4;
    return pv, wind, demand, heatdemand, pars
end

function load_max_boegl_district_heating(timesteps; offset = 0, heat = false, heat_timeseries = "district")
    pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "pv_Halle18.csv"), DataFrame, header = false)
    wind_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "wind_Karholz.csv"), DataFrame, header = false)
    demand_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "demand_Industriepark.csv"), DataFrame, header = false)
    if heat
        heatdemand_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_heatdemand_apartment_block(KFW40)_300km2_12300MWhperyear.csv"), DataFrame, header = false)
        heatdemand = heatdemand_data[timesteps .+ offset, 1]
    else
        heatdemand = zeros(length(timesteps))
    end
    pv = pv_data[timesteps .+ offset, 1]
    wind = wind_data[timesteps .+ offset, 1]
    demand = demand_data[timesteps .+ offset, 1]

    pv_data = nothing;
    wind_data = nothing;
    demand_data = nothing; # Free the memory

    pars = copy(default_es_pars)

    pars[:recovery_time] = 12
    pars[:c_storage] = 300.
    pars[:c_pv] = 800.
    pars[:c_wind] = 1150.
    pars[:c_in] = 0.165
    pars[:c_out] = 0.02
    pars[:asset_lifetime] = 5.
    pars[:inv_budget] = 10000000.
    pars[:sto_ef_ch] = 0.97
    pars[:sto_ef_dis] = 0.97
    pars[:feedincap] = 1e7
    pars[:max_pv] = 10^4
    pars[:max_wind] = 10^4;
    return pv, wind, demand, heatdemand, pars
end

function load_basic_example(timesteps; offset = 0)
    data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)
    heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)

    pv = data[timesteps .+ offset, 3]
    wind = data[timesteps .+ offset, 4]
    demand = data[timesteps .+ offset, 2]
    heatdemand = heatdemand_data[timesteps .+ offset, 1]

    data = nothing; 
    heatdemand_data = nothing; # Free the memory

    pars = copy(default_es_pars)

    pars[:recovery_time] = 12
    pars[:c_storage] = 100.
    pars[:c_pv] = 300.
    pars[:c_wind] = 550.
    pars[:penalty] = 1000000.;
    return pv, wind, demand, heatdemand, pars
end

function load_quartier(timesteps; offset = 0)
    pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_pv_pik.csv"), DataFrame, header = false)
    wind_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_wind_pik.csv"), DataFrame, header = false)
    demand_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "demand_Industriepark.csv"), DataFrame, header = false)

    heatdemand = zeros(length(timesteps))
    pv = pv_data[timesteps .+ offset, 1]
    wind = wind_data[timesteps .+ offset, 1]
    demand = demand_data[timesteps .+ offset, 1]

    pv_data = nothing;
    wind_data = nothing;
    demand_data = nothing; # Free the memory

    pars = copy(default_es_pars)

    pars[:recovery_time] = 12
    pars[:c_storage] = 100.
    pars[:c_pv] = 300.
    pars[:c_wind] = 550.
    pars[:penalty] = 1000000.;

    return pv, wind, demand, heatdemand, pars
end

function load_costs(run_id, n_samples, scen_freq; savepath = joinpath(basepath, "results"))
    return CSV.read(joinpath(savepath, run_id, "costs$(n_samples)_$(scen_freq).csv"), DataFrame)[!,1]
end

function load_costs(run_id, param::String; savepath = joinpath(basepath, "results"))
    return CSV.read(joinpath(savepath, run_id, "costs"*param*".csv"), DataFrame)[!,1]
end

function load_invs(run_id, n_samples, scen_freq; savepath = joinpath(basepath, "results"))
    inv_data = CSV.read(joinpath(savepath, run_id, "inv$(n_samples)_$(scen_freq).csv"), DataFrame)
    inv = Dict(pairs(eachcol(inv_data)))
    inv_data = nothing;
    mean_inv = Dict((
        [v => mean(inv[v]) for v in keys(inv)]
    ))
    return inv, mean_inv
end

function load_ops(run_id, n_samples, scen_freq; savepath = joinpath(basepath, "results"))
    op_data = CSV.read(joinpath(savepath, run_id, "operation$(n_samples)_$(scen_freq).csv"), DataFrame)
    op = Dict(pairs(eachcol(op_data)))
    op_data = nothing;
    mean_op = Dict((
        [v => mean(op[v]) for v in keys(op)]
    ))
    return op, mean_op
end

function load_runtime(run_id, n_samples, scen_freq; savepath = joinpath(basepath, "results"))
    return CSV.read(joinpath(savepath, run_id, "runtime$(n_samples)_$(scen_freq).csv"), DataFrame)[!,1]
end