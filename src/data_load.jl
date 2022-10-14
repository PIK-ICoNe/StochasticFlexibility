using DataFrames
using CSV

function load_max_boegl(timesteps; offset = 0)
    pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "pv_Halle18.csv"), DataFrame, header = false)
    wind_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "wind_Karholz.csv"), DataFrame, header = false)
    demand_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "demand_Industriepark.csv"), DataFrame, header = false)
    heatdemand = zeros(length(timesteps))
    pv = pv_data[timesteps .+ offset, 1]
    wind = wind_data[timesteps .+ offset, 1]
    demand = demand_data[timesteps .+ offset, 1]

    pv_data = nothing;
    wind_data = nothing;
    demand_data = nothing; # Free the memory
    return pv, wind, demand, heatdemand
end

function load_basic_example(timesteps; offset = 0)
    data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)
    heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)

    pv = data[timesteps, 3]
    wind = data[timesteps, 4]
    demand = data[timesteps, 2]
    heatdemand = heatdemand_data[timesteps, 1]

    data = nothing; 
    heatdemand_data = nothing; # Free the memory
    return pv, wind, demand, heatdemand
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
    return pv, wind, demand, heatdemand
end
