# Everything runs in the Project environment on the basepath

basepath = realpath(joinpath((@__DIR__, "..")))

using Pkg
Pkg.activate(basepath)
Pkg.instantiate()

using DataFrames
using CSV
using Clp

include(joinpath((basepath, "src", "sp_model.jl")))
include(joinpath((basepath, "src", "plot_utils.jl")))
include(joinpath((basepath, "src", "evaluation_utils.jl")))

offset = 6531
timesteps = 1:168

pv = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps .+ offset, 3]
wind = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps .+ offset, 4]
demand = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps .+ offset, 2]

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

