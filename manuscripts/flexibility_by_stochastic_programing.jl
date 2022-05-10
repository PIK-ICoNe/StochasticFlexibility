## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath((@__DIR__, "..")))

using Pkg
Pkg.activate(basepath)
Pkg.instantiate()

using DataFrames
using CSV
using Clp

#-

#=
# Flexibility analysis and optimization

We take an energy system and analyze and optimize its flexibility potential.

Te model is defined in the file sp_model.jl, analysis and utility functions are in the other two files.
=#

#- 

include(joinpath((basepath, "src", "sp_model.jl")))
include(joinpath((basepath, "src", "plot_utils.jl")))
include(joinpath((basepath, "src", "evaluation_utils.jl")))

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#

offset = 6531
timesteps = 1:168

pv = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps .+ offset, 3]
wind = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps .+ offset, 4]
demand = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps .+ offset, 2]

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#

