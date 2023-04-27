## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
Pkg.instantiate()

using DataFrames
using CSV
using Clp
using JSON
using BSON
using Statistics;


using Random
Random.seed!(1);

#-

#=
# Flexibility analysis and optimization

We take an energy system and analyze and optimize its flexibility potential.

The model is defined in the file sp_model.jl, analysis and utility functions are in the other two files.

The parametrization is not completely done yet, and running the analysis for a whole year is very time intensive.
The main optimization runs fine for a whole year, but analyzing the flexibility potential for many hours
is computationally intense.

We need to either start caching more computations, or find a way to speed up the evaluation further.

This file gives a conceptual demonstration of what we do on a two week horizon.
=#

#- 

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#
t_max_offset = 24 #not using the last day of the year, it is not allowed to have a request in the last day  
offset = 0
timesteps = 1:(24*365)

#loading a timeseries for photovoltaic (pv) and wind potial and the demnad for electricity and heat"
pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_pv_pik.csv"), DataFrame, header = false)
wind_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_wind_pik.csv"), DataFrame, header = false)
demand_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_apartment_block(KFW40)_300km2_6600MWhperyear.csv"), DataFrame, header = false)

heatdemand_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_heatdemand_apartment_block(KFW40)_300km2_12300MWhperyear.csv"), DataFrame, header = false)

#defining the offset of each timesiries
pv = pv_data[timesteps .+ offset, 1]
wind = wind_data[timesteps .+ offset, 1]
demand = demand_data[timesteps .+ offset, 1]
heatdemand = heatdemand_data[timesteps .+ offset, 1]
#heatdemand = zeros(length(timesteps)) #if heatdemand is not needed 

t_max = minimum((length(pv), length(wind), length(demand), length(heatdemand))) - t_max_offset #using the longest possible timeintervall 

data = nothing; # Free the memory
heatdemand_data = nothing;
pv_data = nothing;
wind_data = nothing;
demand_data = nothing; # Free the memory

plt = plot(timesteps, pv .* (mean(demand) / mean(pv)), label="PV (unitless)")
plot!(plt, timesteps, wind.* (mean(demand) / mean(wind)), label="Wind (unitless)")
plot!(plt, timesteps, heatdemand, label="Heat Demand")
plot!(plt, timesteps, demand, label="Electric Demand")
plt
#-
## Paul: maybe using refreshed one 
#=
Next we continue the set up. Our model comes with default parameters,
which we slightly adjust here. We use some arbitrary values to define a dummy heat demand.
=#

pars = copy(default_es_pars)

average_hourly_demand = mean(demand)

recovery_time = 12

pars[:recovery_time] = recovery_time
pars[:c_storage] = 600.
pars[:c_pv] = 800.
pars[:c_wind] = 2500.
pars[:c_i] = 0.4
pars[:c_o] = 0.04
pars[:c_heat_storage] = 400.
pars[:asset_lifetime] = 20.
pars[:c_heatpump] = 457.; #kW heat -> mayve wrong and it should be 1600 kW electricity


savepath = joinpath(basepath, "results")
if !isdir(savepath)
    mkdir(savepath)
end
#=
The model itself is constructed by the function define_energy_system
=#

es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = false, override_no_event_per_scen = true)
es_no_reg = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = false)


#-

#=
This contains investment variables which we collectively call $I$, and the operational schedule $O^t$.
The total cost given a certain investment $I$ and schedule $O^t$ is denoted $C(I, O^t)$.

We now can optimize the system, initialy while ignoring flexibility:
=#

sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
#no flex sampler 

# Prevent investment in the heat components if there is no heat demand
if maximum(heatdemand) == 0.
    fix(decision_by_name(sp_bkg, 1, "u_heatpump"), 0.)
    fix(decision_by_name(sp_bkg, 1, "u_heat_storage"), 0.)
end

set_silent(sp_bkg) #no direct output from the solver 

optimize!(sp_bkg)

bkg_decision = optimal_decision(sp_bkg)

objective_value(sp_bkg)

all_data = get_all_data(sp_bkg, rec = false, scen = false)
all_data = merge(all_data, Dict((:cost => objective_value(sp_bkg))))
open(joinpath(savepath, "bkg.json"), "w") do f
    JSON.print(f,all_data)
end

#-

plot_results(all_data, pv, wind, demand)
#-
plot_heat_layer(all_data, heatdemand)
#-
all_data = nothing #clean the memory 
#=
We can analyze individual scenarios by using the evaluate decision function.
This allows us to test the system on scenarios that were not present in the original problem definition

To get a model that corresponds to a single flexibility event (e.g. to investigate how the system 
is reacting to any one event) we can construct outcome models as below. These are normal JuMP models.

=#
t_max = length(pv) - 24
F_max = 10000.
F_min = 3000.
delta_t = 3*24 # Flex event every week
pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1);
n_samples = [collect(15:5:35); collect(40:10:90); collect(100:25:175)]
for ns in n_samples
    n = round(Int, ns * pars[:event_per_scen])
    scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
    es_st = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = true) #defining the stochastic energy system with flexibility
    sp_st = instantiate(es_st, scens, optimizer = Clp.Optimizer)
    set_silent(sp_st)
    optimize!(sp_st)
    println("termination_status=", termination_status(sp_st)) #sucessfull or infeasible

    opt_params = Dict((:F_min => F_min, :F_max => F_max, :t_max_offset => t_max_offset, :n_samples => n_samples, :scen_freq => scen_freq))
    all_data = get_all_data(sp_st)
    all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp_st))))
    bson("run_$(n_samples)_$(scen_freq).bson", all_data)
    if ns == n_samples[end]
        ns = n_samples[end]
        es_st = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = true) #defining the stochastic energy system with flexibility
        sp_st = instantiate(es_st, scens, optimizer = Clp.Optimizer)
        fix_investment!(sp_st, get_investments(sp_bkg))
        set_silent(sp_st)
        optimize!(sp_st)
        opt_params = Dict((:F_min => F_min, :F_max => F_max, :t_max_offset => t_max_offset, :n_samples => n_samples, :scen_freq => scen_freq))
        all_data = get_all_data(sp_st)
        all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp_st))))
        bson("run_$(n_samples)_$(scen_freq)_fixed_inv.bson", all_data)
    end
end



costs = zeros(length(n_samples))
u_pv = zeros(length(n_samples))
u_wind = zeros(length(n_samples))
u_storage = zeros(length(n_samples))
u_heat_storage = zeros(length(n_samples))
u_heatpump = zeros(length(n_samples))

for (i,ns) in enumerate(n_samples) 
    opt = BSON.load("run_$(n_samples)_$(scen_freq).bson")
    costs[i]= opt[:cost]
    u_pv[i]= opt[:inv][:u_pv]
    u_wind = 
end


#finishing the loop above 
#creating the plots (see notes)
#flex potential function -> read two bson files naive_flex_potential -> look in the parameters, it returns f_pos and f_neg



