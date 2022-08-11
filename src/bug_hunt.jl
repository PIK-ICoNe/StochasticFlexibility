# Why are there times in which our algo finds 0 flexibility in either direction?

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using DataFrames
using CSV
using Clp
using Statistics;

using Random
Random.seed!(1);

#-

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#

offset = 24*7*14
timesteps = 1:(24*7*2)

data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)
heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"), DataFrame)

pv = data[timesteps .+ offset, 3]
wind = data[timesteps .+ offset, 4]
demand = data[timesteps .+ offset, 2]
heatdemand = heatdemand_data[timesteps .+ offset, 1]

data = nothing; # Free the memory
heatdemand_data = nothing;

plt = plot(timesteps, pv .* (mean(demand) / mean(pv)), label="PV (unitless)")
plot!(plt, timesteps, wind.* (mean(demand) / mean(wind)), label="Wind (unitless)")
plot!(plt, timesteps, heatdemand, label="Heat Demand")
plot!(plt, timesteps, demand, label="Electric Demand")
plt
#-
pars = copy(default_es_pars)

average_hourly_demand = mean(demand)

pars[:recovery_time] = 24
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 450.
pars[:c_sto_op] = 0.00001;
#-

#=
The model itself is constructed by the function define_energy_system
=#

es = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true)

#-

sp_no_flex = instantiate(es, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)
set_silent(sp_no_flex)

optimize!(sp_no_flex)

no_flex_decision = optimal_decision(sp_no_flex)

objective_value(sp_no_flex)

#-

analysis_window = 190+1:190+48

cost_pos, pot_pos, cost_neg, pot_neg = analyze_flexibility_potential(sp_no_flex, analysis_window)

plot_flexibility(analysis_window, cost_pos, pot_pos, cost_neg, pot_neg)

#- No flexibility between 215 and 220:

analysis_window = 214:220

cost_pos, pot_pos, cost_neg, pot_neg = analyze_flexibility_potential(sp_no_flex, analysis_window)

plot_flexibility(analysis_window, cost_pos, pot_pos, cost_neg, pot_neg)

#-

@show find_f_max(sp_no_flex, 216, 1, no_flex_decision; tol = 10., maxiter = 100)
@show find_f_max(sp_no_flex, 216, -1, no_flex_decision; tol = 10., maxiter = 100);

#-

# Here is the crux of it: The cost of the scenario with F_xi is Infinity even though this scenario shoudl be identical
# to the base case that was evaluated above.

scen = @scenario t_xi = 217 s_xi = -1 F_xi = 0. probability = 1.
cost_a = evaluate_decision(sp_no_flex,no_flex_decision,scen)

#-

costs_neg = []
costs_pos = []
for t in 1:300
    scen = @scenario t_xi = t s_xi = -1 F_xi = 0. probability = 1.
    push!(costs_neg, evaluate_decision(sp_no_flex,no_flex_decision,scen))
    scen = @scenario t_xi = t s_xi = +1 F_xi = 0. probability = 1.
    push!(costs_pos, evaluate_decision(sp_no_flex,no_flex_decision,scen))
end

#-

# no_flex_pseudo_sampler()[1] corresponds to
# @scenario t_xi = 1 s_xi = 1 F_xi = 0. probability = 1.
# so costs_pos[1]
# All these costs should be the same, yet:

@show count(isinf, costs_pos)
histogram(costs_pos)

#- at least they are the same, so the sign doesn not matter.

@show all(costs_pos .== costs_neg)

#-

# From reading the evaluate_decision code:

scen = @scenario t_xi = 218 s_xi = 1 F_xi = 0. probability = 1.
outcome = outcome_model(sp_no_flex, optimal_decision(sp_no_flex), scen; optimizer = subproblem_optimizer(sp_no_flex))
optimize!(outcome)
@show objective_value(outcome)
@show termination_status(outcome)

#- let's write a function that plots us the recovery of a scenario:

function plot_outcome(sp_base, t_xi, s_xi, F_xi; window_start=-2, window_end=2)
    scen = @scenario t_xi = t_xi s_xi = s_xi F_xi = F_xi probability = 1.
    sp = outcome_model(sp_base, optimal_decision(sp_base), scen; optimizer = subproblem_optimizer(sp_base))
    optimize!(sp)
    if termination_status(sp) != JuMP.MathOptInterface.OPTIMAL
        println("No optimum")
        return termination_status(sp)
    end

    recovery_window = t_xi:t_xi+length(sp[:gco2])-1
    plot_window = t_xi+1+window_start:t_xi+length(sp[:gco2])+window_end

    plt_gb = plot(title="grid buy", legend=:outertopright)
    plot!(plt_gb, plot_window, value.(sp[:gci][plot_window]) .- value.(sp[:gco][plot_window]), label = "Base Case")
    plot!(plt_gb, recovery_window, value.(sp[:gci2]) .- value.(sp[:gco2]), label = "Event")

    plt_soc = plot(title="soc", legend=:outertopright)
    plot!(plt_soc, plot_window, value.(sp[:sto_soc][plot_window]), label = "Base Case")
    plot!(plt_soc, recovery_window, value.(sp[:sto_soc2]), label = "Event")

    plt_h_soc = plot(title="heat soc", legend=:outertopright)
    plot!(plt_h_soc, plot_window, value.(sp[:heat_sto_soc][plot_window]), label = "Base Case")
    plot!(plt_h_soc, recovery_window, value.(sp[:heat_sto_soc2]), label = "Event")

    plt_e2h = plot(title="e2h", legend=:outertopright)
    plot!(plt_e2h, plot_window, value.(sp[:flow_energy2heat][plot_window]), label = "Base Case")
    plot!(plt_e2h, recovery_window, value.(sp[:flow_energy2heat2]), label = "Event")

    plot(plt_gb, plt_h_soc, plt_soc, plt_e2h, layout=(4,1))
end

#-

plot_outcome(sp_no_flex, 218, -1, 50000)

# Reading the heat balance equations I spotted that possibly heat demand in the scenario is not shifted appropriately...
# fixed this and rerun.... Seems to have fixed it mostly...