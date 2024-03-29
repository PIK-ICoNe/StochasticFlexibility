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
include(joinpath(basepath, "src", "evaluation_utils.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

#-

#=
We load timeseries for photovoltaic (pv) and wind potential as well as demand.
=#

timesteps = 1:(24*365)
pv, wind, demand, heatdemand = load_basic_example(timesteps);

plt = plot(timesteps, pv .* (mean(demand) / mean(pv)), label="PV (unitless)")
plot!(plt, timesteps, wind.* (mean(demand) / mean(wind)), label="Wind (unitless)")
plot!(plt, timesteps, heatdemand, label="Heat Demand")
plot!(plt, timesteps, demand, label="Electric Demand")
plt
#-

#=
Next we continue the set up. Our model comes with default parameters,
which we slightly adjust here. We use some arbitrary values to define a dummy heat demand.
=#

# pars = copy(default_es_pars)

# average_hourly_demand = mean(demand)

# pars[:recovery_time] = 24
# pars[:c_storage] = 100.
# pars[:c_pv] = 300.
# pars[:c_wind] = 550.
# pars[:c_sto_op] = 0.00001;

# pars[:]
# #-


pars = copy(default_es_pars)

average_hourly_demand = mean(demand)

recovery_time = 12

pars[:recovery_time] = recovery_time
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 550.
pars[:penalty] = 1000000.;

t_max = length(pv) - 24
F_max = 10000.
F_min = 3000.
delta_t = 3*24 # Flex event every week
pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1);
n = round(Int, 10 * pars[:event_per_scen])

scens = poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)

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

set_silent(sp_bkg)

optimize!(sp_bkg)

bkg_decision = optimal_decision(sp_bkg)

objective_value(sp_bkg)

#-

#=
```math
\overline I_{nf} , \overline O^t_{nf} = \argmin_{I, O^t} C(I, O^t)
```

We can now analyse how much flexibility is available at each point in time with this investment and this schedule.

To do so we distinguish between parts of the operational schedule that can react to a flexibility demand in the moment that it occurs, and those that can only adapt later to help recover the system back towards its previously planned schedule.
Essentially into slow components $O^{s,t}$ and fast components $O^{f,t}$

```math
O^t = (O^{f,t}, O^{s,t})
```

We then consider the cost of the system when confronted with an additional negative or positive demand $F$ at some timepoint $t^F$.
At the timepoint $t^F$ only the fast reacting components can adjust their schedule, the slow components still have to stick to the previous schedule.
After the event we give the system a time $t_r$ to recover back to its original schedule. Here the slow components can contribute. Thus we have:

```math
\begin{aligned}
c(F,t^F) &= \min_{O^t} C'(\overline I_{nf}, O^t, F, t_f)  \\
O^{s,t} &= \overline O^t_{nf} \;\; \forall \;\;\; t \leq t^F , \;\;\; t > t^F + t_r \\
O^{f,t} &= \overline O^t_{nf} \;\; \forall  \;\;\; t < t^F , \;\;\; t > t^F + t_r
\end{aligned}
```

This is the cost of flexibility in the sense of Harder et.al. This optimization is not always feasible.
The maximum and minimum $F$ for which it is feasible is called the positive and negative flexibility potential $pot_\pm$.
The marginal cost of flexibility is $cost_{\pm}(t) = c(pot_\pm(t), t) / pot_\pm(t)$.
Using our model above we can analyze this in the following way:
=#

# analysis_window = 190+1:190+48

# cost_pos, pot_pos, cost_neg, pot_neg = analyze_flexibility_potential(sp_no_flex, analysis_window)

# plot_flexibility(analysis_window, cost_pos, pot_pos, cost_neg, pot_neg)

#-

plot_results(sp_bkg, pv, wind, demand)
#-
plot_heat_layer(sp_bkg, heatdemand)
#-

#=
We can analyze individual scenarios by using the evaluate decision function.
This allows us to test the system on scenarios that were not present in the original problem definition

To get a model that corresponds to a single flexibility event (e.g. to investigate how the system 
is reacting to any one event) we can construct outcome models as below. These are normal JuMP models.

=#

prob_scen = @scenario t_xi = 195 s_xi = 1. F_xi = 20. probability = 1.
evaluate_decision(sp_bkg, bkg_decision, prob_scen)

#-

outcome = outcome_model(sp_bkg, bkg_decision, prob_scen; optimizer = subproblem_optimizer(sp_bkg))
set_silent(outcome)
optimize!(outcome)
termination_status(outcome)

#-

#=
Under the hood this uses the evaluate_decision function of stochastic programs that evaluates the cost of a specified scenario given a decision.

We can go further if we introduce a certain distribution of flexibility demands.
Then we can choose to operate the system in such a way that the expected cost, including the expected cost of flexibility, is minimized.
This is exactly a two stage stochastic programming problem with $c(F,t)$ as the second stage.

We can start by taking a simple uniform distribution between some maximum flexibility demand that can occur at any time that leaves enough space for the recovery window:
=#

# n = 100
# F_max = average_hourly_demand * 0.1 # Have unaticipated demand equal to 10% of our typical demand
# t_max = length(pv) - es.parameters[2].defaults[:recovery_time]

# scens = simple_flex_sampler(n, F_max, t_max);

#-

#=
We can now evaluate the expected cost of running the system determined above with this flexibility distribution.
=#

sp_no_reg = instantiate(es_no_reg, scens, optimizer = Clp.Optimizer)
set_silent(sp_no_reg)

evaluate_decision(sp_no_reg, bkg_decision)

#-
#=
This is infinite as the system as built above can not actually provide the desired flexibility at all times.
One way to deal with this problem is to regularize the problem, by allowing a heavily penalized deviation from satisfying the extra demand.
=#

es_reg = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = true)

sp_reg = instantiate(es_reg, scens, optimizer = Clp.Optimizer)
set_silent(sp_reg)

#-
#=
Now we can evaluate the expected cost of flexibility in the system: 
=#

evaluate_decision(sp_reg, bkg_decision)

#-

#=
TODO / BUG: This is infinite for some time windows/seeds! This should not happen. (SOLVED?!)

This cost is of course completely dominated by the regularizer. However, using stochastic programming we can optimize the operational schedule to directly optimize this expected cost.
To do so we fix the investment to the system we have and then optimize the remaining variables:
=#


investments_nf = get_investments(sp_bkg)
fix_investment!(sp_no_reg, investments_nf)
fix_investment!(sp_reg, investments_nf)

optimize!(sp_no_reg)
optimize!(sp_reg)

termination_status(sp_no_reg) # Infeasible with F_max 10000 and 1000
termination_status(sp_reg)

# no_reg_no_invest_decision = optimal_decision(sp_no_reg)
reg_no_invest_decision = optimal_decision(sp_reg);

#-

# prob = 1/length(scens)
# is_optimal = zeros(Bool, length(scens))

# Threads.@threads for i in eachindex(scens)
#     println(i)
#     s = scens[i]
#     s_new = @scenario t_xi = s.data[:t_xi] s_xi = s.data[:s_xi] F_xi = 10000. probability = prob
#     scenario_model = outcome_model(sp_no_reg, optimal_decision(sp_reg), s_new; optimizer = subproblem_optimizer(sp_reg))
#     set_silent(scenario_model)
#     optimize!(scenario_model) # Primal infeasible
#     # println(termination_status(scenario_model))
#     # println(s_new)
#     is_optimal[i] = termination_status(scenario_model) == MOI.TerminationStatusCode(1)
# end
#-
#=
Then evaluating the decision we find much more reasonable values:
=#

evaluate_decision(sp_reg, no_reg_no_invest_decision)

#-

evaluate_decision(sp_reg, reg_no_invest_decision)

#-

no_invest_relative_cost_of_flex = evaluate_decision(sp_reg, reg_no_invest_decision) / objective_value(sp_bkg) - 1.

#-

#=
It is interesting to note that this way servicing the flexibility is actually not very expensive
compared to the baseline cost of the system if there were no flexibility demands, with the relative cost of flexibility coming in at 0.2% here.

We can also again evaluate the amount of flexibility available at each point in time, and the associated cost:
=#

# To find flexibility potentials we right now have to use the unregularized model:
# cost_pos_flex, pot_pos_flex, cost_neg_flex, pot_neg_flex = analyze_flexibility_potential(sp_flex, analysis_window)

# plot_flexibility(analysis_window, cost_pos_flex, pot_pos_flex, cost_neg_flex, pot_neg_flex)

#-
#=
This shows that the model is not actually able to guarantee that there is flexibility at all times. However, it dramatically increases the amount of available flexibility:
=#

# plt_av = plot();
# flexibility_availability!(plt_av, pot_pos, label = "positive flexbility unaware", c = :red);
# flexibility_availability!(plt_av, pot_pos_flex, label = "positive flexibility aware", c = :green);
# flexibility_availability!(plt_av, pot_neg, label = "negative flexibility unaware", c = :red);
# flexibility_availability!(plt_av, pot_neg_flex, label = "negative flexibility aware", c = :green);
# plt_av

#-

#=
We can of course also optimize the overall system investment to take flexibility into account.
=#

unfix_investment!(sp_no_reg, investments_nf)
unfix_investment!(sp_reg, investments_nf)

optimize!(sp_no_reg)
optimize!(sp_reg)

termination_status(sp_no_reg)
termination_status(sp_reg)

no_reg_invest_decision = optimal_decision(sp_no_reg)
reg_invest_decision = optimal_decision(sp_reg);

#=
Then the relative cost of the system exposed to flexibility is further reduced:
=#

evaluate_decision(sp_reg, no_reg_invest_decision)
#-

evaluate_decision(sp_reg, reg_invest_decision)

#-

invest_relative_cost_of_flex = evaluate_decision(sp_reg, reg_invest_decision) / objective_value(sp_bkg) - 1.

#-

#=
Considering the flexiblity demands at investment time, rather than only during operations lowers the cost of flexibility by:
=#

invest_relative_cost_of_flex / no_invest_relative_cost_of_flex - 1.

#=
27% for this toy model.

# Further thoughts

There are numerous directions to go from here.

This document focused on (a regularized version of) the expected price of flexibility.
Another possibility would be to make a minimum amount of available flexibility a constraint.
There are algorithms taht concern stochastic constraints that might be useful to explore then.

A general question with all of this is: We are sampling the space of possible events.
The operational schedule we obtain will probably contain moments at which there is no flexibility.
The sample might miss the few hours at which flexibility is hardest to come by.

It is unclear (and requires further analysis) how much of a problem that is. E.g. maybe the flex aware schedule
based on a sample will not have flexibility at some hours, but a very minor adjustment would.
In partiuclar I consider it plausible that the investment decisions are not strongly affected by these gaps,
as long as sufficiently many representative scenarios are sampled.

In evaluating the quality of the sampling approach it might also be appropriate to evaluate the assumption of
a schedule based on perfect foresight for everything _but_ the felxiblity. A proper validation set up would also
consider the weather and demand uncertainties.

A general choice we have is whether to think of the flexibility as for the system itself, or for selling as an auxilliary service.
In the former case we can think of the regularizer as buying flexibility on the open market once providing it ourselfs becomes to expensive.
For the latter we would need to think about potential products that can reasonably be offered/modelled in the stochastic program.
E.g. we can not simply say we are offering the full flexibility potential as that is a quantity that is hard to evaluate.

Given multiple energy systems that trade with each other we could try to analyze a market of such systems as well

=#

