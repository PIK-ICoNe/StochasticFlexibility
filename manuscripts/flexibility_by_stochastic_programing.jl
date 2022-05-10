## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath((@__DIR__, "..")))

using Pkg
Pkg.activate(basepath)
Pkg.instantiate()

using DataFrames
using CSV
using Clp
using Statistics

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

##

#=
Next we continue the set up. Our model comes with default parameters,
which we slightly adjust here. We use some arbitrary values to define a dummy heat demand.
=#

pars = copy(default_es_pars)

average_hourly_demand = mean(demand)

pars[:recovery_time] = 24
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 800.
pars[:c_sto_op] = 0.00001

heatdemand = copy(demand)./300.
heatdemand0 = zeros(length(demand))

##

#=
The model itself is constructed by the function define_energy_system
=#

es = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true)

##

#=
This contains investment variables which we collectively call $I$, and the operational schedule $O^t$.
The total cost given a certain investment $I$ and schedule $O^t$ is denoted $C(I, O^t)$.

We now can optimize the system, initialy while ignoring flexibility:
=#

sp_no_flex = instantiate(es, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)

optimize!(sp_no_flex)

no_flex_decision = optimal_decision(sp_no_flex)

objective_value(sp_no_flex)

##

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
c(F,t^F) = \min_{O^t} C'(\overline I_{nf}, O^t, F, t_f)
O^{s,t} = \overline O^t_{nf} \forall t \leq t^F t > t^F + t_r
O^{f,t} = \overline O^t_{nf} \forall t < t^F t > t^F + t_r
```

This is the cost of flexibility in the sense of Harder et.al. This optimization is not always feasible.
The maximum and minimum $F$ for which it is feasible is called the positive and negative flexibility potential $pot_\pm$.
The marginal cost of flexibility is $cost_{\pm}(t) = c(pot_\pm(t), t) / pot_\pm(t)$.
Using our model above we can analyze this in the following way:
=#

cost_pos, pot_pos, cost_neg, pot_neg = analyze_flexibility_potential(sp_no_flex, 20:50)

plot_flexibility(20:50, cost_pos, pot_pos, cost_neg, pot_neg)

##

#=
Under the hood this uses the evaluate_decision function of stochastic programs that evaluates the cost of a specified scenario given a decision.

We can go further if we introduce a certain distribution of flexibility demands.
Then we can choose to operate the system in such a way that the expected cost, including the expected cost of flexibility, is minimized.
This is exactly a two stage stochastic programming problem with $c(F,t)$ as the second stage.

We can start by taking a simple uniform distribution between some maximum flexibility demand that can occur at any time that leaves enough space for the recovery window:
=#

n = 100
F_max = average_hourly_demand * 0.1 # Have unaticipated demand equal to 10% of our typical demand
t_max = length(pv) - es.parameters[2].defaults[:recovery_time]

scens = simple_flex_sampler(n, F_max, t_max)

##

#=
We can now evaluate the expected cost of running the system determined above with this flexibility distribution.
=#

sp_flex = instantiate(es, scens, optimizer = Clp.Optimizer)

evaluate_decision(sp_flex, no_flex_decision)

##
#=
This is infinite as the system as built above can not actually provide the desired flexibility at all times.
One way to deal with this problem is to regularize the problem, by allowing a heavily penalized deviation from satisfying the extra demand.
=#

es_reg = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = false)

sp_reg_flex = instantiate(es_reg, scens, optimizer = Clp.Optimizer)


##
#=
Now we can evaluate the expected cost of flexibility in the system: 
=#

evaluate_decision(sp_reg_flex, no_flex_decision)

##

#=
This cost is of course completely dominated by the regularizer. However, using stochastic programming we can optimize the operational schedule to directly optimize this expected cost.
To do so we fix the investment to the system we have and then optimize the remaining variables:
=#

investments_nf = get_investments(sp_no_flex)
fix_investment!(sp_flex, investments_nf)
fix_investment!(sp_reg_flex, investments_nf)

optimize!(sp_flex)
optimize!(sp_reg_flex)

flex_no_invest_decision = optimal_decision(sp_flex)
reg_flex_no_invest_decision = optimal_decision(sp_reg_flex)

##
#=
Then evaluating the decision we find much more reasonable values:
=#

evaluate_decision(sp_reg_flex, flex_no_invest_decision)

##
evaluate_decision(sp_reg_flex, reg_flex_no_invest_decision)

##

relative_flex_cost = evaluate_decision(sp_reg_flex, flex_no_invest_decision) / objective_value(sp_no_flex) - 1.

##

#=
It is interesting to note that this way servicing the flexibility is actually not very expensive
compared to the baseline cost of the system if there were no flexibility demands, with the relative cost of flexibility coming in at 0.2% here.

We can also again evaluate the amount of flexibility available at each point in time, and the associated cost:
=#

# To find flexibility potentials we right now have to use the unregularized model:
cost_pos_flex, pot_pos_flex, cost_neg_flex, pot_neg_flex = analyze_flexibility_potential(sp_flex, 20:50)

plot_flexibility(20:50, cost_pos_flex, pot_pos_flex, cost_neg_flex, pot_neg_flex)

##
#=
This shows that the model is not actually able to guarantee that there is flexibility at all times. However, it dramatically increases the amount of available flexibility:
=#

plt_av = plot();
flexibility_availability!(plt_av, pot_pos, label = "positive flexbility unaware");
flexibility_availability!(plt_av, pot_pos_flex, label = "positive flexibility aware");
flexibility_availability!(plt_av, pot_neg, label = "negative flexibility unaware");
flexibility_availability!(plt_av, pot_neg_flex, label = "negative flexibility aware");
display(plt_av)

##

#=
We can of course also optimize the overall system investment to take flexibility into account.
=#

unfix_investment!(sp_flex, investments_nf)
unfix_investment!(sp_reg_flex, investments_nf)

optimize!(sp_flex)
optimize!(sp_reg_flex)

flex_invest_decision = optimal_decision(sp_flex)
reg_flex_invest_decision = optimal_decision(sp_reg_flex)

#=
Then the relative cost of the system exposed to flexibility is further reduced:
=#

evaluate_decision(sp_reg_flex, flex_invest_decision)

##
evaluate_decision(sp_reg_flex, reg_flex_invest_decision)

##

relative_flex_cost_inv = evaluate_decision(sp_reg_flex, flex_invest_decision) / objective_value(sp_no_flex) - 1.
