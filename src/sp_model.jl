using StochasticPrograms
using Random
using Distributions

##
"""
Default values of system parameters.
- c_i - price of buying energy from the grid, Euro/kW
- c_o - price of selling energy to the grid, Euro/kW
- asset_lifetime - expected lifetime of a component in years. Needed to bring investment and operational costs to the same timescale
- c_pv, c_wind - price of installing pv and wind components per kW peak, Euro/kWp
- c_storage - price of storage, Euro/kW
- c_heat_storage - TODO: heat layer units
- c_heatpump
- inv_budget - maximum investment budget, Euro
- recovery_time - time after the event, in which the system is allowed to deviate from background optimal schedule, h
- COP - heatpump efficiency coefficient
- heat_losses - losses in the heat storage
- penalty - price of non-delivery of flexibility, Euro
"""
default_es_pars = Dict((
    :c_i => .3,
    :c_o => .05,
    :asset_lifetime => 10.,
    :c_pv => 700.,
    :c_wind => 2000.,
    :c_storage => 600.,
    :c_heat_storage => 30.,
    :c_heatpump => 20.,
    :inv_budget => 500000000.,
    :recovery_time => 12,
    :COP => 3.5,
    :feedincap => 10^7,
    :heat_losses => 0.2,
    :heat_eff => 0.9,
    :sto_ef_ch => 0.95,
    :sto_ef_dis => 0.95,
    :penalty => 10000.,
    :event_per_scen => 1,
    :max_sto_flow => 0.2,
    :max_pv => 10^3,
    :max_wind => 10^3
))

"""
Draw a sample of n scenarios
"""
function simple_flex_sampler(n, F_max, t_max; F_min = 0.)
    [@scenario t_xi = rand(1:t_max) s_xi = rand([-1, 1]) F_xi = F_min + rand() * (F_max - F_min) probability = 1/n 
        for i in 1:n]
end

"""
Exponentially distributed inter event times, with the inter event time counting
from after the recovery time. Average inter event time should be delta_t + recovery_time + 1
"""
function poisson_events_with_offset(n, delta_t, recovery_time, F_max, t_max; F_min = 0.)
    @assert t_max > delta_t
    times = []
    t = 1
    interval_distribution = Exponential(delta_t)
    while length(times) < n
        t += recovery_time + 1
        t += round(Int, rand(interval_distribution))
        push!(times, t)
    end
    times = times .% t_max .+ 1
    [@scenario t_xi = t s_xi = rand([-1, 1]) F_xi = F_min + rand() * (F_max - F_min) probability = 1/n 
    for t in times]
end


"""
Get an array with a single pseudo scenario with no flexiblity. 
This is useful for optimizing the system as if no flexibilty was introduced.
"""
function no_flex_pseudo_sampler()
    [@scenario t_xi = 1 s_xi = 1 F_xi = 0. probability = 1.
    ,]
end

"""
Draw a sample of n scenarios which have different probaility at day and night 
"""

function timedependent_flex_sampler(n, F_max, t_day, t_night, prob_day, prob_night)
    sday = [@scenario t_xi = rand(t_day) s_xi = rand([-1, 1]) F_xi = rand() * F_max probability = prob_day/n 
        for i in 1:(n/2)]
    snight = [@scenario t_xi = rand(t_night) s_xi = rand([-1, 1]) F_xi = rand() * F_max probability = prob_night/n 
        for i in 1:(n/2)]
    return vcat(sday,snight)
end

"""
Define energy system.
Parameters:
- pv, wind - weather timeseries
- demand, heatdemand - demand timeseries
- p - dictionary with system parameters, such as component costs, losses and recovery time window
- regularized - bool, if true, finite penalty is used
"""
function define_energy_system(pv, wind, demand, heatdemand; p = default_es_pars, regularized = true, debug_cap = 10^9, reg_lossy_flows = 0.0000001, override_no_event_per_scen = false, guaranteed_flex = false, F_pos = 1000., F_neg = -1000., cap_constraint = "naive")
    number_of_hours = minimum([length(pv), length(demand), length(wind)])
    if override_no_event_per_scen
        event_per_scen = 0
    else
        event_per_scen = p[:event_per_scen]
    end
    if cap_constraint == "naive"
        energy_system = @stochastic_model begin 
            @stage 1 begin
                @parameters begin
                    # Expected lifetime of components, years
                    asset_lifetime = p[:asset_lifetime]
                    # Costs in Euro/1kWp
                    c_pv = p[:c_pv]
                    c_wind = p[:c_wind]
                    c_storage = p[:c_storage]
                    c_i = p[:c_i]
                    c_o = p[:c_o]
                    c_heat_storage = p[:c_heat_storage]
                    c_heatpump = p[:c_heatpump]
                    COP = p[:COP]
                    heat_losses = p[:heat_losses]
                    heat_eff = p[:heat_eff]
                    sto_ef_ch = p[:sto_ef_ch] # efficiency of storage charge (from bus)
                    sto_ef_dis = p[:sto_ef_dis] # efficiency of storage discharge
                    feedincap = p[:feedincap]
                    max_sto_flow = p[:max_sto_flow] # relative cap of charge/discharge in one hour
                    event_per_scen = event_per_scen
                    # Euro
                    inv_budget = p[:inv_budget] # Make the problem bounded
                    max_pv = p[:max_pv]
                    max_wind = p[:max_wind]
                    regularize_lossy_flows = reg_lossy_flows
                end
                lifetime_factor = asset_lifetime * 365 * 24 / number_of_hours
                # Component units to be invested in, kWp
                @decision(model, u_pv >= 0)
                @decision(model, u_wind >= 0)
                @decision(model, u_storage >= 0)

                @constraint(model, u_pv <= max_pv)
                @constraint(model, u_wind <= max_wind)

                # Curtailment
                @decision(model, 0 <= pv_cur[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= wind_cur[t in 1:number_of_hours] <= debug_cap)
                @constraint(model, [t in 1:number_of_hours], pv_cur[t] <= u_pv * pv[t])
                @constraint(model, [t in 1:number_of_hours], wind_cur[t] <= u_wind * wind[t])

                # Grid connection
                @decision(model, 0 <= gci[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= gco[t in 1:number_of_hours] <= debug_cap)
                
                # Feedin cap
                @constraint(model, [t in 1:number_of_hours], gci[t] <= feedincap)

                # Storage model
                @decision(model, 0 <= sto_to_bus[t in 1:number_of_hours] <= debug_cap) # into the bus from storage
                @decision(model, 0 <= sto_from_bus[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= sto_soc[t in 1:number_of_hours] <= debug_cap)
                @constraint(model, [t in 1:number_of_hours-1], sto_soc[t+1] == sto_soc[t] + sto_from_bus[t] * sto_ef_ch - sto_to_bus[t] / sto_ef_dis)
                @constraint(model, [t in 1:number_of_hours], sto_soc[t] <= u_storage)
                @constraint(model, [t in 1:number_of_hours], sto_from_bus[t] <= max_sto_flow * u_storage)
                @constraint(model, [t in 1:number_of_hours], sto_to_bus[t] <= max_sto_flow * u_storage)
                # Start and end condition
                @constraint(model, sto_soc[1] == u_storage / 2)
                @constraint(model, sto_soc[number_of_hours] + sto_from_bus[number_of_hours] * sto_ef_ch- sto_to_bus[number_of_hours] / sto_ef_dis == sto_soc[1])

                # Heat model
                @decision(model, u_heatpump >= 0)
                @decision(model, u_heat_storage >= 0)
                @decision(model, 0 <= heat_sto_to_bus[t in 1:number_of_hours] <= debug_cap) # to the heat storage
                @decision(model, 0 <= heat_sto_from_bus[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= heat_sto_soc[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= flow_energy2heat[t in 1:number_of_hours] <= debug_cap)
                @constraint(model, [t in 1:number_of_hours-1], heat_sto_soc[t+1] == heat_sto_soc[t] + heat_sto_from_bus[t] * heat_eff - heat_sto_to_bus[t] / heat_eff)
                @constraint(model, [t in 1:number_of_hours], heat_sto_soc[t] <= u_heat_storage)
                @constraint(model, [t in 1:number_of_hours], flow_energy2heat[t] <= (1/COP)*u_heatpump)
                # Start and end condition
                @constraint(model, heat_sto_soc[1] == u_heat_storage / 2)
                @constraint(model, heat_sto_soc[number_of_hours] + heat_sto_from_bus[number_of_hours] * heat_eff - heat_sto_to_bus[number_of_hours] / heat_eff == heat_sto_soc[1])

                # Investment budget
                @constraint(model, inv_budget, c_pv * u_pv + c_wind * u_wind + c_storage * u_storage 
                + c_heat_storage * u_heat_storage + c_heatpump * u_heatpump <= inv_budget)

                # Energy balance
                @constraint(model, [t in 1:number_of_hours], 
                gci[t] - gco[t] + u_pv * pv[t] - pv_cur[t] + u_wind * wind[t] - wind_cur[t] - demand[t] + sto_to_bus[t] - sto_from_bus[t] - flow_energy2heat[t] == 0)
                
                # Heat balance
                @constraint(model, [t in 1:number_of_hours], -heatdemand[t] + heat_sto_to_bus[t] - heat_sto_from_bus[t] + COP*flow_energy2heat[t] - heat_losses*heat_sto_soc[t] == 0)
                if guaranteed_flex
                    @constraint(model, [t in 1:number_of_hours-1], pv_cur[t] + wind_cur[t] + sto_soc[t+1] + heat_sto_soc[t+1]/COP >= F_pos) # TODO change to flow_energy2heat[t]
                    @constraint(model, [t in 1:number_of_hours-1], pv_cur[t] + wind_cur[t] - pv[t]*u_pv - wind[t]*u_wind + sto_soc[t+1] - u_storage + (heat_sto_soc[t+1] - u_heat_storage)/COP <= F_neg)
                    # should include heatpump constraint, but we do not consider it here
                end
                # Investment costs ...
                # ... and background operational schedule
                @objective(model, Min, (u_pv * c_pv + u_wind * c_wind + u_storage * c_storage
                + u_heat_storage * c_heat_storage + u_heatpump * c_heatpump) / lifetime_factor
                + c_i * sum(gci) - c_o * sum(gco) + regularize_lossy_flows * (sum(heat_sto_from_bus) + sum(heat_sto_to_bus) + sum(sto_from_bus) + sum(sto_to_bus)))
            end
            @stage 2 begin
                @parameters begin
                    recovery_time = p[:recovery_time]
                    c_i = p[:c_i]
                    c_o = p[:c_o]
                    penalty = p[:penalty]
                    heat_losses = p[:heat_losses]
                    heat_eff = p[:heat_eff]
                    sto_ef_ch = p[:sto_ef_ch]
                    sto_ef_dis = p[:sto_ef_dis]
                    COP = p[:COP]
                    max_sto_flow = p[:max_sto_flow]
                    event_per_scen = event_per_scen
                    regularize_lossy_flows = reg_lossy_flows
                    feedincap = p[:feedincap]
                end
                @uncertain t_xi s_xi F_xi # t_xi the time of flexibility demand, s_xi - sign (±1 or 0)
                t_xi_final = t_xi + recovery_time
                # t_xi is the event time.
                # at t_xi + recovery_time all variables have to match again.

                # Curtailment
                @recourse(model, 0 <= pv_cur2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= wind_cur2[t in 1:1+recovery_time] <= debug_cap)
                @constraint(model, [t in 1:1+recovery_time], pv_cur2[t] <= u_pv * pv[t + t_xi - 1])
                @constraint(model, [t in 1:1+recovery_time], wind_cur2[t] <= u_wind * wind[t + t_xi - 1])

                @constraint(model, pv_cur2[1 + recovery_time] == pv_cur[t_xi + recovery_time])
                @constraint(model, wind_cur2[1 + recovery_time] == wind_cur[t_xi + recovery_time])

                # Grid connection
                @recourse(model, 0 <= gci2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= gco2[t in 1:1+recovery_time] <= debug_cap)

                # final time equality
                @constraint(model, gci2[1 + recovery_time] == gci[t_xi + recovery_time])
                @constraint(model, gco2[1 + recovery_time] == gco[t_xi + recovery_time])

                # Feedin cap
                @constraint(model, [t in 2:1+recovery_time], gci2[t] <= feedincap)
                
                # initial time equality is more complex
                # if we have strict flex 
                if ! regularized
                    @constraint(model, strict_flex_in, gci2[1] == gci[t_xi])
                    @constraint(model, strict_flex_out, gco2[1] == gco[t_xi])
                end

                ## Utility variables to linearize min|gci[t_xi]-gci2[1]|
                @recourse(model, 0 <= gi1 <= debug_cap)
                @recourse(model, 0 <= gi2 <= debug_cap)
                @constraint(model, gci2[1]-gci[t_xi] == gi1-gi2)
                @recourse(model, 0 <= go1 <= debug_cap)
                @recourse(model, 0 <= go2 <= debug_cap)
                #@recourse(model, penalty_taken)
                #@constraint(model, penalty_taken == gi1 + gi2 + go1 + go2)
                
                @constraint(model, gco2[1]-gco[t_xi] == go1-go2)

                # Storage model
                @recourse(model, 0 <= sto_to_bus2[t in 1:1+recovery_time] <= debug_cap) # into the bus from storage
                @recourse(model, 0 <= sto_from_bus2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= sto_soc2[t in 1:1+recovery_time] <= debug_cap)
                @constraint(model, [t in 1:1+recovery_time-1], sto_soc2[t+1] == sto_soc2[t] + sto_from_bus2[t] * sto_ef_ch - sto_to_bus2[t] / sto_ef_dis)
                @constraint(model, [t in 1:1+recovery_time], sto_soc2[t] <= u_storage)
                @constraint(model, [t in 1:1+recovery_time], sto_from_bus2[t] <= max_sto_flow * u_storage)
                @constraint(model, [t in 1:1+recovery_time], sto_to_bus2[t] <= max_sto_flow * u_storage)

                # Start and end condition 
                @constraint(model, sto_soc2[1] == sto_soc[t_xi])
                @constraint(model, sto_soc2[1 + recovery_time] == sto_soc[t_xi + recovery_time])
                @constraint(model,  sto_from_bus2[1 + recovery_time] == sto_from_bus[t_xi + recovery_time])
                @constraint(model,  sto_to_bus2[1 + recovery_time] == sto_to_bus[t_xi + recovery_time])

                # Heat model
                @recourse(model, 0 <= heat_sto_to_bus2[t in 1:1+recovery_time] <= debug_cap) # from the heat storage
                @recourse(model, 0 <= heat_sto_from_bus2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= heat_sto_soc2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= flow_energy2heat2[t in 1:1+recovery_time] <= debug_cap)
                @constraint(model, [t in 1:1+recovery_time-1], heat_sto_soc2[t+1] == heat_sto_soc2[t] + heat_sto_from_bus2[t] * heat_eff - heat_sto_to_bus2[t] / heat_eff)
                @constraint(model, [t in 1:1+recovery_time], heat_sto_soc2[t] <= u_heat_storage)
                @constraint(model, [t in 1:1+recovery_time], flow_energy2heat2[t] <= (1/COP)*u_heatpump)
                # Start and end condition
                @constraint(model, heat_sto_soc2[1] == heat_sto_soc[t_xi])
                @constraint(model, heat_sto_soc2[1 + recovery_time] == heat_sto_soc[t_xi + recovery_time])
                @constraint(model,  heat_sto_from_bus2[1 + recovery_time] == heat_sto_from_bus[t_xi + recovery_time])
                @constraint(model,  heat_sto_to_bus2[1 + recovery_time] == heat_sto_to_bus[t_xi + recovery_time])
                @constraint(model,  flow_energy2heat2[1 + recovery_time] == flow_energy2heat[t_xi + recovery_time])

                # Event energy balance
                # The storage and other fast acting components use the recourse variables here.
                # They provide the balance. Grid connection is not allowed, as we are suporting the grid here. 
                @constraint(model, gci2[1] - gco2[1] + u_pv * pv[t_xi] - pv_cur2[1] + u_wind * wind[t_xi] - wind_cur2[1]
                - demand[t_xi] + sto_to_bus2[1] - sto_from_bus2[1]
                - flow_energy2heat2[1] - F_xi * s_xi == 0) # TODO CHeck that our sign convention on positive and negative flexibility agrees with literature

                # Post event energy balance
                @constraint(model, [t in 2:1+recovery_time],
                gci2[t] - gco2[t]
                + u_pv * pv[t + t_xi - 1] - pv_cur2[t] + u_wind * wind[t + t_xi - 1] - wind_cur2[t]
                - demand[t + t_xi - 1] + sto_to_bus2[t] - sto_from_bus2[t] - flow_energy2heat2[t] == 0)

                # Heat balance
                @constraint(model, [t in 1:1+recovery_time], -heatdemand[t + t_xi - 1] + heat_sto_to_bus2[t] - heat_sto_from_bus2[t] + COP*flow_energy2heat2[t] - heat_losses*heat_sto_soc2[t] == 0)

                # The objective function is the difference between the adjusted schedule and the final schedule,
                # scaled by the number of events per year.
                # plus the penalty.
                @objective(model, Min, event_per_scen*(
                c_i * (sum(gci2) - sum(gci[t_xi:t_xi_final]))
                - c_o * (sum(gco2) - sum(gco[t_xi:t_xi_final]))
                + penalty * (gi1 + gi2) + penalty * (go1 + go2)
                - regularize_lossy_flows * (
                    sum(heat_sto_from_bus[t_xi:t_xi_final]) + 
                    sum(heat_sto_to_bus[t_xi:t_xi_final]) + 
                    sum(sto_from_bus[t_xi:t_xi_final]) + 
                    sum(sto_to_bus[t_xi:t_xi_final]))
                + regularize_lossy_flows * (
                    sum(heat_sto_from_bus2) + 
                    sum(heat_sto_to_bus2) + 
                    sum(sto_from_bus2) + 
                    sum(sto_to_bus2))
                ))
            end
        end
    elseif cap_constraint == "improved_heat"
        energy_system = @stochastic_model begin 
            @stage 1 begin
                @parameters begin
                    # Expected lifetime of components, years
                    asset_lifetime = p[:asset_lifetime]
                    # Costs in Euro/1kWp
                    c_pv = p[:c_pv]
                    c_wind = p[:c_wind]
                    c_storage = p[:c_storage]
                    c_i = p[:c_i]
                    c_o = p[:c_o]
                    c_heat_storage = p[:c_heat_storage]
                    c_heatpump = p[:c_heatpump]
                    COP = p[:COP]
                    heat_losses = p[:heat_losses]
                    heat_eff = p[:heat_eff]
                    sto_ef_ch = p[:sto_ef_ch] # efficiency of storage charge (from bus)
                    sto_ef_dis = p[:sto_ef_dis] # efficiency of storage discharge
                    feedincap = p[:feedincap]
                    max_sto_flow = p[:max_sto_flow] # relative cap of charge/discharge in one hour
                    event_per_scen = event_per_scen
                    # Euro
                    inv_budget = p[:inv_budget] # Make the problem bounded
                    max_pv = p[:max_pv]
                    max_wind = p[:max_wind]
                    regularize_lossy_flows = reg_lossy_flows
                end
                lifetime_factor = asset_lifetime * 365 * 24 / number_of_hours
                # Component units to be invested in, kWp
                @decision(model, u_pv >= 0)
                @decision(model, u_wind >= 0)
                @decision(model, u_storage >= 0)

                @constraint(model, u_pv <= max_pv)
                @constraint(model, u_wind <= max_wind)

                # Curtailment
                @decision(model, 0 <= pv_cur[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= wind_cur[t in 1:number_of_hours] <= debug_cap)
                @constraint(model, [t in 1:number_of_hours], pv_cur[t] <= u_pv * pv[t])
                @constraint(model, [t in 1:number_of_hours], wind_cur[t] <= u_wind * wind[t])

                # Grid connection
                @decision(model, 0 <= gci[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= gco[t in 1:number_of_hours] <= debug_cap)
                
                # Feedin cap
                @constraint(model, [t in 1:number_of_hours], gci[t] <= feedincap)

                # Storage model
                @decision(model, 0 <= sto_to_bus[t in 1:number_of_hours] <= debug_cap) # into the bus from storage
                @decision(model, 0 <= sto_from_bus[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= sto_soc[t in 1:number_of_hours] <= debug_cap)
                @constraint(model, [t in 1:number_of_hours-1], sto_soc[t+1] == sto_soc[t] + sto_from_bus[t] * sto_ef_ch - sto_to_bus[t] / sto_ef_dis)
                @constraint(model, [t in 1:number_of_hours], sto_soc[t] <= u_storage)
                @constraint(model, [t in 1:number_of_hours], sto_from_bus[t] <= max_sto_flow * u_storage)
                @constraint(model, [t in 1:number_of_hours], sto_to_bus[t] <= max_sto_flow * u_storage)
                # Start and end condition
                @constraint(model, sto_soc[1] == u_storage / 2)
                @constraint(model, sto_soc[number_of_hours] + sto_from_bus[number_of_hours] * sto_ef_ch- sto_to_bus[number_of_hours] / sto_ef_dis == sto_soc[1])

                # Heat model
                @decision(model, u_heatpump >= 0)
                @decision(model, u_heat_storage >= 0)
                @decision(model, 0 <= heat_sto_to_bus[t in 1:number_of_hours] <= debug_cap) # to the heat storage
                @decision(model, 0 <= heat_sto_from_bus[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= heat_sto_soc[t in 1:number_of_hours] <= debug_cap)
                @decision(model, 0 <= flow_energy2heat[t in 1:number_of_hours] <= debug_cap)
                @constraint(model, [t in 1:number_of_hours-1], heat_sto_soc[t+1] == heat_sto_soc[t] + heat_sto_from_bus[t] * heat_eff - heat_sto_to_bus[t] / heat_eff)
                @constraint(model, [t in 1:number_of_hours], heat_sto_soc[t] <= u_heat_storage)
                @constraint(model, [t in 1:number_of_hours], flow_energy2heat[t] <= (1/COP)*u_heatpump)
                # Start and end condition
                @constraint(model, heat_sto_soc[1] == u_heat_storage / 2)
                @constraint(model, heat_sto_soc[number_of_hours] + heat_sto_from_bus[number_of_hours] * heat_eff - heat_sto_to_bus[number_of_hours] / heat_eff == heat_sto_soc[1])

                # Investment budget
                @constraint(model, inv_budget, c_pv * u_pv + c_wind * u_wind + c_storage * u_storage 
                + c_heat_storage * u_heat_storage + c_heatpump * u_heatpump <= inv_budget)

                # Energy balance
                @constraint(model, [t in 1:number_of_hours], 
                gci[t] - gco[t] + u_pv * pv[t] - pv_cur[t] + u_wind * wind[t] - wind_cur[t] - demand[t] + sto_to_bus[t] - sto_from_bus[t] - flow_energy2heat[t] == 0)
                
                # Heat balance
                @constraint(model, [t in 1:number_of_hours], -heatdemand[t] + heat_sto_to_bus[t] - heat_sto_from_bus[t] + COP*flow_energy2heat[t] - heat_losses*heat_sto_soc[t] == 0)
                # Naive flex potential constraint
                # for the positive constraint it's easy: we stop pumping energy into heat
                # for the negative, we define utility variables to linearize max(flow_energy2heat[t]-u_heatpump, (heat_sto_soc[t+1]-u_heat_storage)/COP)
                @decision(model, -debug_cap*100. <= hc[t in 1:number_of_hours-1] <= 0.) # this variable will be the minimum
                @constraint(model, [t in 1:number_of_hours-1], hc[t] >= flow_energy2heat[t]-u_heatpump)
                @constraint(model, [t in 1:number_of_hours-1], hc[t] >= (heat_sto_soc[t+1]-u_heat_storage)/COP)
                @decision(model, hc_decision[t in 1:number_of_hours-1], binary = true)
                #set_binary(hc_decision)
                @constraint(model, [t in 1:number_of_hours-1], flow_energy2heat[t]-u_heatpump-(heat_sto_soc[t+1]-u_heat_storage)/COP<=debug_cap*100*hc_decision[t])
                @constraint(model, [t in 1:number_of_hours-1], -flow_energy2heat[t]+u_heatpump+(heat_sto_soc[t+1]-u_heat_storage)/COP<=debug_cap*100*(1-hc_decision[t]))
                @constraint(model, [t in 1:number_of_hours-1], hc[t] <= flow_energy2heat[t]-u_heatpump + debug_cap*100. *(1-hc_decision[t]))
                @constraint(model, [t in 1:number_of_hours-1], hc[t] <= (heat_sto_soc[t+1]-u_heat_storage)/COP + debug_cap*100. * hc_decision[t])
                # these represent constraints by heat pump transport and heat storage capacity
                # we need two new dummy variables to linearize the minimum of the two
                if guaranteed_flex
                    @constraint(model, [t in 1:number_of_hours-1], pv_cur[t] + wind_cur[t] + sto_soc[t+1] + flow_energy2heat[t] >= F_pos)
                    @constraint(model, [t in 1:number_of_hours-1], pv_cur[t] + wind_cur[t] - pv[t]*u_pv - wind[t]*u_wind + sto_soc[t+1] - u_storage + hc[t] <= F_neg)
                    # should include heatpump constraint, but we do not consider it here
                end
                # Investment costs ...
                # ... and background operational schedule
                @objective(model, Min, (u_pv * c_pv + u_wind * c_wind + u_storage * c_storage
                + u_heat_storage * c_heat_storage + u_heatpump * c_heatpump) / lifetime_factor
                + c_i * sum(gci) - c_o * sum(gco) + regularize_lossy_flows * (sum(heat_sto_from_bus) + sum(heat_sto_to_bus) + sum(sto_from_bus) + sum(sto_to_bus)))
            end
            @stage 2 begin
                @parameters begin
                    recovery_time = p[:recovery_time]
                    c_i = p[:c_i]
                    c_o = p[:c_o]
                    penalty = p[:penalty]
                    heat_losses = p[:heat_losses]
                    heat_eff = p[:heat_eff]
                    sto_ef_ch = p[:sto_ef_ch]
                    sto_ef_dis = p[:sto_ef_dis]
                    COP = p[:COP]
                    max_sto_flow = p[:max_sto_flow]
                    event_per_scen = event_per_scen
                    regularize_lossy_flows = reg_lossy_flows
                    feedincap = p[:feedincap]
                end
                @uncertain t_xi s_xi F_xi # t_xi the time of flexibility demand, s_xi - sign (±1 or 0)
                t_xi_final = t_xi + recovery_time
                # t_xi is the event time.
                # at t_xi + recovery_time all variables have to match again.

                # Curtailment
                @recourse(model, 0 <= pv_cur2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= wind_cur2[t in 1:1+recovery_time] <= debug_cap)
                @constraint(model, [t in 1:1+recovery_time], pv_cur2[t] <= u_pv * pv[t + t_xi - 1])
                @constraint(model, [t in 1:1+recovery_time], wind_cur2[t] <= u_wind * wind[t + t_xi - 1])

                @constraint(model, pv_cur2[1 + recovery_time] == pv_cur[t_xi + recovery_time])
                @constraint(model, wind_cur2[1 + recovery_time] == wind_cur[t_xi + recovery_time])

                # Grid connection
                @recourse(model, 0 <= gci2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= gco2[t in 1:1+recovery_time] <= debug_cap)

                # final time equality
                @constraint(model, gci2[1 + recovery_time] == gci[t_xi + recovery_time])
                @constraint(model, gco2[1 + recovery_time] == gco[t_xi + recovery_time])

                # Feedin cap
                @constraint(model, [t in 2:1+recovery_time], gci2[t] <= feedincap)
                
                # initial time equality is more complex
                # if we have strict flex 
                if ! regularized
                    @constraint(model, strict_flex_in, gci2[1] == gci[t_xi])
                    @constraint(model, strict_flex_out, gco2[1] == gco[t_xi])
                end

                ## Utility variables to linearize min|gci[t_xi]-gci2[1]|
                @recourse(model, 0 <= gi1 <= debug_cap)
                @recourse(model, 0 <= gi2 <= debug_cap)
                @constraint(model, gci2[1]-gci[t_xi] == gi1-gi2)
                @recourse(model, 0 <= go1 <= debug_cap)
                @recourse(model, 0 <= go2 <= debug_cap)
                #@recourse(model, penalty_taken)
                #@constraint(model, penalty_taken == gi1 + gi2 + go1 + go2)
                
                @constraint(model, gco2[1]-gco[t_xi] == go1-go2)

                # Storage model
                @recourse(model, 0 <= sto_to_bus2[t in 1:1+recovery_time] <= debug_cap) # into the bus from storage
                @recourse(model, 0 <= sto_from_bus2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= sto_soc2[t in 1:1+recovery_time] <= debug_cap)
                @constraint(model, [t in 1:1+recovery_time-1], sto_soc2[t+1] == sto_soc2[t] + sto_from_bus2[t] * sto_ef_ch - sto_to_bus2[t] / sto_ef_dis)
                @constraint(model, [t in 1:1+recovery_time], sto_soc2[t] <= u_storage)
                @constraint(model, [t in 1:1+recovery_time], sto_from_bus2[t] <= max_sto_flow * u_storage)
                @constraint(model, [t in 1:1+recovery_time], sto_to_bus2[t] <= max_sto_flow * u_storage)

                # Start and end condition 
                @constraint(model, sto_soc2[1] == sto_soc[t_xi])
                @constraint(model, sto_soc2[1 + recovery_time] == sto_soc[t_xi + recovery_time])
                @constraint(model,  sto_from_bus2[1 + recovery_time] == sto_from_bus[t_xi + recovery_time])
                @constraint(model,  sto_to_bus2[1 + recovery_time] == sto_to_bus[t_xi + recovery_time])

                # Heat model
                @recourse(model, 0 <= heat_sto_to_bus2[t in 1:1+recovery_time] <= debug_cap) # from the heat storage
                @recourse(model, 0 <= heat_sto_from_bus2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= heat_sto_soc2[t in 1:1+recovery_time] <= debug_cap)
                @recourse(model, 0 <= flow_energy2heat2[t in 1:1+recovery_time] <= debug_cap)
                @constraint(model, [t in 1:1+recovery_time-1], heat_sto_soc2[t+1] == heat_sto_soc2[t] + heat_sto_from_bus2[t] * heat_eff - heat_sto_to_bus2[t] / heat_eff)
                @constraint(model, [t in 1:1+recovery_time], heat_sto_soc2[t] <= u_heat_storage)
                @constraint(model, [t in 1:1+recovery_time], flow_energy2heat2[t] <= (1/COP)*u_heatpump)
                # Start and end condition
                @constraint(model, heat_sto_soc2[1] == heat_sto_soc[t_xi])
                @constraint(model, heat_sto_soc2[1 + recovery_time] == heat_sto_soc[t_xi + recovery_time])
                @constraint(model,  heat_sto_from_bus2[1 + recovery_time] == heat_sto_from_bus[t_xi + recovery_time])
                @constraint(model,  heat_sto_to_bus2[1 + recovery_time] == heat_sto_to_bus[t_xi + recovery_time])
                @constraint(model,  flow_energy2heat2[1 + recovery_time] == flow_energy2heat[t_xi + recovery_time])

                # Event energy balance
                # The storage and other fast acting components use the recourse variables here.
                # They provide the balance. Grid connection is not allowed, as we are suporting the grid here. 
                @constraint(model, gci2[1] - gco2[1] + u_pv * pv[t_xi] - pv_cur2[1] + u_wind * wind[t_xi] - wind_cur2[1]
                - demand[t_xi] + sto_to_bus2[1] - sto_from_bus2[1]
                - flow_energy2heat2[1] - F_xi * s_xi == 0) # TODO CHeck that our sign convention on positive and negative flexibility agrees with literature

                # Post event energy balance
                @constraint(model, [t in 2:1+recovery_time],
                gci2[t] - gco2[t]
                + u_pv * pv[t + t_xi - 1] - pv_cur2[t] + u_wind * wind[t + t_xi - 1] - wind_cur2[t]
                - demand[t + t_xi - 1] + sto_to_bus2[t] - sto_from_bus2[t] - flow_energy2heat2[t] == 0)

                # Heat balance
                @constraint(model, [t in 1:1+recovery_time], -heatdemand[t + t_xi - 1] + heat_sto_to_bus2[t] - heat_sto_from_bus2[t] + COP*flow_energy2heat2[t] - heat_losses*heat_sto_soc2[t] == 0)

                # The objective function is the difference between the adjusted schedule and the final schedule,
                # scaled by the number of events per year.
                # plus the penalty.
                @objective(model, Min, event_per_scen*(
                c_i * (sum(gci2) - sum(gci[t_xi:t_xi_final]))
                - c_o * (sum(gco2) - sum(gco[t_xi:t_xi_final]))
                + penalty * (gi1 + gi2) + penalty * (go1 + go2)
                - regularize_lossy_flows * (
                    sum(heat_sto_from_bus[t_xi:t_xi_final]) + 
                    sum(heat_sto_to_bus[t_xi:t_xi_final]) + 
                    sum(sto_from_bus[t_xi:t_xi_final]) + 
                    sum(sto_to_bus[t_xi:t_xi_final]))
                + regularize_lossy_flows * (
                    sum(heat_sto_from_bus2) + 
                    sum(heat_sto_to_bus2) + 
                    sum(sto_from_bus2) + 
                    sum(sto_to_bus2))
                ))
            end
        end
    else
        println("Unknown constraint type!")
        energy_system = nothing
    end
    energy_system
end

inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
op_vars = [:gci, :gco, :sto_soc, :sto_to_bus, :sto_from_bus,
            :heat_sto_soc, :heat_sto_to_bus, :heat_sto_from_bus,
            :flow_energy2heat, :pv_cur, :wind_cur]
"""
Get decision variables associated with investment rather than system operation.
"""
get_investments(sp) = Dict((
    [var => value.(sp[1, var]) for var in inv_vars]
))

"""
Fix the investment variables.
"""
function fix_investment!(sp, investments)
    for (var_sym, value) in zip(keys(investments), values(investments))
        fix(decision_by_name(sp, 1, string(var_sym)), value)
    end
    println("The investments are fixed.")
end

function unfix_investment!(sp)
    for var_sym in inv_vars
        unfix(decision_by_name(sp, 1, string(var_sym)))
    end
    println("The investments are released")
end

"""
Get decision variables associated with 1st stage.
"""
get_operation(sp) = Dict((
    :pv_cur => value.(sp[1,:pv_cur]),
    :wind_cur => value.(sp[1,:wind_cur]),
    :gci => value.(sp[1,:gci]),
    :gco => value.(sp[1,:gco]),
    :sto_to_bus => value.(sp[1,:sto_to_bus]),
    :sto_from_bus => value.(sp[1,:sto_from_bus]),
    :sto_soc => value.(sp[1,:sto_soc]),
    :heat_sto_to_bus => value.(sp[1,:heat_sto_to_bus]),
    :heat_sto_from_bus => value.(sp[1,:heat_sto_from_bus]),
    :heat_sto_soc => value.(sp[1,:heat_sto_soc]),
    :flow_energy2heat => value.(sp[1,:flow_energy2heat])
))

"""
Fix the investment variables.
"""
function fix_operation!(sp, operation, number_of_hours)
    for (var_sym, value) in zip(keys(operation), values(operation))
        for i in 1:number_of_hours
            fix(decision_by_name(sp, 1, string(var_sym)*"[$i]"), value[i])
        end
    end
    println("Operational schedule is fixed")
end

function unfix_operation!(sp, number_of_hours)
    for var_sym in keys(op_vars)
        for i in 1:number_of_hours
            unfix(decision_by_name(sp, 1, string(var_sym)*"[$i]"))
        end
    end
    println("Operational schedule is released")
end

get_recovery(sp, n)= Dict((
    :pv_cur2 => value.(sp[2,:pv_cur2], n),
    :wind_cur2 => value.(sp[2,:wind_cur2], n),
    :gci2 => value.(sp[2,:gci2], n),
    :gco2 => value.(sp[2,:gco2], n),
    :sto_to_bus2 => value.(sp[2,:sto_to_bus], n),
    :sto_from_bus2 => value.(sp[2,:sto_from_bus2], n),
    :sto_soc2 => value.(sp[2,:sto_soc2], n),
    :heat_sto_to_bus2 => value.(sp[2,:heat_sto_to_bus2], n),
    :heat_sto_from_bus2 => value.(sp[2,:heat_sto_from_bus2], n),
    :heat_sto_soc2 => value.(sp[2,:heat_sto_soc2], n),
    :flow_energy2heat2 => value.(sp[2,:flow_energy2heat2], n)
))

"""
Get array of penalties taken in each scenario. 
0 means that scenario was satisfied without using grid connection.
"""
function get_penalty_array(sp; scens = nothing)
    if isnothing(scens)
        scens = scenarios(sp)
    end
    penalties_gci = []
    penalties_gco = []
    for i in eachindex(scens)
        t_xi = scenarios(sp)[i].data.t_xi
        push!(penalties_gci,value.(sp[2,:gci2],i)[1] - value.(sp[1,:gci])[t_xi])
        push!(penalties_gco,value.(sp[2,:gco2],i)[1] - value.(sp[1,:gco])[t_xi])
    end
    return penalties_gci, penalties_gco
end

function get_penalized_scenarios(sp; scens = nothing)
    if isnothing(scens)
        scens = scenarios(sp)
    end
    penalized = []
    for i in eachindex(scens)
        t_xi = scenarios(sp)[i].data.t_xi
        if abs(value.(sp[2,:gci2],i)[1] - value.(sp[1,:gci])[t_xi]) >= 1e-9 || abs(value.(sp[2,:gco2],i)[1] - value.(sp[1,:gco])[t_xi]) >= 1e-9
        #if Bool(value.(sp[2,:penalty_taken],i))
            push!(penalized,i)
        end
    end
    return penalized
end

function get_total_investment(sp; number_of_hours = 24*365)    
    total_inv = sp.stages[1].parameters[:c_pv]*value.(sp[1, :u_pv]) + 
                sp.stages[1].parameters[:c_wind]*value.(sp[1, :u_wind]) +
                sp.stages[1].parameters[:c_storage]*value.(sp[1, :u_storage]) +
                sp.stages[1].parameters[:c_heat_storage]*value.(sp[1, :u_heat_storage]) +
                sp.stages[1].parameters[:c_heatpump]*value.(sp[1, :u_heatpump])
    #lifetime_factor = sp.stages[1].parameters[:asset_lifetime] * 365 * 24 / number_of_hours
    return total_inv
end

function get_total_investment(data_dict::Dict{String, Any}; number_of_hours = 24*365)    
    total_inv = sum([data_dict["params"]["c_"*var]]*data_dict["inv"]["u_"*var] 
        for var in ["pv","wind","storage","heat_storage","heatpump"])
    lifetime_factor = data_dict["params"]["asset_lifetime"] * 365 * 24 / number_of_hours
    return total_inv[1]
end

function get_total_investment(data_dict::Dict{Symbol, Any}; number_of_hours = 24*365)    
    total_inv = sum([data_dict[:params][Symbol("c_"*var)]]*data_dict[:inv][Symbol("u_"*var)] 
        for var in ["pv","wind","storage","heat_storage","heatpump"])
    lifetime_factor = data_dict[:params][:asset_lifetime]* 365 * 24 / number_of_hours
    return total_inv[1]
end
function get_operation_cost(sp)
    op_cost = sp.stages[1].parameters[:c_i]*sum(value.(sp[1, :gci])) - sp.stages[1].parameters[:c_o]*sum(value.(sp[1, :gco])) +
    sp.stages[1].parameters[:regularize_lossy_flows]*(sum(value.(sp[1, :heat_sto_from_bus]))+sum(value.(sp[1, :heat_sto_to_bus]))+
        sum(value.(sp[1, :sto_from_bus])) + sum(value.(sp[1, :sto_to_bus])))
    return op_cost
end
function get_operation_cost(data_dict::Dict{String, Any})
    op_cost = data_dict["params"]["c_i"]sum(data_dict["op"]["gci"]) - data_dict["params"]["c_o"]sum(data_dict["op"]["gco"]) +
        data_dict["params"]["regularize_lossy_flows"]*(sum(data_dict["op"]["heat_sto_from_bus"])+sum(data_dict["op"]["heat_sto_to_bus"])+
        sum(data_dict["op"]["sto_from_bus"]) + sum(data_dict["op"]["sto_to_bus"]))
    return op_cost
end
function get_operation_cost(data_dict::Dict{Symbol, Any})
    op_cost = data_dict[:params][:c_i]sum(data_dict[:op][:gci]) - data_dict[:params][:c_o]sum(data_dict[:op][:gco]) +
        data_dict[:params][:regularize_lossy_flows]*(sum(data_dict[:op][:heat_sto_from_bus])+sum(data_dict[:op][:heat_sto_to_bus])+
        sum(data_dict[:op][:sto_from_bus]) + sum(data_dict[:op][:sto_to_bus]))
    return op_cost
end
function get_servicing_cost(data_dict::Dict{Symbol, Any}; number_of_hours = 24*365)
    total_inv = get_total_investment(data_dict, number_of_hours=number_of_hours)
    op_cost = get_operation_cost(data_dict)
    return data_dict[:cost]-total_inv-op_cost
end

function get_servicing_cost(sp; number_of_hours = 24*365)
    total_inv = get_total_investment(sp, number_of_hours = number_of_hours)
    op_cost = get_operation_cost(sp)
    return objective_value(sp)-total_inv-op_cost
end


function get_all_data(sp; rec = true, scen = true)
    inv_d = get_investments(sp)
    op_d = get_operation(sp)
    n = length(scenarios(sp))
    if rec
        rec_d = [get_recovery(sp, i) for i in 1:n]
    else rec_d = []
    end
    if scen
        scen_d = [Dict((:probability => s.probability, :t_xi => s.data[:t_xi], 
        :F_xi => s.data[:F_xi], :s_xi => s.data[:s_xi])) for s in scenarios(sp)]
    else 
        scen_d = []
    end
    params = merge(sp.stages[1].parameters, sp.stages[2].parameters)
    return Dict((:inv => inv_d, :op => op_d, :rec => rec_d, :scen => scen_d, :params => params))
end

"""This function reconstructs the sample of scenarios back from the format saved to JSON,
where dict keys are replaced with strings instead of symbols
"""
function reconstruct_sample(scens)
    return [@scenario t_xi = s["t_xi"] s_xi = s["s_xi"] F_xi = s["F_xi"] probability = s["probability"]["π"] for s in scens]
end