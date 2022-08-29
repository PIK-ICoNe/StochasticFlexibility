using StochasticPrograms
using Random

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
    :recovery_time => 72,
    :COP => 3.5,
    :heat_losses => 0.2,
    :sto_ef_ch => 0.95,
    :sto_ef_dis => 0.95,
    :storage_losses => 0.03,
    :penalty => 10000.,
    :feedincap => 1000000.,
    :scens_in_year => 1.

))

"""
Draw a sample of n scenarios
"""
function simple_flex_sampler(n, F_max, t_max)
    [@scenario t_xi = rand(1:t_max) s_xi = rand([-1, 1]) F_xi = rand() * F_max probability = 1/n 
        for i in 1:n]
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

#=
function gap_sampler(indices)
    [@scenario t_xi = 1 s_xi = 1 F_xi = 0. probability = 1.
     for i in indices]
end
=#

"""
Define energy system.
Parameters:
- pv, wind - weather timeseries
- demand, heatdemand - demand timeseries
- p - dictionary with system parameters, such as component costs, losses and recovery time window
- strict_flex - bool, if false, finite penalty is used
"""
function define_energy_system(pv, wind, demand, heatdemand; p = default_es_pars, strict_flex=false)
    number_of_hours = minimum([length(pv), length(demand), length(wind)])
    c_i = p[:c_i]
    c_o = p[:c_o]
    recovery_time = p[:recovery_time]
    COP = p[:COP]
    heat_losses = p[:heat_losses]
    storage_losses = p[:storage_losses]
    sto_ef_ch = p[:sto_ef_ch] # efficiency of storage charge (from bus)
    sto_ef_dis = p[:sto_ef_dis] # efficiency of storage discharge
    max_sto_flow = 0.2 # relative cap of charge/discharge in one hour
    energy_system = @stochastic_model begin 
        @stage 1 begin
            @parameters begin
                # Expected lifetime of components, years
                asset_lifetime = p[:asset_lifetime]
                # Costs in Euro/1kWp
                c_pv = p[:c_pv]
                c_wind = p[:c_wind]
                c_storage = p[:c_storage]
                c_i = c_i
                c_o = c_o
                c_heat_storage = p[:c_heat_storage]
                c_heatpump = p[:c_heatpump]
                COP = COP
                heat_losses = heat_losses
                storage_losses = storage_losses
                sto_ef_ch = sto_ef_ch
                sto_ef_dis = sto_ef_dis
                feedincap = p[:feedincap]
                # Euro
                inv_budget = p[:inv_budget] # Make the problem bounded
            end
            lifetime_factor = asset_lifetime * 365 * 24 / number_of_hours
            # Component units to be invested in, kWp
            @decision(model, u_pv >= 0)
            @decision(model, u_wind >= 0)
            @decision(model, u_storage >= 0)
            @constraint(model, c_pv * u_pv + c_wind * u_wind + c_storage * u_storage <= inv_budget)
            # Grid connection
            @decision(model, gci[t in 1:number_of_hours] >= 0)
            @decision(model, gco[t in 1:number_of_hours] >= 0)
            @constraint(model, sum(gco) <= feedincap)
            # Storage model
            @decision(model, sto_to_bus[t in 1:number_of_hours] >= 0) # into the bus from storage
            @decision(model, sto_from_bus[t in 1:number_of_hours] >= 0)
            @decision(model, sto_soc[t in 1:number_of_hours] >= 0)
            @constraint(model, [t in 1:number_of_hours], sto_from_bus[t] <= max_sto_flow*u_storage)
            @constraint(model, [t in 1:number_of_hours], sto_to_bus[t] <= max_sto_flow*u_storage)
            @constraint(model, [t in 1:number_of_hours-1], sto_soc[t+1] == sto_soc[t] + sto_from_bus[t] - sto_to_bus[t])
            @constraint(model, [t in 1:number_of_hours], sto_soc[t] <= u_storage)
            @constraint(model, sto_soc[1] == u_storage / 2)
            @constraint(model, sto_soc[number_of_hours] + sto_from_bus[number_of_hours] - sto_to_bus[number_of_hours] == sto_soc[1])
            # Heat model
            @decision(model, u_heatpump >= 0)
            @decision(model, u_heat_storage >= 0)
            @decision(model, heat_sto_to_bus[t in 1:number_of_hours] >= 0) # to the heat storage
            @decision(model, heat_sto_from_bus[t in 1:number_of_hours] >= 0)
            @decision(model, heat_sto_soc[t in 1:number_of_hours] >= 0)
            @decision(model, flow_energy2heat[t in 1:number_of_hours] >= 0)
            if maximum(heatdemand) > 0.
                @constraint(model, [t in 1:number_of_hours-1], heat_sto_soc[t+1] == heat_sto_soc[t] + heat_sto_from_bus[t] - heat_sto_to_bus[t])
                @constraint(model, [t in 1:number_of_hours], heat_sto_soc[t] <= u_heat_storage)
                @constraint(model, heat_sto_soc[1] == u_heat_storage / 2)
                @constraint(model, heat_sto_soc[number_of_hours] + heat_sto_from_bus[number_of_hours] - heat_sto_to_bus[number_of_hours] == heat_sto_soc[1])
                @constraint(model, [t in 1:number_of_hours], flow_energy2heat[t] <= 1/COP*u_heatpump)
                # Heat balance
                @constraint(model, [t in 1:number_of_hours], -heatdemand[t] + heat_sto_to_bus[t] - heat_sto_from_bus[t] + COP*flow_energy2heat[t] - heat_losses*heat_sto_soc[t] == 0)
            else
                @constraint(model, u_heat_storage == 0.)
                @constraint(model, u_heatpump == 0.)
                @constraint(model, [t in 1:number_of_hours], flow_energy2heat[t] == 0.)
            end
            # Energy balance
            @constraint(model, [t in 1:number_of_hours], 
            gci[t] - gco[t] + u_pv * pv[t] + u_wind * wind[t] - demand[t] + sto_to_bus[t]*sto_ef_dis - sto_from_bus[t]/sto_ef_ch - flow_energy2heat[t] == 0)

            # Investment costs ...
            # ... and background operational schedule
            @objective(model, Min, (u_pv * c_pv + u_wind * c_wind + u_storage * c_storage
            + u_heat_storage * c_heat_storage + u_heatpump * c_heatpump) / lifetime_factor
            + c_i * sum(gci) - c_o * sum(gco))
        end
        @stage 2 begin
            @parameters begin
                recovery_time = recovery_time
                c_i = c_i
                c_o = c_o
                penalty = p[:penalty]
                heat_losses = heat_losses
                COP = COP
                storage_losses = storage_losses
                sto_ef_ch = sto_ef_ch
                sto_ef_dis = sto_ef_dis
                scens_in_year = p[:scens_in_year]
            end
            @uncertain t_xi s_xi F_xi # t_xi the time of flexibility demand, s_xi - sign (±1 or 0)
            t_xi_final = t_xi + recovery_time - 1
            @known(model, u_pv)
            @known(model, u_wind)
            @known(model, u_storage)
            @known(model, gci)
            @known(model, gco)
            @known(model, sto_to_bus)
            @known(model, sto_from_bus)
            # Post event components
            # Grid connection
            @recourse(model, gci2[t in 1:recovery_time] >= 0)
            @recourse(model, gco2[t in 1:recovery_time] >= 0)
            # Storage model
            @recourse(model, sto_to_bus2[t in 1:recovery_time] >= 0) # into the bus from storage
            @recourse(model, sto_from_bus2[t in 1:recovery_time] >= 0)
            @recourse(model, sto_soc2[t in 1:recovery_time] >= 0)
            @constraint(model, [t in 1:recovery_time], sto_from_bus2[t] <= max_sto_flow*u_storage)
            @constraint(model, [t in 1:recovery_time], sto_to_bus2[t] <= max_sto_flow*u_storage)
            @constraint(model, [t in 1:recovery_time-1], sto_soc2[t+1] == sto_soc2[t] + sto_from_bus2[t] - sto_to_bus2[t])
            @constraint(model, [t in 1:recovery_time], sto_soc2[t] <= u_storage)
            @constraint(model, sto_soc2[1] == sto_soc[t_xi])
            @constraint(model, sto_soc2[recovery_time] + sto_from_bus2[recovery_time] - sto_to_bus2[recovery_time] == sto_soc[t_xi+recovery_time])
            # Heat model
            @recourse(model, heat_sto_to_bus2[t in 1:recovery_time] >= 0) # from the heat storage
            @recourse(model, heat_sto_from_bus2[t in 1:recovery_time] >= 0)
            @recourse(model, heat_sto_soc2[t in 1:recovery_time] >= 0)
            @recourse(model, flow_energy2heat2[t in 1:recovery_time] >= 0)
            if maximum(heatdemand) > 0.
                @constraint(model, [t in 1:recovery_time-1], heat_sto_soc2[t+1] == heat_sto_soc2[t] + heat_sto_from_bus2[t] - heat_sto_to_bus2[t])
                @constraint(model, [t in 1:recovery_time], heat_sto_soc2[t] <= u_heat_storage)
                @constraint(model, heat_sto_soc2[1] == heat_sto_soc[t_xi])
                @constraint(model, heat_sto_soc2[recovery_time] + heat_sto_from_bus2[recovery_time] - heat_sto_to_bus2[recovery_time] == heat_sto_soc[t_xi+recovery_time])
                @constraint(model, [t in 1:recovery_time], flow_energy2heat2[t] <= 1/COP*u_heatpump)
                # Heat balance
                @constraint(model, [t in 1:recovery_time], -heatdemand[t + t_xi - 1] + heat_sto_to_bus2[t] - heat_sto_from_bus2[t] + COP*flow_energy2heat2[t] - heat_losses*heat_sto_soc2[t] == 0)
            else
                @constraint(model, [t in 1:recovery_time], flow_energy2heat2[t] == 0.)
            end
            # Event energy balance
            # The storage and other fast acting components use the recourse variables here.
            # They provide the balance. Grid connection is not allowed, as we are suporting the grid here. 
            @constraint(model, gci2[1] - gco2[1] + u_pv * pv[t_xi] + u_wind * wind[t_xi]
             - demand[t_xi] + sto_to_bus2[1]*sto_ef_dis - sto_from_bus2[1]/sto_ef_ch
             - flow_energy2heat2[1] - F_xi * s_xi == 0) # TODO CHeck that our sign convention on positive and negative flexibility agrees with literature
            if strict_flex
                @constraint(model, gci2[1] == gci[t_xi])
                @constraint(model, gco2[1] == gco[t_xi])
            end
            ## Utility variables to linearize min|gci[t_xi]-gci2[1]|
            @recourse(model, gi1>=0)
            @recourse(model, gi2>=0)
            @constraint(model, gci2[1]-gci[t_xi] == gi1-gi2)
            @recourse(model, go1>=0)
            @recourse(model, go2>=0)
            @constraint(model, gco2[1]-gco[t_xi] == go1-go2)
            # Post event energy balance
            @constraint(model, [t in 2:recovery_time],
            gci2[t] - gco2[t]
            + u_pv * pv[t + t_xi - 1] + u_wind * wind[t + t_xi - 1]
            - demand[t + t_xi - 1] + sto_to_bus2[t]*sto_ef_dis - sto_from_bus2[t]/sto_ef_ch - flow_energy2heat2[t] == 0)
            # The objective function is the difference between the adjusted schedule and the final schedule
            # plus the penalty. TODO: We only evaluate the cost of individual events, so we should multiply the expectation
            # value with the expected number of events, e.g. one per week.
            # I tried to implement this by scaling the objective function here, but that made the problem unbounded.
            @objective(model, Min, scens_in_year*(
              c_i * (sum(gci2) - sum(gci[t_xi:t_xi_final]))
            - c_o * (sum(gco2) - sum(gco[t_xi:t_xi_final]))
            + penalty * (gi1 + gi2) + penalty * (go1 + go2)))
        end
    end
    energy_system
end

"""
Get decision variables associated with investment rather than system operation.
"""
get_investments(sp) = Dict((
    :u_pv => value.(sp[1, :u_pv]),
    :u_wind => value.(sp[1, :u_wind]),
    :u_storage => value.(sp[1, :u_storage]),
    :u_heatpump => value.(sp[1, :u_heatpump]),
    :u_heat_storage => value.(sp[1, :u_heat_storage])
))

"""
Fix the investment variables.
"""
function fix_investment!(sp, investments)
    for (var_sym, value) in zip(keys(investments), values(investments))
        fix(decision_by_name(sp, 1, string(var_sym)), value)
    end
end

function unfix_investment!(sp, investments)
    for var_sym in keys(investments)
        unfix(decision_by_name(sp, 1, string(var_sym)))
    end
end