using StochasticPrograms
using Random

##
"""
Default values of system parameters.
"""
default_es_pars = Dict((
    :c_i => .3,
    :c_o => .05,
    :c_sto_op => 0.0001,
    :asset_lifetime => 10.,
    :c_pv => 700.,
    :c_wind => 2000.,
    :c_storage => 600.,
    :c_heat_storage => 30.,
    :c_heatpump => 20.,
    :inv_budget => 500000000.,
    :recovery_time => 72,
    :COP => 0.8,
    :heat_losses => 0.01
    :penalty => 10000.
))

"""
Sample n scenarios.
"""
function simple_flex_sampler(n, F_max, t_max)
    [@scenario t_xi = rand(1:t_max) s_xi = rand([-1, 1]) F_xi = rand() * F_max probability = 1/n 
        for i in 1:n]
end

"""
Get an array with a single pseudo scenario with no flexiblity. This is useful for optimizing the system as if no flexibilty was introduced.
"""
function no_flex_pseudo_sampler()
    [@scenario t_xi = 1 s_xi = 1 F_xi = 0. probability = 1.
    ,]
end

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
    c_sto_op = p[:c_sto_op]
    c_i = p[:c_i]
    c_o = p[:c_o]
    recovery_time = p[:recovery_time]
    COP = p[:COP]
    heat_losses = p[:heat_losses]
    energy_system = @stochastic_model begin 
        @stage 1 begin
            @parameters begin
                # Expected lifetime of components, years
                asset_lifetime = p[:asset_lifetime]
                # Costs in Euro/1kWp
                c_pv = p[:c_pv]
                c_wind = p[:c_wind]
                c_storage = p[:c_storage]
                c_sto_op = c_sto_op
                c_i = c_i
                c_o = c_o
                c_heat_storage = p[:c_heat_storage]
                c_heatpump = p[:c_heatpump]
                COP = COP
                heat_losses = heat_losses
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
            # Storage model
            @decision(model, sto_in[t in 1:number_of_hours] >= 0) # into the bus from storage
            @decision(model, sto_out[t in 1:number_of_hours] >= 0)
            @decision(model, sto_soc[t in 1:number_of_hours] >= 0)
            @constraint(model, [t in 1:number_of_hours-1], sto_soc[t+1] == sto_soc[t] - sto_out[t] + sto_in[t])
            @constraint(model, [t in 1:number_of_hours], sto_soc[t] <= u_storage)
            @constraint(model, sto_soc[1] == u_storage / 2)
            @constraint(model, sto_soc[number_of_hours] - sto_out[number_of_hours] + sto_in[number_of_hours] == sto_soc[1])
            # Heat model
            @decision(model, u_heatpump >= 0)
            @decision(model, u_heat_storage >= 0)
            @decision(model, heat_sto_in[t in 1:number_of_hours] >= 0) # from the heat storage
            @decision(model, heat_sto_out[t in 1:number_of_hours] >= 0)
            @decision(model, heat_sto_soc[t in 1:number_of_hours] >= 0)
            @decision(model, heatpumpflow[t in 1:number_of_hours] >= 0)

            @constraint(model, [t in 1:number_of_hours-1], heat_sto_soc[t+1] == heat_sto_soc[t] - heat_sto_out[t] + heat_sto_in[t])
            @constraint(model, [t in 1:number_of_hours], heat_sto_soc[t] <= u_heat_storage)
            @constraint(model, heat_sto_soc[1] == u_heat_storage / 2)
            @constraint(model, heat_sto_soc[number_of_hours] - heat_sto_out[number_of_hours] + heat_sto_in[number_of_hours] == heat_sto_soc[1])
            # Energy balance
            @constraint(model, [t in 1:number_of_hours], 
            gci[t] - gco[t] + u_pv * pv[t] + u_wind * wind[t] - demand[t] + sto_in[t] - sto_out[t] - heatpumpflow[t] == 0)
            # Heat balance
            @constraint(model, [t in 1:number_of_hours], -heatdemand[t] - heat_sto_in[t] + heat_sto_out[t] + COP*heatpumpflow[t] - heat_losses*heat_sto_soc[t] == 0)
            # Investment costs
            @objective(model, Min, (u_pv * c_pv + u_wind * c_wind + u_storage * c_storage
            + u_heat_storage * c_heat_storage + u_heatpump * c_heatpump) / lifetime_factor
            + c_i * sum(gci) - c_o * sum(gco) +
            c_sto_op * sum(sto_in) + c_sto_op * sum(sto_out))
        end
        @stage 2 begin
            @parameters begin
                recovery_time = recovery_time
                c_sto_op = c_sto_op
                c_i = c_i
                c_o = c_o
                penalty = p[:penalty]
                heat_losses = heat_losses
                COP = COP
            end
            @uncertain t_xi s_xi F_xi # t_xi the time of flexibility demand, s_xi - sign (±1 or 0)
            t_xi_final = t_xi + recovery_time - 1
            @known(model, u_pv)
            @known(model, u_wind)
            @known(model, u_storage)
            @known(model, gci)
            @known(model, gco)
            @known(model, sto_in)
            @known(model, sto_out)
            # Post event components
            # Grid connection
            @recourse(model, gci2[t in 1:recovery_time] >= 0)
            @recourse(model, gco2[t in 1:recovery_time] >= 0)
            # Storage model
            @recourse(model, sto_in2[t in 1:recovery_time] >= 0) # into the bus from storage
            @recourse(model, sto_out2[t in 1:recovery_time] >= 0)
            @recourse(model, sto_soc2[t in 1:recovery_time] >= 0)
            @constraint(model, [t in 1:recovery_time-1], sto_soc2[t+1] == sto_soc2[t] - sto_out2[t] + sto_in2[t])
            @constraint(model, [t in 1:recovery_time], sto_soc2[t] <= u_storage)
            @constraint(model, sto_soc2[1] == sto_soc[t_xi])
            @constraint(model, sto_soc2[recovery_time] - sto_out2[recovery_time] + sto_in2[recovery_time] == sto_soc[t_xi+recovery_time])
            # Heat model
            @recourse(model, heat_sto_in2[t in 1:recovery_time] >= 0) # from the heat storage
            @recourse(model, heat_sto_out2[t in 1:recovery_time] >= 0)
            @recourse(model, heat_sto_soc2[t in 1:recovery_time] >= 0)
            @recourse(model, heatpumpflow2[t in 1:recovery_time] >= 0)
            @constraint(model, [t in 1:recovery_time-1], heat_sto_soc2[t+1] == heat_sto_soc2[t] - heat_sto_out2[t] + heat_sto_in2[t])
            @constraint(model, [t in 1:recovery_time], heat_sto_soc2[t] <= u_heat_storage)
            @constraint(model, heat_sto_soc2[1] == heat_sto_soc[t_xi])
            @constraint(model, heat_sto_soc2[recovery_time] - heat_sto_out2[recovery_time] + heat_sto_in2[recovery_time] == heat_sto_soc[t_xi+recovery_time])

            # Event energy balance
            # The storage and other fast acting components use the recourse variables here.
            # They provide the balance. Grid connection is not allowed, as we are suporting the grid here. 
            @constraint(model, gci2[1] - gco2[1] + u_pv * pv[t_xi] + u_wind * wind[t_xi]
             - demand[t_xi] + sto_in2[1] - sto_out2[1]
             - heatpumpflow2[1] - F_xi * s_xi == 0) # TODO CHeck that our sign convention on positive and negative flexibility agrees with literature
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
            - demand[t + t_xi - 1] + sto_in2[t] - sto_out2[t] - heatpumpflow2[t] == 0)
            # Heat balance
            @constraint(model, [t in 1:recovery_time], -heatdemand[t] - heat_sto_in2[t] + heat_sto_out2[t] + COP*heatpumpflow2[t] - heat_losses*heat_sto_soc2[t] == 0)

            @objective(model, Min,
            + c_i * (sum(gci2) - sum(gci[t_xi:t_xi_final]))
            - c_o * (sum(gco2) - sum(gco[t_xi:t_xi_final]))
            + penalty * (gi1 + gi2) + penalty * (go1 + go2)
            + c_sto_op * (sum(sto_in2) + sum(sto_out2) 
            - sum(sto_in[t_xi:t_xi_final]) - sum(sto_out[t_xi:t_xi_final])))
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