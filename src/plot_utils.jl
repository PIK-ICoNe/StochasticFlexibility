using Plots;gr()
using SankeyPlots

function plot_results(sp_data::Dict{String, Any}, pv, w, el_d; plot_window = 1:length(pv), s=1, el_balance_vars=["gci", "gco"], storage_vars=[])
    u_pv = sp_data["inv"]["u_pv"]
    u_wind = sp_data["inv"]["u_wind"]
    plt_sto = plot(; legend = :outertopright)
    plt_invest = plot(; legend = :outertopright)
    plt = plot(; legend = :outertopright)
    plot!(plt_invest, plot_window, pv[plot_window] .* u_pv, label="PV")
    plot!(plt_invest, plot_window, w[plot_window] .* u_wind, label="Wind")
    plot!(plt_invest, plot_window, el_d[plot_window], label="Electrical demand")
    t_xi = sp_data["scen"][s]["t_xi"]
    recovery_time = sp_data["params"]["recovery_time"]

    stor_charge = sp_data["op"]["sto_soc"]
    plot!(plt_sto, plot_window, stor_charge[plot_window], label="global storage charge")
    if t_xi in plot_window
        plot!(plt_sto, (t_xi):(t_xi+recovery_time), sp_data["rec"][s]["sto_soc2"], label=string("stochastic storage charge")*string(s), linestyle=:dash, linewidth=2)
    end

    for var in el_balance_vars
        plot!(plt, plot_window, sp_data["op"][var][plot_window], label=var)
        if t_xi in plot_window
            var2 = var*"2"
            plot!(plt, (t_xi+1):(t_xi+recovery_time + 1), sp_data["rec"][s][var2], label=var*" 2nd stage, s = "*string(s), linestyle=:dash, linewidth=2)
        end
    end
    return plot(plt_invest, plt, plt_sto, layout = (3,1))
end

function plot_results(sp_data::Dict{Symbol, Any}, pv, w, el_d; plot_window = 1:length(pv), s=1, el_balance_vars=[:gci, :gco], storage_vars=[])
    u_pv = sp_data[:inv][:u_pv]
    u_wind = sp_data[:inv][:u_wind]
    plt_sto = plot(; legend = :outertopright)
    plt_invest = plot(; legend = :outertopright)
    plt = plot(; legend = :outertopright)
    plot!(plt_invest, plot_window, pv[plot_window] .* u_pv, label="PV")
    plot!(plt_invest, plot_window, w[plot_window] .* u_wind, label="Wind")
    plot!(plt_invest, plot_window, el_d[plot_window], label="Electrical demand")
    t_xi = sp_data[:scen][s][:t_xi]
    recovery_time = sp_data[:params][:recovery_time]

    stor_charge = sp_data[:op][:sto_soc]
    plot!(plt_sto, plot_window, stor_charge[plot_window], label="global storage charge")
    if t_xi in plot_window
        plot!(plt_sto, (t_xi):(t_xi+recovery_time), sp_data[:rec][s][:sto_soc2], label=string("stochastic storage charge")*string(s), linestyle=:dash, linewidth=2)
    end

    for var in el_balance_vars
        plot!(plt, plot_window, sp_data[:op][var][plot_window], label=string(var))
        if t_xi in plot_window
            var2 = string(var)*"2"
            plot!(plt, (t_xi+1):(t_xi+recovery_time + 1), sp_data[:rec][s][Symbol(var2)], label=string(var)*" 2nd stage, s = "*string(s), linestyle=:dash, linewidth=2)
        end
    end
    return plot(plt_invest, plt, plt_sto, layout = (3,1))
end

function plot_heat_layer(sp_data, heatdemand; plot_window = 1:length(heatdemand), s = 1)
    plt_heat = plot(; legend = :outertopright)
    COP = sp_data["params"]["COP"]

    plot!(plt_heat, plot_window, heatdemand[plot_window], label = "heat demand")
    plot!(plt_heat, plot_window, COP*sp_data["op"]["flow_energy2heat"][plot_window], label = "heatpump")
    plot!(plt_heat, plot_window, sp_data["op"]["heat_sto_soc"][plot_window], label = "heat storage SOC")
    plot!(plt_heat, plot_window, sp_data["op"]["heat_sto_to_bus"][plot_window]-sp_data["op"]["heat_sto_from_bus"][plot_window], label = "heat storage use")
    # plot!(plt_heat, (t_xi+1):(t_xi+1+recovery_time), COP*sp_data["rec"][s]["flow_energy2heat2"], label = "heatpump$s", linestyle=:dash, linewidth=2)
    # plot!(plt_heat, (t_xi+1):(t_xi+1+recovery_time), sp_data["rec"][s]["heat_sto_to_bus2"]-sp_data["rec"][2]["heat_sto_from_bus2"], label = "heat storage use $s", linestyle=:dash, linewidth=2)

end

function plot_raw_heat_layer(sp_data, heatdemand; plot_span = 1:length(heatdemand), s = 1)
    plt_heat = plot(; legend = :outertopright)
    COP = sp_data["params"]["COP"]

    #plot!(plt_heat, plot_span, heatdemand[plot_span], label = "heat demand")
    #plot!(plt_heat, plot_span, COP*value.(sp[1, :flow_energy2heat])[plot_span], label = "heatpump")
    plot!(plt_heat, plot_span, sp_data["op"]["heat_sto_soc"][plot_span], label = "heat storage SOC")
    plot!(plt_heat, plot_span, sp_data["op"]["heat_sto_to_bus"][plot_span], label = "heat storage to bus")
    plot!(plt_heat, plot_span, sp_data["op"]["heat_sto_from_bus"][plot_span], label = "heat storage from bus")
    # plot!(plt_heat, (t_xi+1):(t_xi+1+recovery_time), COP*value.(sp[2, :flow_energy2heat2], s), label = "heatpump$s", linestyle=:dash, linewidth=2)
    # plot!(plt_heat, (t_xi+1):(t_xi+1+recovery_time), value.(sp[2, :heat_sto_in2], s)-value.(sp[2, :heat_sto_out2], s), label = "heat storage use $s", linestyle=:dash, linewidth=2)

end

function plot_recovery_window_deviation(sp_data::Dict{String, Any}; s = 1)
    plt_gc = plot()
    plt_sto = plot()
    t_xi = sp_data["scen"][s]["t_xi"]
    recovery_time = sp_data["params"]["recovery_time"]
    plot!(plt_gc, COP*sp_data["op"]["flow_energy2heat"][t_xi+1:t_xi+1+recovery_time], label = "global flow_energy2heat")
    plot!(plt_gc, COP*sp_data["rec"][s]["flow_energy2heat2"], label = "stochastic flow_energy2heat")
    plot!(plt_gc, -sp_data["op"]["gci"][t_xi+1:t_xi+1+recovery_time]+sp_data["op"]["gco"][t_xi+1:t_xi+1+recovery_time], label = "global grid connection use")
    plot!(plt_gc, sp_data["rec"][s]["gco2"] - sp_data["rec"][s]["gci2"], xlabel = "time after the event, h", label = "stochastic grid connection use")
    #plot!(plt_gc, -value.(sp[1, :gco])[t_xi+1:t_xi+recovery_time]+value.(sp[2, :gco2], s), xlabel = "time after the event, h", label = "gco2-gco")
    
    stor_charge = sp_data["op"]["sto_soc"][t_xi+1:t_xi+1+recovery_time]
    stor_charge2 = sp_data["rec"][s]["sto_soc2"]
    
    plot!(plt_sto, stor_charge2-stor_charge, label = "soc2-soc")


    plot(plt_gc, plt_sto, layout = (2, 1))
end

function plot_recovery_window_deviation(sp; s = 1)
    plt_gc = plot()
    plt_sto = plot()
    t_xi = scenarios(sp)[s].data.t_xi
    recovery_time = sp.stages[2].parameters[:recovery_time]
    plot!(plt_gc, value.(sp[1, :flow_energy2heat])[t_xi+1:t_xi+1+recovery_time], label = "global flow_energy2heat")
    plot!(plt_gc, value.(sp[2, :flow_energy2heat2], s), label = "stochastic flow_energy2heat")
    plot!(plt_gc, -value.(sp[1, :gci])[t_xi+1:t_xi+1+recovery_time]+value.(sp[2, :gci2], s), label = "global grid connection use")
    plot!(plt_gc, value.(sp[2, :gco2], s) - value.(sp[1, :gco])[t_xi+1:t_xi+1+recovery_time], xlabel = "time after the event, h", label = "stochastic grid connection use")
    #plot!(plt_gc, -value.(sp[1, :gco])[t_xi+1:t_xi+recovery_time]+value.(sp[2, :gco2], s), xlabel = "time after the event, h", label = "gco2-gco")
    
    stor_charge = value.(sp[1, :sto_soc])[t_xi+1:t_xi+1+recovery_time]
    stor_charge2 = value.(sp[2, :sto_soc2], s)
    
    plot!(plt_sto, stor_charge2-stor_charge, label = "soc2-soc")

    plot(plt_gc, plt_sto, layout = (2, 1))
end

function plot_outcome(sp_base, t_xi, s_xi, F_xi; window_start=-2, window_end=2)
    scen = @scenario t_xi = t_xi s_xi = s_xi F_xi = F_xi probability = 1.
    sp = outcome_model(sp_base, optimal_decision(sp_base), scen; optimizer = subproblem_optimizer(sp_base))
    optimize!(sp)
    if termination_status(sp) != JuMP.MathOptInterface.OPTIMAL
        println("No optimum")
        return termination_status(sp)
    end

    recovery_window = t_xi:t_xi+length(sp[:gco2])-1
    plot_window = t_xi+window_start:t_xi+length(sp[:gco2])+window_end

    # Some of these should probably be shifted by ±1

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

function plot_scenario_debug(sp_data, s; window_start=-2, window_end=2, vars = ["gci", "gco", "sto_soc", "heat_sto_soc", "flow_energy2heat"])

    t_xi = sp_data["scen"][s]["t_xi"]
    recovery_time = sp_data["params"]["recovery_time"]

    recovery_window = t_xi:t_xi+recovery_time
    plot_window = t_xi+window_start:t_xi+recovery_time+window_end

    plots = []

    for var in vars
        plt_var = plot(title=var, legend=:outertopright)


        plot!(plt_var, plot_window, sp_data["op"][var][plot_window], label = "Base Case")
        plot!(plt_var, recovery_window, sp_data["rec"][s][var*"2"][plot_window], label = "Event")
    
        push!(plots, plt_var)
    end

    plot(plots...; layout=(length(plots),1), size = (800, 150*length(plots)))
end

function plot_base_case_raw(sp_data; plot_window = nothing, vars = ["gci", "gco", "sto_soc", "heat_sto_soc", "flow_energy2heat"])

    if isnothing(plot_window)
        plot_window = 1:length(sp_data["op"]["gci"])
    end

    plots = []

    for var in vars
        plt_var = plot(title=var, legend=:none)

        plot!(plt_var, plot_window, sp_data["op"][var][plot_window])
    
        push!(plots, plt_var)
    end

    plot(plots...; layout=(length(plots),1), size = (800, 150*length(plots)))
end


function plot_outcome_debug(sp_base, t_xi, s_xi, F_xi; window_start=-2, window_end=2, vars = ["gci", "gco", "sto_soc", "heat_sto_soc", "flow_energy2heat"])
    scen = @scenario t_xi = t_xi s_xi = s_xi F_xi = F_xi probability = 1.
    sp = outcome_model(sp_base, optimal_decision(sp_base), scen; optimizer = subproblem_optimizer(sp_base))
    optimize!(sp)
    if termination_status(sp) != JuMP.MathOptInterface.OPTIMAL
        println("No optimum")
        return termination_status(sp)
    end
    recovery_time = length(sp[:gco2])-1

    recovery_window = t_xi:t_xi+recovery_time
    plot_window = t_xi+window_start:t_xi+recovery_time+window_end

    plots = []

    for var in vars
        plt_var = plot(title=var, legend=:outertopright)
        symbol_base_case = Symbol(var) # e.g. :sto_soc
        symbol_event_case = Symbol(var * "2") # e.g. :sto_soc2

        plot!(plt_var, plot_window, value.(sp[symbol_base_case][plot_window]), label = "Base Case")
        plot!(plt_var, recovery_window, value.(sp[symbol_event_case]), label = "Event")

        push!(plots, plt_var)
    end

    plot(plots...; layout=(length(plots),1), size = (800, 150*length(plots)))
end

function sankey_results(sp, pv, w, el_d, timesteps)
    total_pv = value.(sp[1, :u_pv])*sum(pv[timesteps])
    total_wind = value.(sp[1, :u_wind])*sum(w[timesteps])
    total_demand = sum(el_d[timesteps])
    pv_cur = sum(value.(sp[1,:pv_cur])[timesteps])
    wind_cur = sum(value.(sp[1,:wind_cur])[timesteps])
    st_in = sum(value.(sp[1,:sto_to_bus])[timesteps])
    st_out = sum(value.(sp[1,:sto_from_bus])[timesteps])
    grid_in = sum(value.(sp[1,:gci])[timesteps])
    grid_out = sum(value.(sp[1,:gco])[timesteps])
    energy2heat = sum(value.(sp[1,:flow_energy2heat])[timesteps])
    losses_charge = 1/sp.stages[1].parameters[:sto_ef_ch] - 1.
    losses_discharge = 1. - sp.stages[1].parameters[:sto_ef_dis]
    storage_losses = losses_charge*sum(value.(sp[1,:sto_from_bus])[timesteps]) + losses_discharge*sum(value.(sp[1,:sto_to_bus])[timesteps])
    labels = ["PV", "Wind", "Storage (input)", "Storage (output)", "Demand", "Grid input", "Grid output", "Bus", "Losses", "Energy to heat", "Curtailment"]
    src = [1,1,2,2,4,6,8,8,8,8,8]
    trg = [8,11,8,11,8,8,3,5,7,9,10]
    weights = [total_pv-pv_cur, pv_cur, total_wind-wind_cur, wind_cur, st_in, grid_in, st_out, total_demand, grid_out, storage_losses, energy2heat]
    total_in = total_pv+total_wind+st_in+grid_in-pv_cur-wind_cur
    total_out = total_demand+st_out+grid_out+energy2heat
    println("relative mismatch = $((total_in-total_out)/total_in)")
    #println(storage_losses/(st_in+st_out))
    sankey(src, trg, weights, node_labels = labels)
end

function plot_scenario_distribution(scenarios; by_sign = false)
    if !by_sign
        scen_matrix = zeros((24,365))
        for s in scenarios
            t_xi = s.data.t_xi
            scen_matrix[mod(t_xi,24)+1,t_xi÷24+1] = 1
        end
        heatmap(scen_matrix, leg=false,
        xlabel = "Day of the year",
        ylabel = "Hour of the day", title = "Ancillary requests sample")
    else
        scen_matrix_pos = zeros((24,365))
        scen_matrix_neg = zeros((24,365))
    
        for s in scenarios
            t_xi = s.data.t_xi
            s_xi = s.data.s_xi
            F_xi = s.data.F_xi
            if s_xi == 1
                scen_matrix_pos[mod(t_xi,24)+1,t_xi÷24+1] = F_xi
            else
                scen_matrix_neg[mod(t_xi,24)+1,t_xi÷24+1] = F_xi
            end
        end
        p = plot(layout = (1,2))
        plot!(p[1],scen_matrix_pos, seriestype = :heatmap, ratio = 365/24, framestyle = :none)
        plot!(p[2],scen_matrix_neg, seriestype = :heatmap, ratio = 365/24, framestyle = :none)
    end 
end

using VegaLite
function plot_flex_sources(sp_data, pv, wind; timesteps = 1:24*365)
    pos_cur = sp_data["op"]["pv_cur"][timesteps[begin:end-1]] .+ sp_data["op"]["wind_cur"][timesteps[begin:end-1]]
    neg_cur = pos_cur .- pv[timesteps[begin:end-1]].*sp_data["inv"]["u_pv"] .- wind[timesteps[begin:end-1]].*sp_data["inv"]["u_wind"]
    pos_sto = sp_data["op"]["sto_soc"][timesteps[begin+1:end]]
    neg_sto = pos_sto .- sp_data["inv"]["u_storage"]
    COP = COP = sp_data["params"]["COP"]
    pos_heat_sto = sp_data["op"]["heat_sto_soc"][timesteps[begin+1:end]]/COP
    neg_heat_sto = pos_heat_sto .- sp_data["inv"]["u_heat_storage"]/COP
    #=plt = plot(layout = (1,2))
    plot!(plt[1], pos_cur, fill = (0, 0.5))
    plot!(plt[1], pos_sto, fill = (0, 0.5))
    plot!(plt[1], pos_heat_sto, fill = (0, 0.5))
    plot!(plt[2], neg_cur, fill = (0, 0.5))
    plot!(plt[2], neg_sto, fill = (0, 0.5))
    plot!(plt[2], neg_heat_sto, fill = (0, 0.5))=#
    p1 = areaplot(timesteps[begin:end-1], [pos_cur, pos_sto, pos_heat_sto], fillalpha = [0.5 0.5 0.5])
    p2 = areaplot(timesteps[begin:end-1], [-neg_cur, -neg_sto, -neg_heat_sto], fillalpha = [0.5 0.5 0.5])
    plot(p1, p2, layout = (1,2))
end

function plot_flex(sp_data, pv, wind; timesteps = 1:24*365)
    pos_cur = sp_data["op"]["pv_cur"][timesteps[begin:end-1]] .+ sp_data["op"]["wind_cur"][timesteps[begin:end-1]]
    neg_cur = pos_cur .- pv[timesteps[begin:end-1]].*sp_data["inv"]["u_pv"] .- wind[timesteps[begin:end-1]].*sp_data["inv"]["u_wind"]
    pos_sto = sp_data["op"]["sto_soc"][timesteps[begin+1:end]]
    neg_sto = pos_sto .- sp_data["inv"]["u_storage"]
    COP = COP = sp_data["params"]["COP"]
    pos_heat_sto = sp_data["op"]["heat_sto_soc"][timesteps[begin+1:end]]/COP
    neg_heat_sto = pos_heat_sto .- sp_data["inv"]["u_heat_storage"]/COP
    #=plt = plot(layout = (1,2))
    plot!(plt[1], pos_cur, fill = (0, 0.5))
    plot!(plt[1], pos_sto, fill = (0, 0.5))
    plot!(plt[1], pos_heat_sto, fill = (0, 0.5))
    plot!(plt[2], neg_cur, fill = (0, 0.5))
    plot!(plt[2], neg_sto, fill = (0, 0.5))
    plot!(plt[2], neg_heat_sto, fill = (0, 0.5))=#

    df_pos = DataFrame(time = timesteps[begin:end-1], cur = pos_cur, sto = pos_sto, heat_sto = pos_heat_sto)
    df_neg = DataFrame(time = timesteps[begin:end-1], cur = neg_cur, sto = neg_sto, heat_sto = neg_heat_sto)
    stack(df_pos,[:cur, :sto, :heat_sto]) |> @vlplot(:area, x = :time, y = {value, stack = :zero}, color="variable:n")
    return df_pos, df_neg
end
