using Plots
using SankeyPlots

function plot_results(sp, pv, w, el_d; plot_span = 1:length(pv), s=1, el_balance_vars=[:gci, :gco], storage_vars=[], inv_dict = 0)
    if inv_dict != 0
        u_pv = inv_dict[:u_pv]
        u_wind = inv_dict[:u_wind]
        u_storage = inv_dict[:u_storage]
        u_heatpump = inv_dict[:u_heatpump]
        u_heat_storage = inv_dict[:u_heat_storage]
    else
        u_pv = value.(sp[1, :u_pv])
        u_wind = value.(sp[1, :u_wind])
        u_storage = value.(sp[1, :u_storage])
        u_heatpump = value.(sp[1, :u_heatpump])
        u_heat_storage = value.(sp[1, :u_heat_storage])
    end
    plt_sto = plot(; legend = :outertopright)
    plt_invest = plot(; legend = :outertopright)
    plt = plot(; legend = :outertopright)
    plot!(plt_invest, plot_span, pv[plot_span] .* u_pv, label="PV")
    plot!(plt_invest, plot_span, w[plot_span] .* u_wind, label="Wind")
    plot!(plt_invest, plot_span, el_d[plot_span], label="Electrical demand")
    t_xi = scenarios(sp)[s].data.t_xi
    recovery_time = sp.stages[2].parameters[:recovery_time]

    stor_charge = value.(sp[1, :sto_soc])
    plot!(plt_sto, plot_span, stor_charge[plot_span], label="global storage charge")
    plot!(plt_sto, (t_xi):(t_xi+recovery_time-1), value.(sp[2, :sto_soc2],s), label=string("stochastic storage charge")*string(s), linestyle=:dash, linewidth=2)

    for var in el_balance_vars
        plot!(plt, plot_span, value.(sp[1, var])[plot_span], label=string(var))
        var2 = Symbol(string(var)*"2")
        plot!(plt, (t_xi+1):(t_xi+recovery_time), value.(sp[2, var2], s), label=string(var)*" 2nd stage, s = "*string(s), linestyle=:dash, linewidth=2)

    end
    return plot(plt_invest, plt, plt_sto, layout = (3,1))
end

function plot_heat_layer(sp, heatdemand; plot_span = 1:length(heatdemand), s = 1, inv_dict = 0)
    plt_heat = plot(; legend = :outertopright)
    COP = sp.stages[2].parameters[:COP]

    plot!(plt_heat, plot_span, heatdemand[plot_span], label = "heat demand")
    plot!(plt_heat, plot_span, COP*value.(sp[1, :flow_energy2heat])[plot_span], label = "heatpump")
    plot!(plt_heat, plot_span, value.(sp[1, :heat_sto_soc])[plot_span], label = "heat storage SOC")
    plot!(plt_heat, plot_span, value.(sp[1, :heat_sto_to_bus])[plot_span]-value.(sp[1, :heat_sto_from_bus])[plot_span], label = "heat storage use")
    # plot!(plt_heat, (t_xi+1):(t_xi+recovery_time), COP*value.(sp[2, :flow_energy2heat2], s), label = "heatpump$s", linestyle=:dash, linewidth=2)
    # plot!(plt_heat, (t_xi+1):(t_xi+recovery_time), value.(sp[2, :heat_sto_in2], s)-value.(sp[2, :heat_sto_out2], s), label = "heat storage use $s", linestyle=:dash, linewidth=2)

end

function plot_recovery_window_deviation(sp; s = 1)
    plt_gc = plot()
    plt_sto = plot()
    t_xi = scenarios(sp)[s].data.t_xi
    recovery_time = sp.stages[2].parameters[:recovery_time]
    plot!(plt_gc, value.(sp[1, :flow_energy2heat])[t_xi+1:t_xi+recovery_time], label = "global flow_energy2heat")
    plot!(plt_gc, value.(sp[2, :flow_energy2heat2], s), label = "stochastic flow_energy2heat")
    plot!(plt_gc, -value.(sp[1, :gci])[t_xi+1:t_xi+recovery_time]+value.(sp[2, :gci2], s), label = "global grid connection use")
    plot!(plt_gc, value.(sp[2, :gco2], s) - value.(sp[1, :gco])[t_xi+1:t_xi+recovery_time], xlabel = "time after the event, h", label = "stochastic grid connection use")
    #plot!(plt_gc, -value.(sp[1, :gco])[t_xi+1:t_xi+recovery_time]+value.(sp[2, :gco2], s), xlabel = "time after the event, h", label = "gco2-gco")
    
    stor_charge = value.(sp[1, :sto_soc])[t_xi+1:t_xi+recovery_time]
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
    plot_window = t_xi+1+window_start:t_xi+length(sp[:gco2])+window_end

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

function sankey_results(sp, pv, w, el_d, timesteps)
    total_pv = value.(sp[1, :u_pv])*sum(pv[timesteps])
    total_wind = value.(sp[1, :u_wind])*sum(w[timesteps])
    total_demand = sum(el_d[timesteps])
    st_in = sum(value.(sp[1,:sto_to_bus])[timesteps])
    st_out = sum(value.(sp[1,:sto_from_bus])[timesteps])
    grid_in = sum(value.(sp[1,:gci])[timesteps])
    grid_out = sum(value.(sp[1,:gco])[timesteps])
    losses_charge = 1/sp.stages[1].parameters[:sto_ef_ch] - 1.
    losses_discharge = 1. - sp.stages[1].parameters[:sto_ef_dis]
    storage_losses = losses_charge*sum(value.(sp[1,:sto_from_bus])[timesteps]) + losses_discharge*sum(value.(sp[1,:sto_to_bus])[timesteps])
    labels = ["PV", "Wind", "Storage (input)", "Storage (output)", "Demand", "Grid input", "Grid output", "Bus", "Losses"]
    src = [1,2,4,6,8,8,8,8]
    trg = [8,8,8,8,3,5,7,9]
    weights = [total_pv, total_wind, st_in, grid_in, st_out, total_demand, grid_out, storage_losses]
    total_in = total_pv+total_wind+st_in+grid_in
    total_out = total_demand+st_out+grid_out
    println("relative mismatch = $((total_in-total_out)/total_in)")
    #println(storage_losses/(st_in+st_out))
    sankey(src, trg, weights, node_labels = labels)
end