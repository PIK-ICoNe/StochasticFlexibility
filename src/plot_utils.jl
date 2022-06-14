using Plots

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
    plot!(plt_heat, plot_span, value.(sp[1, :heat_sto_in])[plot_span]-value.(sp[1, :heat_sto_out])[plot_span], label = "heat storage use")
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
