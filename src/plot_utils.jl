using Plots

function plot_results(sp, pv, w, d; hd = 0, s=1, stage_1=[:gci, :gco], stage_2=[:gci2, :gco2], inv_dict = 0)
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
    plt_sto = plot()
    plt_invest = plot()
    plt = plot()
    plot!(plt_invest, pv .* u_pv, label="pv")
    plot!(plt_invest, w .* u_wind, label="wind")
    plot!(plt_invest, d, label="demand")
    t_xi = scenarios(sp)[s].data.t_xi
    recovery_time = sp.stages[2].parameters[:recovery_time]

    stor_charge = value.(sp[1, :sto_soc])
    plot!(plt_sto, 1:length(pv), stor_charge, label="global storage charge")
    plot!(plt_sto, (t_xi):(t_xi+recovery_time-1), value.(sp[2, :sto_soc2],s), label=string("stochastic storage charge")*string(s), linestyle=:dash, linewidth=2)

    for var in stage_1
        plot!(plt, 1:length(pv), value.(sp[1, var]), label=string(var))
    end

    for var in stage_2
        plot!(plt, (t_xi+1):(t_xi+recovery_time), value.(sp[2, var], s), label=string(var)*string(s), linestyle=:dash, linewidth=2)
    end
    if hd != 0.
        plt_heat = plot()
        COP = sp.stages[2].parameters[:COP]

        plot!(plt_heat, hd, label = "heat demand")
        plot!(plt_heat, COP*value.(sp[1, :heatpumpflow]), label = "heatpump")
        plot!(plt_heat, value.(sp[1, :heat_sto_in])-value.(sp[1, :heat_sto_out]), label = "heat storage use")
        plot!(plt_heat, (t_xi+1):(t_xi+recovery_time), COP*value.(sp[2, :heatpumpflow2], s), label = "heatpump$s", linestyle=:dash, linewidth=2)
        plot!(plt_heat, (t_xi+1):(t_xi+recovery_time), value.(sp[2, :heat_sto_in2], s)-value.(sp[2, :heat_sto_out2], s), label = "heat storage use $s", linestyle=:dash, linewidth=2)


        return plot(plt_invest, plt, plt_sto, plt_heat, layout = (4,1))
    else
        return plot(plt_invest, plt, plt_sto, layout = (3,1))
    end

end

function plot_recovery_window_deviation(sp; s = 1)
    plt_gc = plot()
    plt_sto = plot()
    t_xi = scenarios(sp)[s].data.t_xi
    recovery_time = sp.stages[2].parameters[:recovery_time]
    plot!(plt_gc, value.(sp[1, :heatpumpflow])[t_xi+1:t_xi+recovery_time], label = "global heatpumpflow")
    plot!(plt_gc, value.(sp[2, :heatpumpflow2], s), label = "stochastic heatpumpflow")
    plot!(plt_gc, -value.(sp[1, :gci])[t_xi+1:t_xi+recovery_time]+value.(sp[2, :gci2], s), label = "global grid connection use")
    plot!(plt_gc, value.(sp[2, :gco2], s) - value.(sp[1, :gco])[t_xi+1:t_xi+recovery_time], xlabel = "time after the event, h", label = "stochastic grid connection use")
    #plot!(plt_gc, -value.(sp[1, :gco])[t_xi+1:t_xi+recovery_time]+value.(sp[2, :gco2], s), xlabel = "time after the event, h", label = "gco2-gco")
    
    stor_charge = value.(sp[1, :sto_soc])[t_xi+1:t_xi+recovery_time]
    stor_charge2 = value.(sp[2, :sto_soc2], s)
    
    plot!(plt_sto, stor_charge2-stor_charge, label = "soc2-soc")


    plot(plt_gc, plt_sto, layout = (2, 1))
end
