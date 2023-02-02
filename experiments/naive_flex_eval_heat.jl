## Everything runs in the Project environment on the basepath

#=
In this experiment we compare flexibility potential of models with and without foreknowldge.

First we load the data for the background model. We fix the 
=#
basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()

using Dates
using Plots

#= Include the file containing all the dependecies and optimization functions.
=#
#-
include(joinpath(basepath, "experiments", "stochastic_optimization_setup.jl"));

#-
warm_up()
#-
bkg_cost = []
F = 5000.:5000.:25000.
append!(bkg_cost, CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg.csv"), DataFrame)[!,1][1])
for f in F
    append!(bkg_cost, CSV.read(joinpath(basepath, "results/naive_flex/bkg", "cost_bkg_$f.csv"), DataFrame)[!,1][1])
end


#n_samples = 10:10:70 #[collect(2:2:10); collect(15:10:25)] #[10]
#scen_freq = 72:48:24*30#24:24:30*24 #[24]

#n_samples = [collect(15:5:35); collect(40:10:80); 100]
#scen_freq = 48:24:120
#-
n_samples = [collect(15:5:40); collect(50:10:100); collect(125:25:150)]
scen_freq = 24*7

flex_cost = zeros(length(n_samples),length(F));
runtime = zeros(length(n_samples),length(F));
investments = zeros(length(n_samples),length(F));
#-
u_wind_data = [];
wind_invs = zeros(length(n_samples),length(F));
storage_invs = zeros(length(n_samples),length(F));

wind_full = [];
storage_full = [];

run_id = "naive_flex/2023-01-18"
for i in eachindex(n_samples)
     for j in eachindex(F)
        event_per_scen = 365*24 / (96 + 1)
        costs = load_costs(run_id, n_samples[i], F[j])
        mean_cost = mean(costs)
        invs, mean_invs = load_invs(run_id, n_samples[i], F[j])
        wind_invs[i,j] = mean_invs[:u_wind]
        storage_invs[i,j] = mean_invs[:u_storage]#
        #=if i == 7
            push!(wind_full,(scen_freq[j], invs[:u_wind]))
            push!(storage_full,(scen_freq[j], invs[:u_storage]))
        end=#
            #plot!(plt_inv, ([scen_freq[i] for n in 1:10],invs[:u_wind]), seriestype = :scatter )
        append!(u_wind_data, [invs[:u_wind]])
        runtimes = load_runtime(run_id, n_samples[i], F[j])
        mean_runtime = mean(runtimes)
        # flex_cost - cost of optimized system with mean investment and new sample TODO
        flex_cost[i,j] = (mean_cost - bkg_cost[1])/round(event_per_scen)#*n_samples[i]
        runtime[i,j] = mean_runtime
    end
end
#-
#=p = plot(heatmap(z=flex_cost, y=n_samples,x=[string(s) for s in scen_freq]),
    Layout(yaxis_title="Scaled sample size", xaxis_title="Average offset between flexibility events, hours",
    xaxis_nticks = 4))
savefig(p, joinpath(basepath,"flex_cost.png"))=#
#
heatmap(F, n_samples, flex_cost,
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", 
    title = "Cost of flexibility")
#-
heatmap(F, n_samples[2:end], runtime[2:end, :], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", 
    title = "Runtime, s")

heatmap(F, n_samples, wind_invs, 
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", 
    title = "Wind, kWp")

heatmap(F, n_samples[2:end],storage_invs[2:end, :], 
    xlabel = "Guaranteed flexibility",
    ylabel = "Scaled sample size", title = "Storage, kWh")

#-
#=for (i, n_samples) in enumerate(n_samples)
    event_per_scen = 365*24 / (scen_freq + 1)
    flex_cost[i] = (objective_value(sps[i]) .- cost_bkg)/event_per_scen
end=#

#-
# Analyze operational schedule for fixed guaranteed flexibility:
#-
ns = 150
F = 5000.

op, = load_ops("naive_flex/2023-01-18", ns, F)
invs, = load_invs(run_id, ns, F)
scen_data = CSV.read("/home/ekatzolo/Nextcloud/code/StochasticFlexibility/results/naive_flex/2023-01-18/scen$(ns)_$F.csv", DataFrame)
scen = Dict(pairs(eachcol(scen_data)))
t_xi = scen[:t_xi] |> sort!

#-

window = 1:length(t_xi)
t_window = t_xi[window[begin] .<= t_xi .<= window[end]]

# plot(collect(window), op[:heat_sto_to_bus][window] .- op[:heat_sto_from_bus][window])
plot(collect(window), op[:sto_soc][window])
scatter!(t_window, zeros(length(t_window)))

#-

include(joinpath(basepath, "src", "plot_utils.jl"));

#plot_scenario_debug(sp,100)

#-
scen_formatted = [@scenario t_xi = scen[:t_xi][i] F_xi = scen[:F_xi][i] s_xi = scen[:s_xi][i] for i in eachindex(scen[:t_xi])]
es_fix = define_energy_system(pv, wind, demand, heatdemand; p = pars)
sp_fix = instantiate(es_fix, scen_formatted, optimizer = Clp.Optimizer)
set_silent(sp_fix)
fix_investment!(sp_fix, invs)
fix_operation!(sp_fix, op, length(timesteps))
optimize!(sp_fix)
#-
savefile_lock = ReentrantLock()

sp, = optimize_sp(pv, wind, demand, heatdemand, pars, 20, 24*7, savefile_lock, F_max = 5000., F_min = 0., F_pos = 25000, F_neg = -25000.)
scens = scenarios(sp)
plot_outcome(sp, scens[1].data[:t_xi], scens[1].data[:s_xi], scens[1].data[:F_xi])
#-
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat=true);

plot_flex_sources(op, invs, pv, wind, pars[:COP], timesteps = 100:110)

