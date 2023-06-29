## Everything runs in the Project environment on the basepath
#This is the manuscipt for the Thesis "Analysing the Benefits of Stochastic Programming for Energy Cell's Ancillary Services Availability"

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
Pkg.instantiate()

using DataFrames
using CSV
using Clp
using JSON
using BSON
using Statistics;
using PlotlyJS

using Random
Random.seed!(1);

#Include the other necessary files: the stochastic programming model, plot functions and analysis tools for more detailed evaluation    

include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "evaluation_utils.jl"));

#The time steps of the system are determined here
t_max_offset = 24 #not using the last day of the year, it is not allowed to have a request in the last day  
offset = 0
timesteps = 1:(24*365)

#loading a timeseries for photovoltaic (pv), wind potial and the demnad for electricity and heat"
pv_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_pv_pik.csv"), DataFrame, header = false)
wind_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_ninja_wind_pik.csv"), DataFrame, header = false)
demand_data =  CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_apartment_block(KFW40)_300km2_6600MWhperyear.csv"), DataFrame, header = false)
heatdemand_data = CSV.read(joinpath(basepath, "timeseries/validation_usecase", "csv_heatdemand_apartment_block(KFW40)_300km2_12300MWhperyear.csv"), DataFrame, header = false)

#defining the offset of each timesiries
pv = pv_data[timesteps .+ offset, 1]
wind = wind_data[timesteps .+ offset, 1]
demand = demand_data[timesteps .+ offset, 1]
heatdemand = heatdemand_data[timesteps .+ offset, 1]
#heatdemand = zeros(length(timesteps)) #if heatdemand is not needed 

t_max = minimum((length(pv), length(wind), length(demand), length(heatdemand))) - t_max_offset #using the longest possible timeintervall 

data = nothing; # Free the memory
heatdemand_data = nothing;
pv_data = nothing;
wind_data = nothing;
demand_data = nothing; # Free the memory

total_demand = heatdemand + demand

# Generate x-axis labels for every second month
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
x_labels = [months[i] for i in 1:1:12]

# Create a vector of indices for x-axis ticks
time = collect(1:8760)
x_ticks = collect(1:730:8760)

y_max = 9000

# Plotting the data
pl = Plots.plot(time, total_demand, xlabel = "Months", ylabel = "Total Demand (kW)", yticks = 0:1000:9000,
     xticks = (x_ticks, x_labels), legend = false, linewidth = 0.5, ylim = (0, y_max), fmt = :png)
plot!(pl, time, heatdemand)     


#plotting the data (unitless for availability timeseries),  
plt = Plots.plot(timesteps, pv .* (mean(demand) / mean(pv)), label="PV (unitless)")
plot!(plt, timesteps, wind.* (mean(demand) / mean(wind)), label="Wind (unitless)")
plot!(plt, timesteps, heatdemand, label="Heat Demand")
plot!(plt, timesteps, demand, label="Electric Demand")
plt

#Here the default parameters of the model are replaced with realistic data for an analysis, more detailed information in the thesis. 
pars = copy(default_es_pars)

#Indication of the mean value for a better idea of the demand series 
average_hourly_demand = mean(demand)
average_hourly_heatdemand = mean(heatdemand)

#Introduction of recovery time for the second stage to prevent perfect foresight of a linear optimisation model 
recovery_time = 12

pars[:recovery_time] = recovery_time
pars[:c_storage] = 600.
pars[:c_pv] = 900.
pars[:c_wind] = 2500.
pars[:c_i] = 0.3 #0.3
pars[:c_o] = 0.03 #0.03
pars[:c_heat_storage] = 200.
pars[:heat_losses] = 0.00416
pars[:heat_eff] = 0.95
pars[:asset_lifetime] = 20.
pars[:c_heatpump] = 533. # cost for 1 kW heat output ->  1600€ for 1 kW electricity in 3.5 KW heat
pars[:max_pv] = 10^4.
pars[:max_wind] = 10^4.

#Here the path is determined for saving and analysing the results
savepath = joinpath(basepath, "results")
if !isdir(savepath)
    mkdir(savepath)
end


#The model itself is constructed by the function define_energy_system: First Energy System background (only first stage), second energy system without Regularization (two stage model)   
es_bkg = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = false, override_no_event_per_scen = true)
es_no_reg = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = false)

#=
This contains investment variables which we collectively call $I$, and the operational schedule $O^t$.
The total cost given a certain investment $I$ and schedule $O^t$ is denoted $C(I, O^t)$.
=#

#First, the Energy System Background is optimised, which does not take into account the need for anyillary services. A "pseudo sampler" is used so that there is no need for generated ancillary service.
sp_bkg = instantiate(es_bkg, no_flex_pseudo_sampler(), optimizer = Clp.Optimizer)

# Prevent investment in the heat components if there is no heat demand
if maximum(heatdemand) == 0.
    fix(decision_by_name(sp_bkg, 1, "u_heatpump"), 0.)
    fix(decision_by_name(sp_bkg, 1, "u_heat_storage"), 0.)
end

set_silent(sp_bkg) #no direct output from the solver 

optimize!(sp_bkg) #the command for optimisation 

bkg_decision = optimal_decision(sp_bkg) #Results

objective_value(sp_bkg) #Objective function value

#Storage the variables in an manageable file format (JSON)
all_data = get_all_data(sp_bkg, rec = false, scen = false)
all_data = merge(all_data, Dict((:cost => objective_value(sp_bkg))))
open(joinpath(savepath, "bkg.json"), "w") do f
    JSON.print(f,all_data)
end


#Plot the results in an hourly dispatch view 
plot_results(all_data, pv, wind, demand)

#Plot the results in an hourly dispatch view 
plot_heat_layer(all_data, heatdemand)


all_data = nothing #clean the memory 


#The following code for calculating the two-stage stochastic optimisation was solved on the cluster of the Potsdam Institute for Climate Impact Research and is presented here for understanding the logical structure of the work 
t_max = length(pv) - 24 #Preventing the occurrence of an Ancillary Service Request in the last 24 hours 
F_max = 5000. #Here the maximum and minimum possible value of the random parameter Ancillary Service Request is defined
F_min = 3000. #300kW to 1000 kW and 1000kW to 3000kW were also analysed in the thesis
delta_t = 3*24 # Flex event every delta_t + recovery time + 1 = 85h (3 and half day)
#pars[:event_per_scen] = t_max / (delta_t + recovery_time + 1); #102 events (ancillary service request per scen) that is used for a further development, but not part of the thesis
n_samples = [collect(15:5:35); collect(40:10:90); collect(100:25:125)] #how big is the sample size (the number of scenarios)


for ns in n_samples
    n = round(Int, ns * pars[:event_per_scen]) #pars[:event_per_scen]=1 (default)
    scens = a(n, delta_t, recovery_time, F_max, t_max, F_min = F_min)
    es_st = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = true) #defining the stochastic energy system with flexibility
    sp_st = instantiate(es_st, scens, optimizer = Clp.Optimizer)
    set_silent(sp_st)
    optimize!(sp_st)
    println("termination_status=", termination_status(sp_st)) #sucessfull or infeasible #all samples should be sucessfull because we are using the reularized model with an penalty 

    opt_params = Dict((:F_min => F_min, :F_max => F_max, :t_max_offset => t_max_offset, :n_samples => n_samples, :scen_freq => scen_freq))
    all_data = get_all_data(sp_st)
    all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp_st))))
    bson("run_$(n_samples)_$(scen_freq).bson", all_data)
    if ns == n_samples[end]
        ns = n_samples[end]
        es_st = define_energy_system(pv, wind, demand, heatdemand; p = pars, regularized = true) #defining the stochastic energy system with flexibility
        sp_st = instantiate(es_st, scens, optimizer = Clp.Optimizer)
        fix_investment!(sp_st, get_investments(sp_bkg))
        set_silent(sp_st)
        optimize!(sp_st)
        opt_params = Dict((:F_min => F_min, :F_max => F_max, :t_max_offset => t_max_offset, :n_samples => n_samples, :scen_freq => scen_freq))
        all_data = get_all_data(sp_st)
        all_data = merge(all_data, opt_params, Dict((:runtime => runtime, :cost => objective_value(sp_st))))
        bson("run_$(n_samples)_$(scen_freq)_fixed_inv.bson", all_data)
    end
end



#From here the analysis of the system begins
#Initialisation of the variables with the number of samples 
costs = zeros(length(n_samples))
u_pv = zeros(length(n_samples))
u_wind = zeros(length(n_samples))
u_storage = zeros(length(n_samples))
u_heat_storage = zeros(length(n_samples))
u_heatpump = zeros(length(n_samples))
scen_freq = delta_t

# Iterate over n_samples, load optimization data, and extract relevant values for analysis
for (i,ns) in enumerate(n_samples) 
    opt = BSON.load((joinpath(basepath, "timeseries/validation_usecase/Paul 2/3000_5000","run_$(n_samples[i])_$(scen_freq).bson")))
    costs[i]= opt[:cost]
    u_pv[i]= opt[:inv][:u_pv]
    u_wind[i] = opt[:inv][:u_wind]
    u_storage[i] = opt[:inv][:u_storage]
    u_heat_storage[i] = opt[:inv][:u_heat_storage]
    u_heatpump[i] = opt[:inv][:u_heatpump]
end

#Re-use of parameters due to free memory
c_storage = 600.
c_pv = 900.
c_wind = 2500.
c_i = 0.3
c_o = 0.03
c_heat_storage = 400.
c_heatpump = 533.

#Third research question
plot_sample_investments(u_pv, u_wind, u_storage, u_heat_storage, u_heatpump, n_samples, c_pv, c_wind, c_storage, c_heat_storage, c_heatpump)
plot_sample_investments_norm(u_pv, u_wind, u_storage, u_heat_storage, u_heatpump, n_samples, c_pv, c_wind, c_storage, c_heat_storage, c_heatpump)


#save in results(savepath) the BSON file and make them readable 
for (i,ns) in enumerate(n_samples) 
    open(joinpath(savepath, "run_$(n_samples[i])_$(scen_freq).bson"), "w") do f
        BSON.print(f,BSON.load((joinpath(basepath, "timeseries/validation_usecase/Paul 2/3000_5000","run_$(n_samples[i])_$(scen_freq).bson"))))
    end
end

#the same for fixed investments, which is not covered by the loop 
open(joinpath(savepath, "run_125_72_fixed_inv.bson"), "w") do f
    BSON.print(f,BSON.load((joinpath(basepath, "timeseries/validation_usecase/Paul 2/3000_5000","run_125_72_fixed_inv.bson"))))
end


#extract background cost 
bkg_cost = JSON.parsefile((joinpath(savepath, "bkg.json")))["cost"]


#First research question
#Comparing the fixed investment with ancillary services request and with the optimised two stage model, where the investment/ bkg is taken into account 
cost_fixed_inv_125_72 = BSON.load((joinpath(basepath, "timeseries/validation_usecase/Paul 2/3000_5000", "run_125_72_fixed_inv.bson")))[:cost]
cost_125_72 = BSON.load((joinpath(basepath, "timeseries/validation_usecase/Paul 2/3000_5000", "run_125_72.bson")))[:cost]
s = ["fixed investment optimization", "two-stage optimization"]
y = [cost_fixed_inv_125_72, cost_125_72]
comp_cost = Plots.bar(s, y, label = "", ylabel = "total system cost") #,yticks = 0:length(y):maximum(y)
for i in eachindex(s, y)
    Plots.annotate!(comp_cost, i, y[i] + 0.5, text(string(y[i]), :black, 8, :left))
end
display(comp_cost)

#second research question
plot_sample_cost_of_flex(n_samples, bkg_cost, costs)


#Fourth research question        
pl = Plots.plot()
n_samples = [collect(15:5:35); collect(40:10:90); collect(100:25:125)]
F_positive = zeros(length(n_samples))
F_negativity = zeros(length(n_samples))

# color codes giving a gradient from green to pink 
plotscolors = ["#00FF00", "#19F30D", "#33E61A", "#4CD927", "#66CE34", "#7FB741", "#99AB4E", "#B2995B", "#CC8668", "#E57375", "#FF617F", "#FF4D8F", "#FF3A9F", "#FF26AF"]

for (i,ns) in enumerate(n_samples)
    analysed_model = BSON.load((joinpath(basepath, "timeseries/validation_usecase/Paul 2/1000_3000","run_$(n_samples[i])_$(scen_freq).bson")))
    F_pos, F_neg, = naive_flex_potential(analysed_model, pv, wind, timesteps)
    Plots.plot!(pl, F_pos, label = "F_pos/neg_$(n_samples[i])_$(scen_freq).bson", linecolor = plotscolors[i], xlabel = "Timesteps in h", ylabel = "Flexibility in kW")
    Plots.plot!(pl, F_neg, label = "", linecolor = plotscolors[i])
    F_positive[i] = sum(F_pos)
    F_negativity[i] = sum(F_neg) 
end

#the size of the window is decided here 
xlims!(pl, 2000, 2500) #1000, 3000 or 2000, 2500
#display(pl)

#with bkg_model 
bkg_model = JSON.parsefile((joinpath(savepath, "bkg.json")))
F_pos, F_neg, = naive_flex_potential(bkg_model, pv, wind, timesteps)
Plots.plot!(pl, F_pos, label = "F_pos/neg_bkg", linecolor = "blue", linealpha = 0.3)
Plots.plot!(pl, F_neg, label = "", linecolor = "blue", linealpha = 0.3)
F_positive_bkg = sum(F_pos)
F_negative_bkg = sum(F_neg)

Plots.plot!(pl, legend = false, linewidth = 0.5, ylim = (-18000,10000), yticks = -18000:2000:10000)  #show the legend 
Plots.yaxis!(pl, formatter = :plain)


#Not used for thesis
#function correlation_coefficient(x, y)
#    n = length(x)
#    mean_x = sum(x) / n
#    mean_y = sum(y) / n
#
#    numerator = sum((x .- mean_x) .* (y .- mean_y))
#    denominator = sqrt(sum((x .- mean_x).^2) * sum((y .- mean_y).^2))
#
#    return numerator / denominator
#end

#Attention: this is only meaningful for a linear relationship and a cardinal sklaierte data, linear relationship may be doubted
correlation = correlation_coefficient(n_samples, F_positive)

#generate the plots for the fourth research question
Plots.plot(n_samples, F_positive, label="Ratio F_positive to sample size", xlabel = "sample size", ylabel = "F_positive") #maybe better without correlation #, title="Correlation: $correlation"
Plots.plot(n_samples, F_negativity, label="Ratio F_negativity to sample size", xlabel = "sample size", ylabel = "F_negativity")

#with bkg 
pushfirst!(F_positive, F_positive_bkg)
pushfirst!(F_negativity, F_negative_bkg)
pushfirst!(n_samples, 0)
Plots.plot(n_samples, F_positive, label="", xlabel = "sample size", ylabel = "F_positive") #maybe better without correlation #, title="Correlation: $correlation" #label="Ratio F_positive to sample size"
Plots.plot(n_samples, F_negativity, label="", xlabel = "sample size", ylabel = "F_negativity") #label="Ratio F_negativity to sample size"


####
#Additional research opportunity from Ekaterina 
testbson = BSON.load((joinpath(basepath, "timeseries/validation_usecase/Paul","run_15_72.bson")))

# Stack area chart for flex potential
flex_sources = [:pos_cur, :pos_sto, :pos_heat, :neg_cur, :neg_sto,:neg_heat]
plot_labels = Dict(("cur" => "Curtailment","sto" => "Electrical storage", "heat" => "Heat storage"))
include(joinpath(basepath, "src", "evaluation_utils.jl"))
p_colors = ["rgb(255, 0, 0)", "rgb(0, 255, 0)", "rgb(0, 0, 255)"]
#p_colors = ["#FF0000", "#0000FF", "#00FF00"]
plot_colors = Dict(("cur"=>p_colors[1], "sto"=>p_colors[2], "heat"=>p_colors[3]))
plot_window = 1:100

        test_model = testbson# filename here BSON.load(joinpath(savepath, "preval","run_$(F)_$(scen_freq)_2_$(m).bson"))
        F_pos, F_neg, F_pos_d, F_neg_d = naive_flex_potential(test_model, pv, wind, timesteps)
        flex_full = DataFrame(t=plot_window, pos_cur=F_pos_d[plot_window,1], pos_sto=F_pos_d[plot_window,2], pos_heat=F_pos_d[plot_window,3],
            neg_cur = F_neg[plot_window,1],neg_sto=F_neg_d[plot_window,2],neg_heat=F_neg_d[plot_window,3])
        pl = PlotlyJS.plot()
        for flex_var in flex_sources
            fl_sign, fl_src = split(String(flex_var),"_")
            @show plot_colors[fl_src]
            add_trace!(pl,PlotlyJS.scatter(
                x=plot_window, y=flex_full[plot_window,flex_var],
                stackgroup=fl_sign, mode="lines", hoverinfo="x+y",
                line=attr(width=0.5), name = plot_labels[fl_src], color=plot_colors[fl_src]
            ))
        end
        PlotlyJS.savefig(pl,joinpath(basepath, "results", "flex_sources2.png"))