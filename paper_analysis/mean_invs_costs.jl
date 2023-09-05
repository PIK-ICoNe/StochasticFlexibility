#=
We focus on three particular pairs of F and scen_freq:
(5000., 48), (25000., 48), (25000., 240)

In this file we gather data from the optimization and 
pack it into a single CSV file for convenience.
=#
## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()
using DataFrames
using CSV
using JSON
using BSON

#-
include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

#-
# Fist we get C^B
run_id = "new_parameters"
base_model = BSON.load(joinpath(basepath, "results", run_id, "baseline/baseline_0.0.bson"))
CB = base_model[:cost]
CIB = get_total_investment(base_model)
COB = get_operation_cost(base_model)
inv_base = base_model[:inv]
P = base_model[:params]
base_model = nothing;
#-
opt_params = [(500., 96, 20)]
timesteps = 1:24*365
savepath = joinpath(basepath, "results", run_id)
#-
df_costs = DataFrame(CSV.File(joinpath(savepath, "val_costs.csv")))
df_inv = DataFrame(CSV.File(joinpath(savepath, "inv.csv")))

#-
# Find mean investments and costs
cost_csv = joinpath(savepath, "mean_costs.csv")
invs_csv = joinpath(savepath, "mean_invs.csv")
if !isfile(cost_csv) || !isfile(invs_csv)
    df_inv_mean = DataFrame()
    df_cost_err = DataFrame()
    cost_vars = [:CF, :CR]
    for p in eachindex(opt_params)
        for m in 1:3
            sel = subset(df_costs, :m => n-> n.==m)
            mean_CF = mean(sel[!,:cost])-CB
            mean_CR = mean(sel[!,:CR])
            err_CF = std(sel[!,:cost].-CB)
            err_CR = std(sel[!, :CR])  
            # Total amount of flexibility requested in a year: F*events_per_year
            total_flex = (length(timesteps)-24) / (opt_params[p][2] + 1)*opt_params[p][1]*opt_params[p][3]*0.4 
            norm_costs = (sel[!, :cost] .- CB)/total_flex
            norm_mean_CF = mean(norm_costs)
            norm_err_CF = std(norm_costs)
            norm_mean_CR = mean_CR/total_flex
            norm_err_CR = std(sel[!, :CR])/total_flex
            append!(df_cost_err, DataFrame(mode = m, 
                mean_CF = mean_CF, mean_CR = mean_CR,
                err_CF = err_CF, err_CR = err_CR,
                norm_mean_CF = norm_mean_CF, norm_mean_CR = norm_mean_CR,
                norm_err_CF = norm_err_CF, norm_err_CR = norm_err_CR,))

            for (var_index, inv_var) in enumerate(inv_vars)
                append!(df_inv_mean, DataFrame(mode = m, var = var_index,
                    mean_value = mean(subset(df_inv, 
                    :m => mode -> mode .== m)[!,inv_vars[var_index]])*P[Symbol(replace(string(inv_var), "u_" => "c_"))]))
            end
        end
    end
    CSV.write(cost_csv, df_cost_err)
    CSV.write(invs_csv, df_inv_mean)
end
#-