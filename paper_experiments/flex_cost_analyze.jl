#=
We focus on three particular pairs of F and scen_freq:
(5000., 48), (7500., 144), (25000., 240)
=#

## Everything runs in the Project environment on the basepath

basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
# Pkg.instantiate()
using DataFrames

using JSON
#-
include(joinpath(basepath, "src", "sp_model.jl"))
include(joinpath(basepath, "src", "plot_utils.jl"))
include(joinpath(basepath, "src", "data_load.jl"));

#-
n_samples = 20
params = [(5000., 48), (7500., 144), (25000., 240)]
n_runs = 20
timesteps = 1:24*365
pv, wind, demand, heatdemand, pars = load_max_boegl(timesteps, heat = true);
pars[:inv_budget] = 10^10;
savepath = joinpath(basepath, "results/flex_cost")
#-
df = DataFrame()
#-
for i in eachindex(params)
    F, sf = params[i]
    for n in 6:15
        for m in ("fixed_fg", "fixed_fg_inv", "unfixed_guar_flex")
            p = joinpath(savepath, "run_$(F)_$(sf)_$(n)_$(m).json")
            if isfile(p)
                opt = JSON.parsefile(p)
                total_cost = opt["cost"]
                CI = get_total_investment(opt)
                CO = get_operation_cost(opt)
                CR = total_cost - CI - CO
                append!(df, DataFrame(p = i, total_cost = total_cost,
                    inv_cost = CI, op_cost = CO, flex_op_cost = CR, 
                    #[var = opt["inv"][var] for var in string.(inv_vars)],
                    sample_num = n, mode = m))
                    opt = nothing;
            else
                println("File $p does not exist")
            end
        end
    end
end
#-
df
