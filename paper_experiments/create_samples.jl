#=
Here we create 20 samples for the full year
for each of the three chosen parameter sets.=#
basepath = realpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(basepath)
include(joinpath(basepath, "src", "sp_model.jl"))
using JSON
#-
# Pkg.instantiate()
params = [(5000., 48, 25), (7500., 144,70), (25000., 240, 120)]

savepath = joinpath(basepath, "samples")
if !isdir(savepath)
    mkdir(savepath)
end

t_max_offset = 24
t_max = 24*365 - t_max_offset
#n_samples = 30
recovery_time = 12#pars[:recovery_time]
#-
for i in eachindex(params)
    F, scen_freq, n_samples = params[i]
    savepath_i = joinpath(savepath, "scen_freq$(scen_freq)_n_$(n_samples)")
    if !isdir(savepath_i)
        mkdir(savepath_i)
    end
    for j in 1:20
        @assert scen_freq > recovery_time
        delta_t = scen_freq - recovery_time
        event_per_scen = t_max / (delta_t + recovery_time + 1)
        n = round(Int, n_samples * event_per_scen)
        scens = poisson_events_with_offset(n, delta_t, recovery_time, F, t_max, F_min = F*0.6)
        scens_dict = [Dict((:probability => s.probability, :t_xi => s.data[:t_xi], 
        :F_xi => s.data[:F_xi], :s_xi => s.data[:s_xi])) for s in scens]
        open(joinpath(savepath_i, "sample$j.json"), "w") do f
            JSON.print(f,scens_dict)
        end
    end
end