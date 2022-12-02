using Plots

"""
Suppress output of evaluate_decision, and deal with solvers that error on infeasibility.
"""
function evaluate_decision_wrapper(p, decision, scenario)
    cost = 0.
    try redirect_stdout((() -> cost = evaluate_decision(p, decision, scenario)), devnull)
    catch e
        cost = Inf
    end
    return cost
end

"""
Find maximum available flexibility and its cost at all points in some time interval.
"""
function analyze_flexibility_potential(sp, timesteps, penalty_cost, no_flex_cost; decision = optimal_decision(sp))
    cost_no_flex = evaluate_decision_wrapper(sp, optimal_decision(sp), no_flex_pseudo_sampler()[1])
    L = length(timesteps)
    cost_pos_flex = zeros(L)
    cost_neg_flex = zeros(L)
    potential_pos_flex = zeros(L)
    potential_neg_flex = zeros(L)
    Threads.@threads for i in eachindex(timesteps)
        t = timesteps[i]
        potential_pos_flex[i], cost_pos_flex[i] = find_f_max(sp,t,1,decision, penalty_cost, no_flex_cost)
        potential_neg_flex[i], cost_neg_flex[i] = find_f_max(sp,t,-1,decision, penalty_cost, no_flex_cost)
    end
    return cost_pos_flex .- cost_no_flex, potential_pos_flex, cost_neg_flex .- cost_no_flex, potential_neg_flex
end

function plot_flexibility(timesteps, cost_pos_flex, potential_pos_flex, cost_neg_flex, potential_neg_flex)
    plt_cost = plot()
    plt_pot = plot()
    plot!(plt_cost, timesteps, cost_pos_flex ./ potential_pos_flex, label = "price of positive flexibility")
    plot!(plt_cost, timesteps, cost_neg_flex ./ abs.(potential_neg_flex), label = "price of negative flexibility")
    plot!(plt_pot, timesteps, potential_pos_flex, fillrange = 0, fillalpha = 0.35, label = "positive flexibility potential")
    plot!(plt_pot, timesteps, potential_neg_flex, fillrange = 0, fillalpha = 0.35, label = "negative flexibility potential")
    plot(plt_cost, plt_pot, layout = (2, 1))
end

"""
Find distribution of maximum available flexibilities.
"""
function flexibility_availability!(plt, flex_potential; plot_options...)
    p_sorted = sort(unique(flex_potential))
    fraction = [sum(abs.(flex_potential[:]).>=abs(f)) for f in p_sorted]./length(p_sorted)
    plot!(plt, p_sorted, fraction)
end

"""
Bisection search for maximum felxibility potential.
Parameters:
- sp - optimized stochastic program
- t - time step
- s - sign of flexibility scenario
- tol - result tolerance
- maxiters - maximum number of iterations
"""
function find_f_max(sp, t, s, od, penalty_cost, no_flex_cost; tol = 10., maxiter = 500)
    a = 0.
    b = 10000.
    cost_a = 0.
    cost_b = 0.
    i = 0
    scen = @scenario t_xi = t s_xi = s F_xi = a probability = 1.
    cost_a = evaluate_decision_wrapper(sp,od,scen)
    if cost_a == Inf
        return 0., cost_a
    else
        while b-a>tol && i<=maxiter && cost_a < no_flex_cost + penalty_cost# gi != 0, ...
            i+=1
            scen = @scenario t_xi = t s_xi = s F_xi = b probability = 1.
            cost_b = evaluate_decision_wrapper(sp,od,scen)
            if cost_b == Inf || cost_b > no_flex_cost + penalty_cost
                b = (a+b)/2
            else
                a = b
                cost_a = cost_b
                b *= 2
            end
        end
        return s*a, cost_a
    end
end

function evaluate_cost(sp, decision, timesteps, F)
    cost_pos = zeros(length(timesteps))
    cost_neg = zeros(length(timesteps))
    Threads.@threads for i in eachindex(timesteps)
        s_pos = @scenario t_xi = timesteps[i] s_xi =  1 F_xi = F
        s_neg = @scenario t_xi = timesteps[i] s_xi = -1 F_xi = F
        cost_pos[i] = evaluate_decision_wrapper(sp, decision, s_pos)
        cost_neg[i] = evaluate_decision_wrapper(sp, decision, s_neg)
    end
    return cost_pos, cost_neg
end

function evaluate_cost_penalty(sp, t, s, F, penalty_cost; decision = optimal_decision(sp))
    scen = @scenario t_xi = t s_xi = s F_xi = F
    outcome = outcome_model(sp, decision, scen; optimizer = subproblem_optimizer(sp))
    set_silent(outcome)
    optimize!(outcome)
    cost = objective_value(outcome)
    penalty = (maximum([value(outcome[:gi1]), value(outcome[:gi2])]) + maximum([value(outcome[:go1]), value(outcome[:go2])]))*penalty_cost
    return cost, penalty
end

#= We use investment and operation dictionaries rather than getting the data directly from the model, 
as this function is meant for post-processing
=#
function naive_flex_potential(invs, ops, pv, wind, timesteps)
    pv_cur = ops[:pv_cur]
    wind_cur = ops[:wind_cur]
    u_storage = invs[:u_storage]
    u_pv = invs[:u_pv]
    u_wind = invs[:u_wind]
    sto_soc = ops[:sto_soc]
    F_pos = zeros(length(timesteps[1:end-12]))
    F_neg = zeros(length(timesteps[1:end-12]))
    for t in timesteps[1:end-12]
        # positive flexibility:
        F_pos[t] = pv_cur[t] + wind_cur[t] + sto_soc[t+1]
        # negative flexibility:
        F_neg[t] = pv_cur[t]  + wind_cur[t] - pv[t]*u_pv - wind[t]*u_wind + sto_soc[t+1] - u_storage
    end
    return F_pos, F_neg.* (-1.)
end