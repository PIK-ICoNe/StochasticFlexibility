using Plots

"""
Suppress output ov evaluate_decision, and deal with solvers that error on infeasibility.
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
Find maximum available fleibility and its cost at all points in some time interval.
"""
function analyze_flexibility_potential(sp, timesteps)
    decision = optimal_decision(sp)
    cost_no_flex = evaluate_decision_wrapper(sp, optimal_decision(sp), no_flex_pseudo_sampler()[1])
    L = length(timesteps)
    cost_pos_flex = zeros(L)
    cost_neg_flex = zeros(L)
    potential_pos_flex = zeros(L)
    potential_neg_flex = zeros(L)
    Threads.@threads for i in 1:length(timesteps)
        t = timesteps[i]
        potential_pos_flex[i], cost_pos_flex[i] = find_f_max(sp,t,1,decision,100.)
        potential_neg_flex[i], cost_neg_flex[i] = find_f_max(sp,t,-1,decision,100.)
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
    display(plot(plt_cost, plt_pot, layout = (2, 1)))
end

"""
Find distribution of maximum available flexibilities.
"""
function flexibility_availability!(plt, flex_potential; plot_options...)
    p_sorted = sort(unique(flex_potential))
    fraction = 1 .-[sum(abs.(flex_potential[:]).<=abs(f)) for f in p_sorted]./length(flex_potential)
    plot!(plt, p_sorted, fraction)
end

"""
Bisection search for maximum felxibility potential.
"""
function find_f_max(p,t,s,od,tol; maxiter = 100)
    a = 0.
    b = 10000.
    cost_a = 0.
    cost_b = 0.
    i = 0
    scen = @scenario t_xi = t s_xi = s F_xi = a probability = 1.
    cost_a = evaluate_decision_wrapper(p,od,scen)
    if cost_a == Inf
        return 0., cost_a
    else
        while b-a>tol && i<=maxiter
            i+=1
            scen = @scenario t_xi = t s_xi = s F_xi = b probability = 1.
            cost_b = evaluate_decision_wrapper(p,od,scen)
            if cost_b == Inf
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

