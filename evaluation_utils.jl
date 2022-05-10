using Plots

function evaluate_decision_wrapper(p, decision, scenario)
    cost = 0.
    try redirect_stdout((() -> cost = evaluate_decision(p, decision, scenario)), devnull)
    catch e
        cost = Inf
    end
    return cost
end

function analyze_flexibility_potential(p, timesteps)
    decision = optimal_decision(p)
    L = length(timesteps)
    cost_pos_flex = zeros(L)
    cost_neg_flex = zeros(L)
    potential_pos_flex = zeros(L)
    potential_neg_flex = zeros(L)
    Threads.@threads for t in timesteps
        potential_pos_flex[t], cost_pos_flex[t] = find_f_max(p,t,1,decision,10.)
        potential_neg_flex[t], cost_neg_flex[t] = find_f_max(p,t,-1,decision,10.)
    end
    return cost_pos_flex, potential_pos_flex, cost_neg_flex, potential_neg_flex
end

function plot_flexibility(timesteps, cost_pos_flex, potential_pos_flex, cost_neg_flex, potential_neg_flex, obj_value)
    plt_cost = plot()
    plt_pot = plot()
    plot!(plt_cost, timesteps, (cost_pos_flex .- obj_value)./ potential_pos_flex, label = "price of positive flexibility")
    plot!(plt_cost, timesteps, (cost_neg_flex .- obj_value)./ potential_neg_flex, label = "price of negative flexibility")
    plot!(plt_pot, timesteps, potential_pos_flex, fillrange = 0, fillalpha = 0.35, label = "positive flexibility potential")
    plot!(plt_pot, timesteps, potential_neg_flex, fillrange = 0, fillalpha = 0.35, label = "negative flexibility potential")
    display(plot(plt_cost, plt_pot, layout = (2, 1)))
end

function flexibility_availability!(plt, flex_potential; plot_options...)
    p_sorted = sort(unique(flex_potential))
    fraction = 1 .-[sum(abs.(flex_potential[:]).<=abs(f)) for f in p_sorted]./length(flex_potential)
    plot!(plt, p_sorted, fraction)
end

function find_f_max(p,t,s,od,tol; maxiter = 100)
    a = 0.
    b = 10000.
    i = 0
    scen = @scenario t_xi = t s_xi = s F_xi = 0. probability = 1.
    cost = evaluate_decision_wrapper(p,od,scen)
    if cost == Inf
        return 0.,cost
    else
        while b-a>tol && i<=maxiter
            i+=1
            scen = @scenario t_xi = t s_xi = s F_xi = b probability = 1.
            cost = evaluate_decision_wrapper(p,od,scen)
            if cost == Inf
                b = (a+b)/2
            else
                a = b
                b *= 2
            end
        end
        if cost == Inf
            scen = @scenario t_xi = t s_xi = s F_xi = b-tol probability = 1.
            cost = evaluate_decision_wrapper(p,od,scen)
        end
        return s*b,cost
    end
end

