function curve_fit(
    x::Vector{Float64},
    y::Vector{Float64},
)
    # Fit a polynomial to the data
    p = Polynomials.fit(x, y, 1)  # linear regression
    return p
end

function fit_objective_values(cluster::Cluster, rep_days_costs::Dict{Int,Float64})
    x_group = [cluster.days[i].fvalue for i in cluster.repdays]
    y_group = [rep_days_costs[i] for i in cluster.repdays]
    polynomial_curve = curve_fit(x_group, y_group)

    return polynomial_curve
end

function cluster_cost_regression(cluster::Cluster, rep_days_costs::Dict{Int,Float64})
    polynomial_curve = fit_objective_values(cluster, rep_days_costs)
    cost_per_day = [polynomial_curve(i.fvalue) for i in values(cluster.days)]
    adjusted_cost_per_day = map(x -> x < 0 ? 0 : x, cost_per_day)
    cost::Float64 = sum(adjusted_cost_per_day)

    return cost
end

function total_year_oper_cost(clusters::Dict{String,Cluster}, costs::Dict{Int,Float64})
    cost = 0.0
    for cluster in values(clusters)
        if length(cluster.days) <= 3
            cost += sum(costs[i] for i in cluster.repdays)
        else 
            cost += cluster_cost_regression(cluster, costs)
        end
    end

    return cost
end