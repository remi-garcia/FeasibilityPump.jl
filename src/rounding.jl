"""
    simplerounding!(x_tilde::Vector{Float64}, x_overline::Vector{Float64},
        time_start::UInt, total_time_limit::Float64, indices::Vector{Int})

The simple rounding round every variable to the nearest binary.
"""
function simple_rounding!(
        x_tilde::Vector{Float64},
        x_overline::Vector{Float64},
        time_start::UInt,
        total_time_limit::Float64,
        indices::Vector{Int},
        opts...
    )
    x_tilde .= x_overline
    for i in indices
        x_tilde[i] = round(x_tilde[i])
    end
    return x_tilde
end


"""
    recursive_rounding!(x_tilde::Vector{Float64}, x_overline::Vector{Float64},
        time_start::UInt, total_time_limit::Float64, _indices::Vector{Int},
        model::Model, indices::Vector{Int},
        objective_function::Vector{Float64})

The recursive rounding round to the nearest. If that's not a feasible solution
we try to found one fixing less and less variables. Then we round it to the
nearest. If that's still not a solution less variables are fixed and try again.
    See https://doi.org/10.1016/j.cor.2013.09.008
"""
function recursive_rounding!(
        x_tilde::Vector{Float64},
        x_overline::Vector{Float64},
        time_start::UInt,
        total_time_limit::Float64,
        _indices::Vector{Int},
        model::Model,
        indices::Vector{Int},
        objective_function::Vector{Float64},
        opts...
    )
    time_start = time_ns()
    x_tilde .= x_overline
    is_feasible = false
    first = 1
    nb_variables = length(indices)
    last = nb_variables
    CPLEX.set_obj!(model, objective_function)
    varLB_init = CPLEX.get_varLB(model)
    varUB_init = CPLEX.get_varUB(model)
    varLB = copy(varLB_init)
    varUB = copy(varUB_init)
    for i in indices
        x_tilde[i] = round(x_overline[i])
    end

    for i in indices[first:last]
        varLB[i] = x_tilde[i]
        varUB[i] = varLB[i]
    end
    while !is_feasible && first < last && (time_ns()-time_start)/1.0e9 < total_time_limit
        CPLEX.set_varLB!(model, varLB)
        CPLEX.set_varUB!(model, varUB)
        CPLEX.set_param!(model.env, "CPX_PARAM_TILIM", maximum([0, total_time_limit - (time_ns() - time_start) / 1.0e9]))
        CPLEX.optimize!(model)

        if CPLEX.get_status(model) == :CPX_STAT_OPTIMAL
            is_feasible = true
        end

        if !is_feasible
            lastBis = last - round(Int, (last - first) / 2, RoundUp)
            for i in indices[(lastBis+1):last]
                varLB[i] = varLB_init[i]
                varUB[i] = varUB_init[i]
            end
            last = lastBis
        elseif last != nb_variables
            is_feasible = false
            first = last+1
            last = nb_variables
            x = CPLEX.get_solution(model)
            for i in indices[(first+1):end]
                x_tilde[i] = round(x[i])
                varLB[i] = x_tilde[i]
                varUB[i] = varLB[i]
            end
        end
    end

    CPLEX.set_varLB!(model, varLB_init)
    CPLEX.set_varUB!(model, varUB_init)
    return x_tilde
end

"""



"""
function get_all_reduced_costs(model::Model)


    return reduced_costs
end

"""
    sort_reduced_costs!(indices::Vector{Int}, model::Model)

Sort indices by decreasing absolute reduced costs for model.
"""
function sort_reduced_costs!(indices::Vector{Int}, model::Model)
    reduced_costs = Vector{Tuple{Int, Float64}}()
    reduced_costs_init = get_all_reduced_costs(model)
    for i in indices
        push!(reduced_costs, (i, abs(reduced_costs_init[i])))
    end
    sort!(reduced_costs, by = x -> x[2], rev=false)
    indices = first.(reduced_costs)
    return indices
end
