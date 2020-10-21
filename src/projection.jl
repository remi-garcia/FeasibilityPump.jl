"""
    minimumdistance(model::CPLEX.Model, time_start::UInt, total_time_limit::Float64,
        x_tilde::Vector{Float64}, indices::Vector{Int})

Replace the objective function to minimize the L1-distance to x_tilde. Then
optimize to return x_overline.
"""
function minimumdistance(
        model::CPLEX.Model,
        time_start::UInt,
        total_time_limit::Float64,
        x_tilde::Vector{Float64},
        indices::Vector{Int},
        opts...
    )
    nb_variables = length(x_tilde)
    objective_function = zeros(Float64, nb_variables)
    for i in indices
        if x_tilde[i] > 0.5
            objective_function[i] = -1.0
        else
            objective_function[i] = 1.0
        end
    end
    CPLEX.set_obj!(model, objective_function)

    CPLEX.set_param!(model.env, "CPX_PARAM_TILIM", maximum([0, total_time_limit - (time_ns() - time_start) / 1.0e9]))
    CPLEX.optimize!(model)

    x_overline = zeros(Float64, nb_variables)
    if CPLEX.get_status(model) == :CPX_STAT_OPTIMAL
        x_overline = CPLEX.get_solution(model)
    end

    return x_overline
end


"""
    minimumdeltaalpha(model::CPLEX.Model, time_start::UInt, total_time_limit::Float64,
        x_tilde::Vector{Float64}, indices::Vector{Int},
        initObjectiveFunction::Vector{Float64}, alpha::Float64, coef::Float64)

Replace the objective function to minimize the L1-distance to x_tilde taking into
account the initial objective function. Then optimize to return x_overline.
"""
function minimumdeltaalpha(
        model::CPLEX.Model,
        time_start::UInt,
        total_time_limit::Float64,
        x_tilde::Vector{Float64},
        indices::Vector{Int},
        initObjectiveFunction::Vector{Float64},
        alpha::Float64,
        coef::Float64,
        opts...
    )
    nb_variables = length(x_tilde)

    # objective_function minimumdistance
    objective_function = zeros(Float64, nb_variables)
    for i in indices
        if x_tilde[i] > 0.5
            objective_function[i] = -1.0
        else
            objective_function[i] = 1.0
        end
    end

    # Add objective function and alpha to it
    objective_function .= (objective_function .* (1-alpha)) .+ ((alpha*coef) .* initObjectiveFunction)
    CPLEX.set_obj!(model, objective_function)

    CPLEX.set_param!(model.env, "CPX_PARAM_TILIM", maximum([0, total_time_limit - (time_ns() - time_start) / 1.0e9]))
    CPLEX.optimize!(model)

    x_overline = zeros(Float64, nb_variables)
    if CPLEX.get_status(model) == :CPX_STAT_OPTIMAL
        x_overline = CPLEX.get_solution(model)
    end

    return x_overline
end


"""
    minimumdistance_reduce(model::CPLEX.Model, time_start::UInt,
        total_time_limit::Float64, x_tilde::Vector{Float64}, indices::Vector{Int})

Replace the objective function to minimize the L1-distance to x_tilde. Multiply
the new objective function by the non zeros reduced costs. Then optimize to
return x_overline.
"""
function minimumdistance_reduce(
        model::CPLEX.Model,
        time_start::UInt,
        total_time_limit::Float64,
        x_tilde::Vector{Float64},
        indices::Vector{Int},
        opts...
    )
    nb_variables = length(x_tilde)
    reducedCosts = Vector{Float64}()
    try
        reducedCosts = CPLEX.get_reduced_costs(model)
    catch
        reducedCosts = zeros(nb_variables)
    end
    objective_function = zeros(Float64, nb_variables)
    for i in indices
        if x_tilde[i] > 0.5
            objective_function[i] = -1.0
        else
            objective_function[i] = 1.0
        end
        if abs(reducedCosts[i]) > 1e-8
            objective_function[i] *= abs(reducedCosts[i])
        end
    end
    CPLEX.set_obj!(model, objective_function)

    CPLEX.set_param!(model.env, "CPX_PARAM_TILIM", maximum([0, total_time_limit - (time_ns() - time_start) / 1.0e9]))
    CPLEX.optimize!(model)

    x_overline = zeros(Float64, nb_variables)
    if CPLEX.get_status(model) == :CPX_STAT_OPTIMAL
        x_overline = CPLEX.get_solution(model)
    end

    return x_overline
end
