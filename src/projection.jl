"""
    minimumdistance(myModel::CPLEX.Model, timeStart::UInt, timeLim::Float64,
        xTilde::Vector{Float64}, indices::Vector{Int})

Replace the objective function to minimize the L1-distance to xTilde. Then
optimize to return xOverline.
"""
function minimumdistance(
        myModel::CPLEX.Model,
        timeStart::UInt,
        timeLim::Float64,
        xTilde::Vector{Float64},
        indices::Vector{Int},
        opts...
        ;
        compile::Bool = false
    )
    if compile
        return [NaN]
    end

    nVars = length(xTilde)
    objectiveFunction = zeros(Float64, nVars)
    for i in indices
        if xTilde[i] > 0.5
            objectiveFunction[i] = -1.0
        else
            objectiveFunction[i] = 1.0
        end
    end
    CPLEX.set_obj!(myModel, objectiveFunction)

    CPLEX.set_param!(myModel.env, "CPX_PARAM_TILIM", maximum([0, timeLim - (time_ns() - timeStart) / 1.0e9]))
    CPLEX.optimize!(myModel)

    xOverline = zeros(Float64, nVars)
    if CPLEX.get_status(myModel) == :CPX_STAT_OPTIMAL
        xOverline = CPLEX.get_solution(myModel)
    end

    return xOverline
end


"""
    minimumdeltaalpha(myModel::CPLEX.Model, timeStart::UInt, timeLim::Float64,
        xTilde::Vector{Float64}, indices::Vector{Int},
        initObjectiveFunction::Vector{Float64}, alpha::Float64, coef::Float64)

Replace the objective function to minimize the L1-distance to xTilde taking into
account the initial objective function. Then optimize to return xOverline.
"""
function minimumdeltaalpha(
        myModel::CPLEX.Model,
        timeStart::UInt,
        timeLim::Float64,
        xTilde::Vector{Float64},
        indices::Vector{Int},
        initObjectiveFunction::Vector{Float64},
        alpha::Float64,
        coef::Float64,
        opts...
        ;
        compile::Bool = false
    )
    if compile
        return [NaN]
    end

    nVars = length(xTilde)

    # objectiveFunction minimumdistance
    objectiveFunction = zeros(Float64, nVars)
    for i in indices
        if xTilde[i] > 0.5
            objectiveFunction[i] = -1.0
        else
            objectiveFunction[i] = 1.0
        end
    end

    # Add objective function and alpha to it
    objectiveFunction .= (objectiveFunction .* (1-alpha)) .+ ((alpha*coef) .* initObjectiveFunction)
    CPLEX.set_obj!(myModel, objectiveFunction)

    CPLEX.set_param!(myModel.env, "CPX_PARAM_TILIM", maximum([0, timeLim - (time_ns() - timeStart) / 1.0e9]))
    CPLEX.optimize!(myModel)

    xOverline = zeros(Float64, nVars)
    if CPLEX.get_status(myModel) == :CPX_STAT_OPTIMAL
        xOverline = CPLEX.get_solution(myModel)
    end

    return xOverline
end


"""
    minimumdistance_reduce(myModel::CPLEX.Model, timeStart::UInt,
        timeLim::Float64, xTilde::Vector{Float64}, indices::Vector{Int})

Replace the objective function to minimize the L1-distance to xTilde. Multiply
the new objective function by the non zeros reduced costs. Then optimize to
return xOverline.
"""
function minimumdistance_reduce(
        myModel::CPLEX.Model,
        timeStart::UInt,
        timeLim::Float64,
        xTilde::Vector{Float64},
        indices::Vector{Int},
        opts...
        ;
        compile::Bool = false
    )
    if compile
        return [NaN]
    end

    nVars = length(xTilde)
    reducedCosts = Vector{Float64}()
    try
        reducedCosts = CPLEX.get_reduced_costs(myModel)
    catch
        reducedCosts = zeros(nVars)
    end
    objectiveFunction = zeros(Float64, nVars)
    for i in indices
        if xTilde[i] > 0.5
            objectiveFunction[i] = -1.0
        else
            objectiveFunction[i] = 1.0
        end
        if abs(reducedCosts[i]) > 1e-8
            objectiveFunction[i] *= abs(reducedCosts[i])
        end
    end
    CPLEX.set_obj!(myModel, objectiveFunction)

    CPLEX.set_param!(myModel.env, "CPX_PARAM_TILIM", maximum([0, timeLim - (time_ns() - timeStart) / 1.0e9]))
    CPLEX.optimize!(myModel)

    xOverline = zeros(Float64, nVars)
    if CPLEX.get_status(myModel) == :CPX_STAT_OPTIMAL
        xOverline = CPLEX.get_solution(myModel)
    end

    return xOverline
end
