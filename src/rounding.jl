"""
    simplerounding!(xTilde::Vector{Float64}, xOverline::Vector{Float64},
        timeStart::UInt, timeLim::Float64, indices::Vector{Int})

The simple rounding round every variable to the nearest binary.
"""
function simplerounding!(
        xTilde::Vector{Float64},
        xOverline::Vector{Float64},
        timeStart::UInt,
        timeLim::Float64,
        indices::Vector{Int},
        opts...
        ;
        compile::Bool = false
    )
    if compile
        return [NaN]
    end

    xTilde .= xOverline
    for i in indices
        xTilde[i] = round(xTilde[i])
    end
    return xTilde
end


"""
    recursiverounding!(xTilde::Vector{Float64}, xOverline::Vector{Float64},
        timeStart::UInt, timeLim::Float64, _indices::Vector{Int},
        myModel::CPLEX.Model, indices::Vector{Int},
        objectiveFunction::Vector{Float64})

The recursive rounding round to the nearest. If that's not a feasible solution
we try to found one fixing less and less variables. Then we round it to the
nearest. If that's still not a solution less variables are fixed and try again.
    See https://doi.org/10.1016/j.cor.2013.09.008
"""
function recursiverounding!(
        xTilde::Vector{Float64},
        xOverline::Vector{Float64},
        timeStart::UInt,
        timeLim::Float64,
        _indices::Vector{Int},
        myModel::CPLEX.Model,
        indices::Vector{Int},
        objectiveFunction::Vector{Float64},
        opts...
        ;
        compile::Bool = false
    )
    if compile
        return [NaN]
    end

    timeStart = time_ns()
    xTilde .= xOverline
    isFeasible = false
    first = 1
    nVars = length(indices)
    last = nVars
    CPLEX.set_obj!(myModel, objectiveFunction)
    varLB_init = CPLEX.get_varLB(myModel)
    varUB_init = CPLEX.get_varUB(myModel)
    varLB = copy(varLB_init)
    varUB = copy(varUB_init)
    for i in indices
        xTilde[i] = round(xOverline[i])
    end

    for i in indices[first:last]
        varLB[i] = xTilde[i]
        varUB[i] = varLB[i]
    end
    while !isFeasible && first < last && (time_ns()-timeStart)/1.0e9 < timeLim
        CPLEX.set_varLB!(myModel, varLB)
        CPLEX.set_varUB!(myModel, varUB)
        CPLEX.set_param!(myModel.env, "CPX_PARAM_TILIM", maximum([0, timeLim - (time_ns() - timeStart) / 1.0e9]))
        CPLEX.optimize!(myModel)

        if CPLEX.get_status(myModel) == :CPX_STAT_OPTIMAL
            isFeasible = true
        end

        if !isFeasible
            lastBis = last - round(Int, (last - first) / 2, RoundUp)
            for i in indices[(lastBis+1):last]
                varLB[i] = varLB_init[i]
                varUB[i] = varUB_init[i]
            end
            last = lastBis
        elseif last != nVars
            isFeasible = false
            first = last+1
            last = nVars
            x = CPLEX.get_solution(myModel)
            for i in indices[(first+1):end]
                xTilde[i] = round(x[i])
                varLB[i] = xTilde[i]
                varUB[i] = varLB[i]
            end
        end
    end

    CPLEX.set_varLB!(myModel, varLB_init)
    CPLEX.set_varUB!(myModel, varUB_init)
    return xTilde
end

"""
    sortReducedCosts!(indices::Vector{Int}, myModel::CPLEX.Model)

sortReducedCosts sort indices by decreasing absolute reduced costs for myModel.
"""
function sortReducedCosts!(
        indices::Vector{Int},
        myModel::CPLEX.Model
        ;
        compile = false
    )
    if compile
        return [0]
    end

    reducedCosts = Vector{Tuple{Int, Float64}}()
    reducedCostsInit = CPLEX.get_reduced_costs(myModel)
    for i in indices
        push!(reducedCosts, (i, abs(reducedCostsInit[i])))
    end
    sort!(reducedCosts, by = x -> x[2], rev=false)
    indices = first.(reducedCosts)
    return indices
end
