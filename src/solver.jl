mutable struct ParametersFP
    rounding_method::Int
    projection_method::Int
    perturbation_method::Int
    restart_method::Int
    presolve::Int
    improve::Bool
    recursive_rounding::Int
    TT::Bool
    TTmin::Int
    TTmax::Int
    feasibility_check::Int # 0: both check, 1: only isfeasible_xOverline, 2: only isfeasible_xTilde
    restart::Int
    check_cycle::Int

    function ParametersFP()
        parameter = new()
        parameter.rounding_method = 0
        parameter.projection_method = 0
        parameter.perturbation_method = 0
        parameter.restart_method = 0
        parameter.presolve = 0
        parameter.improve = false
        parameter.recursive_rounding = 0
        parameter.TT = true
        parameter.TTmin = 10
        parameter.TTmax = 30
        parameter.feasibility_check = 1
        parameter.restart = 100
        parameter.check_cycle = 3
        return parameter
    end
end

#-------------------------------------------------------------------------------
#-----
#-------------------------------------------------------------------------------
function feasibilitypump(
        model::Model,
        parameter::ParametersFP
    )
    starting_time = time_ns()

    rounding_method! = which_rounding[0]
    projection_method = which_projection[0]
    perturbmethod! = which_perturb[0]
    restart_method! = which_restart[0]
    try
        rounding_method! = which_rounding[parameter.rounding_method]
    catch
        @warn("Invalid parameter for rounding_method, default used")
    end
    try
        projection_method = which_projection[parameter.projection_method]
    catch
        @warn("Invalid parameter for projection_method, default used")
    end
    try
        perturbmethod! = which_perturb[parameter.perturbation_method]
    catch
        @warn("Invalid parameter for perturbation_method, default used")
    end
    try
        restart_method! = which_restart[parameter.restart_method]
    catch
        @warn("Invalid parameter for perturbation_method, default used")
    end
    restart = parameter.restart
    restart_step = 0
    check_cycle = parameter.check_cycle
    presolve = parameter.presolve
    improve = parameter.improve
    total_time_limit = time_limit_sec(model)
    if total_time_limit === nothing
        total_time_limit = Inf
    end
    feasibility_check = parameter.feasibility_check

    LPMethod_init = CPLEX.get_param(model.env, "CPX_PARAM_LPMETHOD")
    CPLEX.set_param!(model.env, "CPX_PARAM_LPMETHOD", 1)
    nb_variables = CPLEX.num_var(model)
    varTypes_init = CPLEX.get_vartype(model)
    varLB_init = CPLEX.get_varLB(model)
    varLB = copy(varLB_init)
    varUB_init = CPLEX.get_varUB(model)
    varUB = copy(varUB_init)
    indices = Vector{Int}()
    for i in 1:nb_variables
        if varTypes_init[i] == 'B'
            if varLB[i] != varUB[i]
                push!(indices, i)
            else
                # Fixed variable
            end
        end
    end

    saveSol = Vector{Vector{Float64}}()
    for i in 1:(check_cycle-1)
        push!(saveSol, zeros(Float64, nb_variables))
    end
    saveVect = 0

    initial_objective_function = CPLEX.get_obj(model)
    objective_function = copy(initial_objective_function)
    objSense_init = CPLEX.get_sense(model)
    if objSense_init == :Max
        objective_function = -objective_function
    end

    # This is only useful for recursive_rounding!
    sortBinVariables = false
    sortBinVariablesOnce = false
    objectiveFunctionRR = zeros(Float64, nb_variables)
    if parameter.recursive_rounding == 1
        sortBinVariables = true
        sortBinVariablesOnce = true
    elseif parameter.recursive_rounding == 2
        sortBinVariables = true
        sortBinVariablesOnce = false
    elseif parameter.recursive_rounding == 3
        sortBinVariables = true
        sortBinVariablesOnce = true
        objectiveFunctionRR .= initial_objective_function
    elseif parameter.recursive_rounding == 4
        sortBinVariables = true
        sortBinVariablesOnce = false
        objectiveFunctionRR .= initial_objective_function
    elseif parameter.recursive_rounding == 0
        #Nothing
    else
        @warn "Parameter recursive_rounding has an unauthorized value, parameter ignored"
    end

    #=--------------------------------------------------------------------------
    Variables needed for presolve and to check feasibility
    --------------------------------------------------------------------------=#
    A = CPLEX.get_constr_matrix(model)
    senses = CPLEX.get_constr_senses(model)
    rhs = CPLEX.get_rhs(model)
    nConstr = length(rhs)

    if presolve != 0
        nBin = length(indices)
        indicesTmp = falses(nBin)
        for ind in 1:indices
            i = indices[ind]
            fixed = true
            val = -1
            for j in A.rowval[nzrange(A, i)]
                if senses[j] == 76
                    if val == -1
                        if A[j,i] < 0
                            val = 1
                        elseif A[j,i] > 0
                            val = 0
                        end
                    elseif val == 0
                        if A[j,i] < 0
                            fixed = false
                            break
                        end
                    elseif val == 1
                        if A[j,i] > 0
                            fixed = false
                            break
                        end
                    end
                elseif senses[j] == 71
                    if val == -1
                        if A[j,i] > 0
                            val = 1
                        elseif A[j,i] < 0
                            val = 0
                        end
                    elseif val == 0
                        if A[j,i] > 0
                            fixed = false
                            break
                        end
                    elseif val == 1
                        if A[j,i] < 0
                            fixed = false
                            break
                        end
                    end
                elseif senses[j] == 69
                    fixed = false
                    break
                end
            end
            if fixed
                if presolve == 1 || (val == 0 && objective_function[i] >= 0) || (val == 1 && objective_function[i] <= 0)
                    varLB[i] = val
                    varUB[i] = val
                    indicesTmp[ind] = true
                end
            end
        end
        CPLEX.set_varLB!(model, varLB)
        CPLEX.set_varUB!(model, varUB)
        deleteat!(indices, indicesTmp)
    end

    indicesRR = copy(indices)
    A = sparse(A') # This transpose makes feasibility check much faster

    if parameter.TT
        TTmin = parameter.TTmin
        TTmax = parameter.TTmax
    else
        TTmin = round(Int, length(indices)*parameter.TTmin/100, RoundDown)
        TTmax = round(Int, length(indices)*parameter.TTmax/100, RoundUp)
    end
    if TTmin > length(indices)
        TTmin = length(indices)
    end
    if TTmax > length(indices)
        TTmax = length(indices)
    end
    if TTmin > TTmax
        @warn "TTmin is greater than TTmax, parameters are switched"
        TTmin = TTmax
        TTmax = parameter.TTmin
    end

    problemType_init = CPLEX.get_prob_type(model)
    CPLEX.set_prob_type!(model, 0)
    has_int_init = model.has_int
    model.has_int = false

    #=--------------------------------------------------------------------------
    We initialize variables for variants:
      coef, alpha and alphaChange are parameters of objective
      feasibility pump. freqRound is useful for frequenciesrounding
    --------------------------------------------------------------------------=#
    coef = norm(objective_function)
    if coef > 1e-8
        coef = sqrt(length(indices))/coef
    end
    alphaChange = 0.9
    alpha = 1.0

    freqRound = zeros(Int, nb_variables)

    #=--------------------------------------------------------------------------
    First linear relaxation, rounding and feasibility check
    --------------------------------------------------------------------------=#
    #CPLEX.set_obj!(model, zeros(Float64, nb_variables))
    set_time_limit!(model, total_time_limit, starting_time)
    optimize!(model)

    xOverline = zeros(Float64, nb_variables)
    if CPLEX.get_status(model) == :CPX_STAT_OPTIMAL
        xOverline = CPLEX.get_solution(model)
    end

    isFeasible = false
    if feasibility_check != 2
        isFeasible = isfeasible_xOverline(xOverline, indices, starting_time, total_time_limit)
    end
    if time_limit_sec(model) == 0.0
        return false
    end
    nIter = 0
    xTilde = zeros(Float64, nb_variables)

    if !isFeasible
        if sortBinVariables
            sortReducedCosts!(indicesRR, model)
            if sortBinVariablesOnce
                sortBinVariables = false
            end
        end

        CPLEX.set_sense!(model, :Min)

        rounding_method!(xTilde, xOverline, starting_time, total_time_limit, indices, model, indicesRR, objectiveFunctionRR)
        xTildeBis = copy(xTilde)
        for i in indices
            if xTilde[i] > 0.5
                freqRound[i] += 1
            end
        end
        if feasibility_check != 1
            isFeasible = isfeasible_xTilde(xTilde, A, senses, rhs, starting_time, total_time_limit)
        end
    else
        xTilde .= xOverline
    end

    #=--------------------------------------------------------------------------
    The feasibility pump
    --------------------------------------------------------------------------=#
    while !isFeasible && (time_ns()-starting_time)/1.0e9 < total_time_limit
        nIter += 1
        restart_step += 1
        alpha *= alphaChange
        xOverline = projection_method(model, starting_time, total_time_limit, xTilde, indices, objective_function, alpha, coef)
        feasibility_check != 2
            isFeasible = isfeasible_xOverline(xOverline, indices, starting_time, total_time_limit)
        end

        if !isFeasible
            if sortBinVariables
                sortReducedCosts!(indicesRR, model)
            end

            rounding_method!(xTildeBis, xOverline, starting_time, total_time_limit, indices, model, indicesRR, objectiveFunctionRR)
            if restart_step == restart
                restart_step = 0
                restart_method!(xTildeBis, xOverline, indices)
            end
            if xTildeBis == xTilde
                perturbmethod!(xTildeBis, xOverline, indices, TTmin, TTmax, freqRound, nIter)
            else
                if xTildeBis in saveSol
                    restart_step = 0
                    restart_method!(xTildeBis, xOverline, indices)
                end
            end
            if check_cycle > 1
                saveSol[saveVect+1] .= xTildeBis
                saveVect = mod(saveVect+1, check_cycle-1)
            end
            xTilde = copy(xTildeBis)

            for i in indices
                if xTilde[i] > 0.5
                    freqRound[i] += 1
                end
            end
            if feasibility_check != 1
                isFeasible = isfeasible_xTilde(xTilde, A, senses, rhs, starting_time, total_time_limit)
            end
        else
            xTilde .= xOverline
        end
    end

    objVal = NaN
    status = :TIMEOUT
    if isFeasible
        status = :SOL_FOUND
        objVal = sum(initial_objective_function .* xTilde)
        if improve
            CPLEX.set_sense!(model, objSense_init)
            CPLEX.set_obj!(model, initial_objective_function)
            varLB = CPLEX.get_varLB(model)
            varUB = CPLEX.get_varUB(model)
            for i in indices
                varLB[i] = round(xTilde[i])
                varUB[i] = varLB[i]
            end
            CPLEX.set_varLB!(model, varLB)
            CPLEX.set_varUB!(model, varUB)

            set_time_limit!(model, total_time_limit, starting_time)
            CPLEX.optimize!(model)
            if CPLEX.get_status(model) == :CPX_STAT_OPTIMAL
                xTilde = CPLEX.get_solution(model)
                objVal = CPLEX.get_objval(model)
            end
        end
    else
        xTilde = repeat([NaN], nb_variables)
    end

    #=--------------------------------------------------------------------------
    We changed the model to solve it, now we repair it
    --------------------------------------------------------------------------=#
    CPLEX.set_param!(model.env, "CPX_PARAM_LPMETHOD", LPMethod_init)
    CPLEX.set_vartype!(model, varTypes_init)
    CPLEX.set_vartype!(model, varTypes_init)
    CPLEX.set_varLB!(model, varLB_init)
    CPLEX.set_varUB!(model, varUB_init)
    CPLEX.set_obj!(model, initial_objective_function)
    CPLEX.set_prob_type!(model, problemType_init)
    CPLEX.set_sense!(model, objSense_init)
    model.has_int = has_int_init
    CPLEX.set_param!(model.env, "CPX_PARAM_TILIM", time_limit_init)

    return status, xTilde, objVal
end

feasibility_pump(model::Model) =
    feasibility_pump(model, ParametersFP())



function isfeasible_xTilde(
        xTilde::Vector{Float64},
        A::SparseArrays.SparseMatrixCSC{Float64,Int64},
        senses::Vector{Int8},
        rhs::Vector{Float64},
        starting_time::UInt,
        total_time_limit::Float64
    )
    nb_variables::Int64, nConstr::Int64 = size(A)
    isFeasible = true
    sum::Float64 = 0.0
    i = 1
    while isFeasible && i <= nConstr
    #for i in 1:nConstr
        sum = 0.0
        #for j in 1:nb_variables
        for j in A.rowval[nzrange(A, i)]
            sum += A[j,i] * xTilde[j]
        end
        if senses[i] == 76
            if sum > rhs[i]
                isFeasible = false
            end
        elseif senses[i] == 71
            if sum < rhs[i]
                isFeasible = false
            end
        elseif senses[i] == 69
            if !isapprox(sum, rhs[i], atol=1e-8)
                isFeasible = false
            end
        end
        i += 1
        set_time_limit!(model, total_time_limit, starting_time)
        if time_limit_sec(model) == 0.0
            return false
        end
    end

    return isFeasible
end

function isfeasible_xOverline(
        xOverline::Vector{Float64},
        indices::Vector{Int},
        starting_time::UInt,
        total_time_limit::Float64
    )
    nBin = length(indices)
    ind = 1
    isFeasible = true
    while ind <= nBin && isFeasible
        i = indices[ind]
        if abs(xOverline[i]-round(xOverline[i])) > 1e-8
            isFeasible = false
        end
        set_time_limit!(model, total_time_limit, starting_time)
        if time_limit_sec(model) == 0.0
            return false
        end
        ind += 1
    end

    return isFeasible
end


function initialize_parameters_fp()
    return ParametersFP()
end


function setparameterfp!(parameters_list::ParametersFP, parameter::Int, value::Int)
    if parameter == 10
        parameters_list.rounding_method = value
    elseif parameter == 20
        parameters_list.projection_method = value
    elseif parameter == 30
        parameters_list.perturbation_method = value
    elseif parameter == 32
        parameters_list.TTmin = value
    elseif parameter == 33
        parameters_list.TTmax = value
    elseif parameter == 34
        parameters_list.check_cycle = value
    elseif parameter == 35
        parameters_list.restart = value
    elseif parameter == 2
        parameters_list.time_limit = Float64(value)
    elseif parameter == 11
        parameters_list.recursive_rounding = value
    elseif parameter == 41
        parameters_list.presolve = value
    elseif parameter == 100
        parameters_list.feasibility_check = value
    elseif parameter == 31
        try
            parameters_list.TT = Bool(value)
        catch
            @warn("The parameter $parameter can take only 0 or 1, nothing happened")
        end
    elseif parameter == 42
        try
            parameters_list.improve = Bool(value)
        catch
            @warn("The parameter $parameter can take only 0 or 1, nothing happened")
        end
    else
        @warn("The parameter $parameter does not exist, nothing happened")
    end
    return parameters_list
end

function setparamfp!(parameters_list::ParametersFP, parameter::Int, value::Float64)
    if parameter == 10
        parameters_list.rounding_method = round(Int, value)
    elseif parameter == 20
        parameters_list.projection_method = round(Int, value)
    elseif parameter == 30
        parameters_list.perturbation_method = round(Int, value)
    elseif parameter == 32
        parameters_list.TTmin = round(Int, value)
    elseif parameter == 33
        parameters_list.TTmax = round(Int, value)
    elseif parameter == 34
        parameters_list.check_cycle = round(Int, value)
    elseif parameter == 35
        parameters_list.restart = round(Int, value)
    elseif parameter == 11
        parameters_list.recursive_rounding = round(Int, value)
    elseif parameter == 41
        parameters_list.presolve = round(Int, value)
    elseif parameter == 100
        parameters_list.feasibility_check = round(Int, value)
    elseif parameter == 31
        try
            parameters_list.TT = Bool(round(Int, value))
        catch
            @warn("The parameter $parameter can take only 0 or 1, nothing happened")
        end
    elseif parameter == 42
        try
            parameters_list.improve = Bool(round(Int, value))
        catch
            @warn("The parameter $parameter can take only 0 or 1, nothing happened")
        end
    else
        @warn("The parameter $parameter does not exist, nothing happened")
    end
    return parameters_list
end


function setparamfp!(parameters_list::ParametersFP, parameter::String, value::Union{Float64, Int})
    if haskey(__paramToInt, parameter)
        setparamfp!(parameters_list, __paramToInt[parameter], value)
    else
        @warn("The parameter $parameter does not exist, nothing happened")
    end
    return parameters_list
end


function getparamfp(parameters_list::ParametersFP, parameter::Int)
    value = NaN
    if parameter == 10
        value = parameters_list.rounding_method
    elseif parameter == 20
        value = parameters_list.projection_method
    elseif parameter == 30
        value = parameters_list.perturbation_method
    elseif parameter == 31
        value = parameters_list.TT
    elseif parameter == 32
        value = parameters_list.TTmin
    elseif parameter == 33
        value = parameters_list.TTmax
    elseif parameter == 34
        value = parameters_list.check_cycle
    elseif parameter == 35
        value = parameters_list.restart
    elseif parameter == 11
        value = parameters_list.recursive_rounding
    elseif parameter == 41
        value = parameters_list.presolve
    elseif parameter == 100
        value = parameters_list.feasibility_check
    elseif parameter == 42
        value = parameters_list.improve
    else
        @warn("The parameter $parameter does not exist, nothing happened")
    end
    return value
end


function get_parameter_fp(parameters_list::ParametersFP, parameter::String)
    value = NaN
    if haskey(__paramToInt, parameter)
        value = getparamfp(paramList, __paramToInt[parameter])
    end
    return value
end

function get_parameter_name_fp(parameter::Int)
    paramName = ""
    if haskey(__paramToString, parameter)
        paramName = __paramToString[parameter]
    else
        @warn("No parameter $parameter")
    end
    return paramName
end
