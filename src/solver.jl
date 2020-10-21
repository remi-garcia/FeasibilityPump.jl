mutable struct ParametersFP
    roundingMethod::Int
    projectionMethod::Int
    perturbationMethod::Int
    restartMethod::Int
    presolve::Int
    improve::Bool
    logLevel::Bool
    timelim::Float64
    recursiveRounding::Int
    TT::Bool
    TTmin::Int
    TTmax::Int
    feasibilityCheck::Int # 0: both check, 1: only isfeasible_xOverline, 2: only isfeasible_xTilde
    restart::Int
    checkCycle::Int

    function ParametersFP()
        param = new()
        param.roundingMethod = 0
        param.projectionMethod = 0
        param.perturbationMethod = 0
        param.restartMethod = 0
        param.presolve = 0
        param.improve = false
        param.logLevel = false
        param.timelim = 30.0
        param.recursiveRounding = 0
        param.TT = true
        param.TTmin = 10
        param.TTmax = 30
        param.feasibilityCheck = 1
        param.restart = 100
        param.checkCycle = 3
        return param
    end
end

#-------------------------------------------------------------------------------
#-----
#-------------------------------------------------------------------------------
function feasibilitypump(
        model::CPLEX.Model,
        param::ParametersFP
    )
    timeStart = time_ns()

    roundingmethod! = whichRounding[0]
    projectionmethod = whichProjection[0]
    perturbmethod! = whichPerturb[0]
    restartmethod! = whichRestart[0]
    try
        roundingmethod! = whichRounding[param.roundingMethod]
    catch
        @warn("Invalid parameter for roundingMethod, default used")
    end
    try
        projectionmethod = whichProjection[param.projectionMethod]
    catch
        @warn("Invalid parameter for projectionMethod, default used")
    end
    try
        perturbmethod! = whichPerturb[param.perturbationMethod]
    catch
        @warn("Invalid parameter for perturbationMethod, default used")
    end
    try
        restartmethod! = whichRestart[param.restartMethod]
    catch
        @warn("Invalid parameter for perturbationMethod, default used")
    end
    restart = param.restart
    pasRestart = 0
    checkCycle = param.checkCycle
    presolve = param.presolve
    improve = param.improve
    total_time_limit = param.timelim
    feasibilityCheck = param.feasibilityCheck

    getData = param.logLevel

    # For getData
    myData = Vector{Float64}([
        0.0, # Percentage of variables fixed with presolve      1
        0.0, # Time for presolve                                2
        0.0, # Time for first linear relaxation                 3
        0.0, # Average time for projection                      4
        0.0, # Number of perturbation                           5
        0.0, # Average time for perturbation                    6
        0.0, # Average time for rounding                        7
        0.0, # Average time for scatter search                  8
        0.0, # Average time for feasibility check               9
        0.0, # Number of iterations                            10
        0.0, # Time for improve                                11
        0.0, # Number of feasibility checks                    12
    ])

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
    for i in 1:(checkCycle-1)
        push!(saveSol, zeros(Float64, nb_variables))
    end
    saveVect = 0

    initial_objective_function = CPLEX.get_obj(model)
    objective_function = copy(initial_objective_function)
    objSense_init = CPLEX.get_sense(model)
    if objSense_init == :Max
        objective_function = -objective_function
    end

    # This is only useful for recursiverounding!
    sortBinVariables = false
    sortBinVariablesOnce = false
    objectiveFunctionRR = zeros(Float64, nb_variables)
    if param.recursiveRounding == 1
        sortBinVariables = true
        sortBinVariablesOnce = true
    elseif param.recursiveRounding == 2
        sortBinVariables = true
        sortBinVariablesOnce = false
    elseif param.recursiveRounding == 3
        sortBinVariables = true
        sortBinVariablesOnce = true
        objectiveFunctionRR .= initial_objective_function
    elseif param.recursiveRounding == 4
        sortBinVariables = true
        sortBinVariablesOnce = false
        objectiveFunctionRR .= initial_objective_function
    elseif param.recursiveRounding == 0
        #Nothing
    else
        @warn "Parameter recursiveRounding has an unauthorized value, parameter ignored"
    end

    #=--------------------------------------------------------------------------
    Variables needed for presolve and to check feasibility
    --------------------------------------------------------------------------=#
    A = CPLEX.get_constr_matrix(model)
    senses = CPLEX.get_constr_senses(model)
    rhs = CPLEX.get_rhs(model)
    nConstr = length(rhs)

    if presolve != 0
        timeTmp = @timed begin
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
            if getData
                myData[1] = (nBin-length(indices))*100/nBin
            end
        end
        if getData
            myData[2] = timeTmp[2]
        end
    end

    indicesRR = copy(indices)
    A = sparse(A') # This transpose makes feasibility check much faster

    if param.TT
        TTmin = param.TTmin
        TTmax = param.TTmax
    else
        TTmin = round(Int, length(indices)*param.TTmin/100, RoundDown)
        TTmax = round(Int, length(indices)*param.TTmax/100, RoundUp)
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
        TTmax = param.TTmin
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
    timeLim_init = CPLEX.get_param(model.env, "CPX_PARAM_TILIM")
    timeTmp = @timed begin
        CPLEX.set_param!(model.env, "CPX_PARAM_TILIM", maximum([0, total_time_limit - (time_ns() - timeStart) / 1.0e9]))
        CPLEX.optimize!(model)

        xOverline = zeros(Float64, nb_variables)
        if CPLEX.get_status(model) == :CPX_STAT_OPTIMAL
            xOverline = CPLEX.get_solution(model)
        end
    end
    if getData
        myData[3] = timeTmp[2]
    end

    isFeasible = false
    if feasibilityCheck != 2
        isFeasible = isfeasible_xOverline(xOverline, indices, timeStart, total_time_limit)
    end
    if feasibilityCheck != 2
        if getData
            myData[9] += timeTmp[2]
            myData[12] += 1.0
        end
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

        timeTmp = @timed begin
            roundingmethod!(xTilde, xOverline, timeStart, total_time_limit, indices, model, indicesRR, objectiveFunctionRR)
        end
        xTildeBis = copy(xTilde)
        for i in indices
            if xTilde[i] > 0.5
                freqRound[i] += 1
            end
        end
        if getData
            myData[7] += timeTmp[2]
        end

        timeTmp = @timed if feasibilityCheck != 1
            isFeasible = isfeasible_xTilde(xTilde, A, senses, rhs, timeStart, total_time_limit)
        end
        if feasibilityCheck != 1
            if getData
                myData[9] += timeTmp[2]
                myData[12] += 1.0
            end
        end
    else
        xTilde .= xOverline
    end

    #=--------------------------------------------------------------------------
    The feasibility pump
    --------------------------------------------------------------------------=#
    while !isFeasible && (time_ns()-timeStart)/1.0e9 < total_time_limit
        nIter += 1
        pasRestart += 1
        alpha *= alphaChange
        timeTmp = @timed begin
            xOverline = projectionmethod(model, timeStart, total_time_limit, xTilde, indices, objective_function, alpha, coef)
        end
        if getData
            myData[4] += timeTmp[2]
        end
        timeTmp = @timed if feasibilityCheck != 2
            isFeasible = isfeasible_xOverline(xOverline, indices, timeStart, total_time_limit)
        end
        if feasibilityCheck != 2
            if getData
                myData[9] += timeTmp[2]
                myData[12] += 1.0
            end
        end

        if !isFeasible
            if sortBinVariables
                sortReducedCosts!(indicesRR, model)
            end

            timeTmp = @timed begin
                roundingmethod!(xTildeBis, xOverline, timeStart, total_time_limit, indices, model, indicesRR, objectiveFunctionRR)
            end
            if getData
                myData[7] += timeTmp[2]
            end
            if pasRestart == restart
                pasRestart = 0
                restartmethod!(xTildeBis, xOverline, indices)
            end
            if xTildeBis == xTilde
                timeTmp = @timed begin
                    perturbmethod!(xTildeBis, xOverline, indices, TTmin, TTmax, freqRound, nIter)
                end
                if getData
                    myData[5] += 1
                    myData[6] += timeTmp[2]
                end
            else
                if xTildeBis in saveSol
                    pasRestart = 0
                    restartmethod!(xTildeBis, xOverline, indices)
                end
            end
            if checkCycle > 1
                saveSol[saveVect+1] .= xTildeBis
                saveVect = mod(saveVect+1, checkCycle-1)
            end
            xTilde = copy(xTildeBis)

            for i in indices
                if xTilde[i] > 0.5
                    freqRound[i] += 1
                end
            end
            timeTmp = @timed if feasibilityCheck != 1
                isFeasible = isfeasible_xTilde(xTilde, A, senses, rhs, timeStart, total_time_limit)
            end
            if feasibilityCheck != 1
                if getData
                    myData[9] += timeTmp[2]
                    myData[12] += 1.0
                end
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
        timeTmp = @timed if improve
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

            CPLEX.set_param!(model.env, "CPX_PARAM_TILIM", maximum([0, total_time_limit - (time_ns() - timeStart) / 1.0e9]))
            CPLEX.optimize!(model)
            if CPLEX.get_status(model) == :CPX_STAT_OPTIMAL
                xTilde = CPLEX.get_solution(model)
                objVal = CPLEX.get_objval(model)
            end
            if getData
                myData[11] = timeTmp[2]
            end
        end
    else
        xTilde = repeat([NaN], nb_variables)
    end

    if getData
        myData[4] /= nIter
        myData[6] /= myData[5]
        myData[7] /= nIter+1
        myData[8] /= nIter-1
        myData[9] /= myData[12]
        myData[10] = nIter

        println(myData)
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
    CPLEX.set_param!(model.env, "CPX_PARAM_TILIM", timeLim_init)

    return status, xTilde, objVal
end

feasibilitypump(model::CPLEX.Model) =
    feasibilitypump(model, ParametersFP())



function isfeasible_xTilde(
        xTilde::Vector{Float64},
        A::SparseArrays.SparseMatrixCSC{Float64,Int64},
        senses::Vector{Int8},
        rhs::Vector{Float64},
        timeStart::UInt,
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
        if (time_ns()-timeStart)/1.0e9 > total_time_limit
            isFeasible = false
        end
    end

    return isFeasible
end

function isfeasible_xOverline(
        xOverline::Vector{Float64},
        indices::Vector{Int},
        timeStart::UInt,
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
        if (time_ns()-timeStart)/1.0e9 > total_time_limit
            isFeasible = false
        end
        ind += 1
    end

    return isFeasible
end


function initializeparametersfp()
    return ParametersFP()
end


function setparamfp!(paramList::ParametersFP, param::Int, value::Int)
    if param == 10
        paramList.roundingMethod = value
    elseif param == 20
        paramList.projectionMethod = value
    elseif param == 30
        paramList.perturbationMethod = value
    elseif param == 32
        paramList.TTmin = value
    elseif param == 33
        paramList.TTmax = value
    elseif param == 34
        paramList.checkCycle = value
    elseif param == 35
        paramList.restart = value
    elseif param == 2
        paramList.timelim = Float64(value)
    elseif param == 11
        paramList.recursiveRounding = value
    elseif param == 41
        paramList.presolve = value
    elseif param == 100
        paramList.feasibilityCheck = value
    elseif param == 1
        try
            paramList.logLevel = Bool(value)
        catch
            @warn("The parameter "*string(param)*" can take only 0 or 1")
        end
    elseif param == 31
        try
            paramList.TT = Bool(value)
        catch
            @warn("The parameter "*string(param)*" can take only 0 or 1")
        end
    elseif param == 42
        try
            paramList.improve = Bool(value)
        catch
            @warn("The parameter "*string(param)*" can take only 0 or 1")
        end
    else
        @warn("The parameter "*string(param)*" doesn't exist")
    end
    return paramList
end

function setparamfp!(paramList::ParametersFP, param::Int, value::Float64)
    if param == 10
        paramList.roundingMethod = round(Int, value)
    elseif param == 20
        paramList.projectionMethod = round(Int, value)
    elseif param == 30
        paramList.perturbationMethod = round(Int, value)
    elseif param == 32
        paramList.TTmin = round(Int, value)
    elseif param == 33
        paramList.TTmax = round(Int, value)
    elseif param == 34
        paramList.checkCycle = round(Int, value)
    elseif param == 35
        paramList.restart = round(Int, value)
    elseif param == 2
        paramList.timelim = value
    elseif param == 11
        paramList.recursiveRounding = round(Int, value)
    elseif param == 41
        paramList.presolve = round(Int, value)
    elseif param == 100
        paramList.feasibilityCheck = round(Int, value)
    elseif param == 1
        try
            paramList.logLevel = Bool(round(Int, value))
        catch
            @warn("The parameter "*string(param)*" can take only 0 or 1")
        end
    elseif param == 31
        try
            paramList.TT = Bool(round(Int, value))
        catch
            @warn("The parameter "*string(param)*" can take only 0 or 1")
        end
    elseif param == 42
        try
            paramList.improve = Bool(round(Int, value))
        catch
            @warn("The parameter "*string(param)*" can take only 0 or 1")
        end
    else
        @warn("The parameter "*string(param)*" doesn't exist")
    end
    return paramList
end


function setparamfp!(paramList::ParametersFP, param::String, value::Union{Float64, Int})
    if haskey(__paramToInt, param)
        setparamfp!(paramList, __paramToInt[param], value)
    else
        @warn("The parameter "*string(param)*" doesn't exist")
    end
    return paramList
end


function getparamfp(paramList::ParametersFP, param::Int)
    value = NaN
    if param == 10
        value = paramList.roundingMethod
    elseif param == 20
        value = paramList.projectionMethod
    elseif param == 30
        value = paramList.perturbationMethod
    elseif param == 31
        value = paramList.TT
    elseif param == 32
        value = paramList.TTmin
    elseif param == 33
        value = paramList.TTmax
    elseif param == 34
        value = paramList.checkCycle
    elseif param == 35
        value = paramList.restart
    elseif param == 1
        value = paramList.logLevel
    elseif param == 2
        value = paramList.timelim
    elseif param == 11
        value = paramList.recursiveRounding
    elseif param == 41
        value = paramList.presolve
    elseif param == 100
        value = paramList.feasibilityCheck
    elseif param == 42
        value = paramList.improve
    else
        @warn("The parameter "*string(param)*" doesn't exist")
    end
    return value
end


function getparamfp(paramList::ParametersFP, param::String)
    value = NaN
    if haskey(__paramToInt, param)
        value = getparamfp(paramList, __paramToInt[param])
    end
    return value
end

function getparamnamefp(param::Int)
    paramName = ""
    if haskey(__paramToString, param)
        paramName = __paramToString[param]
    else
        @warn("No parameter "*string(param))
    end
    return paramName
end
