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
        myModel::CPLEX.Model,
        param::ParametersFP
        ;
        compile::Bool = false
    )
    #=---------------------------------
     Allow to precompile this function
    ---------------------------------=#
    if compile
        return :COMPILE, [NaN], NaN
    end

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
    timeLim = param.timelim
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

    LPMethod_init = CPLEX.get_param(myModel.env, "CPX_PARAM_LPMETHOD")
    CPLEX.set_param!(myModel.env, "CPX_PARAM_LPMETHOD", 1)
    nVars = CPLEX.num_var(myModel)
    varTypes_init = CPLEX.get_vartype(myModel)
    varLB_init = CPLEX.get_varLB(myModel)
    varLB = copy(varLB_init)
    varUB_init = CPLEX.get_varUB(myModel)
    varUB = copy(varUB_init)
    indices = Vector{Int}()
    for i in 1:nVars
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
        push!(saveSol, zeros(Float64, nVars))
    end
    saveVect = 0

    objectiveFunction_init = CPLEX.get_obj(myModel)
    objectiveFunction = copy(objectiveFunction_init)
    objSense_init = CPLEX.get_sense(myModel)
    if objSense_init == :Max
        objectiveFunction = -objectiveFunction
    end

    # This is only useful for recursiverounding!
    sortBinVariables = false
    sortBinVariablesOnce = false
    objectiveFunctionRR = zeros(Float64, nVars)
    if param.recursiveRounding == 1
        sortBinVariables = true
        sortBinVariablesOnce = true
    elseif param.recursiveRounding == 2
        sortBinVariables = true
        sortBinVariablesOnce = false
    elseif param.recursiveRounding == 3
        sortBinVariables = true
        sortBinVariablesOnce = true
        objectiveFunctionRR .= objectiveFunction_init
    elseif param.recursiveRounding == 4
        sortBinVariables = true
        sortBinVariablesOnce = false
        objectiveFunctionRR .= objectiveFunction_init
    elseif param.recursiveRounding == 0
        #Nothing
    else
        @warn "Parameter recursiveRounding has an unauthorized value, parameter ignored"
    end

    #=--------------------------------------------------------------------------
    Variables needed for presolve and to check feasibility
    --------------------------------------------------------------------------=#
    A = CPLEX.get_constr_matrix(myModel)
    senses = CPLEX.get_constr_senses(myModel)
    rhs = CPLEX.get_rhs(myModel)
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
                    if presolve == 1 || (val == 0 && objectiveFunction[i] >= 0) || (val == 1 && objectiveFunction[i] <= 0)
                        varLB[i] = val
                        varUB[i] = val
                        indicesTmp[ind] = true
                    end
                end
            end
            CPLEX.set_varLB!(myModel, varLB)
            CPLEX.set_varUB!(myModel, varUB)
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

    problemType_init = CPLEX.get_prob_type(myModel)
    CPLEX.set_prob_type!(myModel, 0)
    has_int_init = myModel.has_int
    myModel.has_int = false

    #=--------------------------------------------------------------------------
    We initialize variables for variants:
      coef, alpha and alphaChange are parameters of objective
      feasibility pump. freqRound is useful for frequenciesrounding
    --------------------------------------------------------------------------=#
    coef = norm(objectiveFunction)
    if coef > 1e-8
        coef = sqrt(length(indices))/coef
    end
    alphaChange = 0.9
    alpha = 1.0

    freqRound = zeros(Int, nVars)

    #=--------------------------------------------------------------------------
    First linear relaxation, rounding and feasibility check
    --------------------------------------------------------------------------=#
    #CPLEX.set_obj!(myModel, zeros(Float64, nVars))
    timeLim_init = CPLEX.get_param(myModel.env, "CPX_PARAM_TILIM")
    timeTmp = @timed begin
        CPLEX.set_param!(myModel.env, "CPX_PARAM_TILIM", maximum([0, timeLim - (time_ns() - timeStart) / 1.0e9]))
        CPLEX.optimize!(myModel)

        xOverline = zeros(Float64, nVars)
        if CPLEX.get_status(myModel) == :CPX_STAT_OPTIMAL
            xOverline = CPLEX.get_solution(myModel)
        end
    end
    if getData
        myData[3] = timeTmp[2]
    end

    isFeasible = false
    if feasibilityCheck != 2
        isFeasible = isfeasible_xOverline(xOverline, indices, timeStart, timeLim)
    end
    if feasibilityCheck != 2
        if getData
            myData[9] += timeTmp[2]
            myData[12] += 1.0
        end
    end
    nIter = 0
    xTilde = zeros(Float64, nVars)

    if !isFeasible
        if sortBinVariables
            sortReducedCosts!(indicesRR, myModel)
            if sortBinVariablesOnce
                sortBinVariables = false
            end
        end

        CPLEX.set_sense!(myModel, :Min)

        timeTmp = @timed begin
            roundingmethod!(xTilde, xOverline, timeStart, timeLim, indices, myModel, indicesRR, objectiveFunctionRR)
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
            isFeasible = isfeasible_xTilde(xTilde, A, senses, rhs, timeStart, timeLim)
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
    while !isFeasible && (time_ns()-timeStart)/1.0e9 < timeLim
        nIter += 1
        pasRestart += 1
        alpha *= alphaChange
        timeTmp = @timed begin
            xOverline = projectionmethod(myModel, timeStart, timeLim, xTilde, indices, objectiveFunction, alpha, coef)
        end
        if getData
            myData[4] += timeTmp[2]
        end
        timeTmp = @timed if feasibilityCheck != 2
            isFeasible = isfeasible_xOverline(xOverline, indices, timeStart, timeLim)
        end
        if feasibilityCheck != 2
            if getData
                myData[9] += timeTmp[2]
                myData[12] += 1.0
            end
        end

        if !isFeasible
            if sortBinVariables
                sortReducedCosts!(indicesRR, myModel)
            end

            timeTmp = @timed begin
                roundingmethod!(xTildeBis, xOverline, timeStart, timeLim, indices, myModel, indicesRR, objectiveFunctionRR)
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
                isFeasible = isfeasible_xTilde(xTilde, A, senses, rhs, timeStart, timeLim)
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
        objVal = sum(objectiveFunction_init .* xTilde)
        timeTmp = @timed if improve
            CPLEX.set_sense!(myModel, objSense_init)
            CPLEX.set_obj!(myModel, objectiveFunction_init)
            varLB = CPLEX.get_varLB(myModel)
            varUB = CPLEX.get_varUB(myModel)
            for i in indices
                varLB[i] = round(xTilde[i])
                varUB[i] = varLB[i]
            end
            CPLEX.set_varLB!(myModel, varLB)
            CPLEX.set_varUB!(myModel, varUB)

            CPLEX.set_param!(myModel.env, "CPX_PARAM_TILIM", maximum([0, timeLim - (time_ns() - timeStart) / 1.0e9]))
            CPLEX.optimize!(myModel)
            if CPLEX.get_status(myModel) == :CPX_STAT_OPTIMAL
                xTilde = CPLEX.get_solution(myModel)
                objVal = CPLEX.get_objval(myModel)
            end
            if getData
                myData[11] = timeTmp[2]
            end
        end
    else
        xTilde = repeat([NaN], nVars)
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
    CPLEX.set_param!(myModel.env, "CPX_PARAM_LPMETHOD", LPMethod_init)
    CPLEX.set_vartype!(myModel, varTypes_init)
    CPLEX.set_vartype!(myModel, varTypes_init)
    CPLEX.set_varLB!(myModel, varLB_init)
    CPLEX.set_varUB!(myModel, varUB_init)
    CPLEX.set_obj!(myModel, objectiveFunction_init)
    CPLEX.set_prob_type!(myModel, problemType_init)
    CPLEX.set_sense!(myModel, objSense_init)
    myModel.has_int = has_int_init
    CPLEX.set_param!(myModel.env, "CPX_PARAM_TILIM", timeLim_init)

    return status, xTilde, objVal
end

feasibilitypump(myModel::CPLEX.Model; compile::Bool = false) =
    feasibilitypump(myModel, ParametersFP(), compile = compile)



function isfeasible_xTilde(
        xTilde::Vector{Float64},
        A::SparseArrays.SparseMatrixCSC{Float64,Int64},
        senses::Vector{Int8},
        rhs::Vector{Float64},
        timeStart::UInt,
        timeLim::Float64
        ;
        compile = false
    )
    if compile
        return false
    end

    nVars::Int64, nConstr::Int64 = size(A)
    isFeasible = true
    sum::Float64 = 0.0
    i = 1
    while isFeasible && i <= nConstr
    #for i in 1:nConstr
        sum = 0.0
        #for j in 1:nVars
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
        if (time_ns()-timeStart)/1.0e9 > timeLim
            isFeasible = false
        end
    end

    return isFeasible
end

function isfeasible_xOverline(
        xOverline::Vector{Float64},
        indices::Vector{Int},
        timeStart::UInt,
        timeLim::Float64
        ;
        compile = false
    )
    if compile
        return false
    end

    nBin = length(indices)
    ind = 1
    isFeasible = true
    while ind <= nBin && isFeasible
        i = indices[ind]
        if abs(xOverline[i]-round(xOverline[i])) > 1e-8
            isFeasible = false
        end
        if (time_ns()-timeStart)/1.0e9 > timeLim
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



function __compilefp__()
    feasibilitypump(
        CPLEX.Model(CPLEX.Env()),
        compile = true
    )

    # --- hybridation.jl --- #

    # --- perturbation.jl --- #
    for perturbmethod! in values(whichPerturb)
        perturbmethod!(
            [NaN],
            [NaN],
            [0],
            0,
            0,
            [0],
            0,
            compile = true
        )
    end

    # --- projection.jl --- #
    for projectionmethod in values(whichProjection)
        projectionmethod(
            CPLEX.Model(CPLEX.Env()),
            UInt(0),
            NaN,
            [NaN],
            [0],
            [NaN],
            NaN,
            NaN,
            compile = true
        )
    end

    # --- rounding.jl --- #
    for roundingmethod! in values(whichRounding)
        roundingmethod!(
            [NaN],
            [NaN],
            UInt(0),
            NaN,
            [0],
            CPLEX.Model(CPLEX.Env()),
            [0],
            [NaN],
            compile = true
        )
    end

    sortReducedCosts!(
        [0],
        CPLEX.Model(CPLEX.Env()),
        compile = true
    )

    isfeasible_xTilde(
        [NaN],
        sparse([NaN NaN; NaN NaN]),
        Vector{Int8}([0]),
        [NaN],
        UInt(0),
        NaN,
        compile = true
    )
    isfeasible_xOverline(
        [NaN],
        [0],
        UInt(0),
        NaN,
        compile = true
    )

    return 0
end
