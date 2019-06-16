__precompile__()

using CPLEX
include("../src/FeasibilityPump.jl")
using .FeasibilityPump
__compilefp__()
using Random


function mainBasic()
    files = "../data/" .* readdir("../data/")

    nSeed = 3

    for file in files
        fileName = split(file, ['/'])[end][1:end-4]

        println("--- "*file*" ---")
        println("- Solving with feasibility pump:")
        for i in 1:nSeed
            GC.gc()
            Random.seed!(i)
            myModel = CPLEX.Model(CPLEX.Env())
            CPLEX.read_model(myModel, file)
            CPLEX.set_param!(myModel.env, "CPX_PARAM_THREADS", 1)

            param = initializeparametersfp()
            #setparamfp!(param, "logLevel", 1)

            timeTmp = @timed begin
                st, x, z = feasibilitypump(myModel, param)
            end
            print("-Time: ")
            println(timeTmp[2])
            println("-Status: ", st)
            println("-ObjectiveValue: ", z)
            println()
            println()
        end
    end
end

mainBasic()
