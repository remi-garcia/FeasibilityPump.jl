#===============================================================================
 This package contains the metaheuristic the feasibility pump
===============================================================================#

module FeasibilityPump

using JuMP
using LinearAlgebra
using Random
using SparseArrays

include("solver.jl")
include("rounding.jl")
include("projection.jl")
include("perturbation.jl")
include("constants.jl")

#export ParametersFP
export feasibility_pump
export set_parameters_fp!
export get_parameter_names_fp
export get_parameters_fp
export initialize_parameters_fp

end # module
