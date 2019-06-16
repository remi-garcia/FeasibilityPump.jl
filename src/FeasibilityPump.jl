__precompile__()

#===============================================================================
 This package contains the metaheuristic the feasibility pump using CPLEX
===============================================================================#

module FeasibilityPump

using CPLEX
using LinearAlgebra
using Random
using SparseArrays

include("solver.jl")
include("rounding.jl")
include("projection.jl")
include("perturbation.jl")
include("constants.jl")

#export ParametersFP
export feasibilitypump
export setparamfp!
export getparamnamefp
export getparamfp
export initializeparametersfp
export __compilefp__

end # module
