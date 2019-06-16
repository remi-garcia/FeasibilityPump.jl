# FeasibilityPump

**Work in progress.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

FeasibilityPump is a Julia package that solve a linear program with binary variables. It uses the Feasibility Pump heuristic or a variant.

This package exports 2 functions: `readmps()` and `mpstomodel()`

### Prerequisites

julia 1.0.0</br>
CPLEX 0.5.0</br>
LinearAlgebra</br>
SparseArrays</br>
Random</br>

### Installing

TODO

## Running the first test

```julia
julia> include("examples/main.jl")
```
It should run over the .mps in data/

### Usage

```julia
julia> myModel = CPLEX.Model(CPLEX.Env())
julia> CPLEX.read_model(myModel, "filepath.mps")
julia> param = initializeparametersfp()
julia> #setparamfp!(param, "paramName", paramValue) # See docs for more info on this
julia> status, solution, objectiveValue = feasibilitypump(myModel, param)
```

### Problem Collections

* MIPLIB: https://miplib.zib.de/

### References

* [The feasibility pump](https://doi.org/10.1007/s10107-004-0570-3), Fischetti, M., Glover, F. & Lodi. A. Math. Program. 104:91. 2005.
* [Improving the feasibility pump](https://doi.org/10.1016/j.disopt.2006.10.004), Achterberg, T. & Berthold, T. Discrete Optimization 4(1) 77:86. 2007.
* TODO
