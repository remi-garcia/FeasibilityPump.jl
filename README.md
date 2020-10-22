# FeasibilityPump

**Work in progress.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

FeasibilityPump is a Julia package that solve a linear program with binary variables. It uses the Feasibility Pump heuristic or a variant.

This package exports a few functions: `feasibilitypump()`, `setparamfp!()`, `getparamnamefp()`, `getparamfp()`, `initializeparametersfp()` and `__compilefp__()`

### Installing

```julia
julia> ] add https://github.com/remi-garcia/FeasibilityPump.jl
julia> using FeasibilityPump
```

## Running the first test

```julia
julia> include("examples/main.jl")
```
It should run over the `.mps` in `/data` using CPLEX.

### Usage

```julia
julia> model = read_from_file("filepath.mps")
julia> parameters = initialize_parameters_fp()
julia> #set_parameters_fp!(parameters, "parameter_name", parameter_value) # See docs for more info on this
julia> status, solution, objective_value = feasibility_pump(model, parameters)
```

### Problem Collections

* MIPLIB: https://miplib.zib.de/

### References

* [The feasibility pump](https://doi.org/10.1007/s10107-004-0570-3), Fischetti, M., Glover, F. & Lodi. A. Math. Program. 104:91. 2005.
* [Improving the feasibility pump](https://doi.org/10.1016/j.disopt.2006.10.004), Achterberg, T. & Berthold, T. Discrete Optimization 4(1) 77:86. 2007.
* TODO
