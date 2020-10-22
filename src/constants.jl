# Constants
const which_rounding = Dict{Int, Function}(
    0 => simple_rounding!,
    1 => recursive_rounding!,
)

const which_projection = Dict{Int, Function}(
    0 => minimum_distance,
    1 => minimum_delta_alpha,
    2 => minimum_distance_reduce,
)

const which_perturbation = Dict{Int, Function}(
    -1 => no_perturbation,
    0 => perturbation!,
    1 => random_perturbation!,
    2 => frequency_perturbation!,

    100 => perturbation_restart!,
)


const which_restart = Dict{Int, Function}(
    -1 => no_restart,
    0 => perturbation_restart!,

    100 => perturbation!,
    101 => random_perturbation!,
    102 => frequency_perturbation!,
)


const __parameters_to_int = Dict{String, Int}([
    "rounding_method" => 10,
    "recursive_rounding" => 11,
    "projection_method" => 20,
    "perturbation_method" => 30,
    "perturbation_TT" => 31,
    "perturbation_TTmin" => 32,
    "perturbation_TTmax" => 33,
    "check_cycle" => 34,
    "restart" => 35,
    "presolve" => 41,
    "improve" => 42,
    "feasibility_check" => 100
])


const __parameters_to_string = Dict{Int, String}([
    10 => "rounding_method",
    11 => "recursive_rounding",
    20 => "projection_method",
    30 => "perturbation_method",
    31 => "perturbation_TT",
    32 => "perturbation_TTmin",
    33 => "perturbation_TTmax",
    34 => "check_cycle",
    35 => "restart",
    41 => "presolve",
    42 => "improve",
    100 => "feasibility_check"
])
