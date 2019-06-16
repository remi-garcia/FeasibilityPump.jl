# Constants
const whichRounding = Dict{Int, Function}(
    0 => simplerounding!,
    1 => recursiverounding!,
)

const whichProjection = Dict{Int, Function}(
    0 => minimumdistance,
    1 => minimumdeltaalpha,
    2 => minimumdistance_reduce,
)

const whichPerturb = Dict{Int, Function}(
    -1 => noperturb,
    0 => perturb!,
    1 => randomperturb!,
    2 => frequencyperturb!,

    100 => perturbRestart!,
)


const whichRestart = Dict{Int, Function}(
    -1 => norestart,
    0 => perturbRestart!,

    100 => perturb!,
    101 => randomperturb!,
    102 => frequencyperturb!,
)


const __paramToInt = Dict{String, Int}([
    "logLevel" => 1,
    "timelim" => 2,
    "roundingMethod" => 10,
    "recursiveRounding" => 11,
    "projectionMethod" => 20,
    "perturbationMethod" => 30,
    "perturbationTT" => 31,
    "perturbationTTmin" => 32,
    "perturbationTTmax" => 33,
    "checkCycle" => 34,
    "restart" => 35,
    "presolve" => 41,
    "improve" => 42,
    "feasibilityCheck" => 100
])


const __paramToString = Dict{Int, String}([
    1 => "logLevel",
    2 => "timelim",
    10 => "roundingMethod",
    11 => "recursiveRounding",
    20 => "projectionMethod",
    30 => "perturbationMethod",
    31 => "perturbationTT",
    32 => "perturbationTTmin",
    33 => "perturbationTTmax",
    34 => "checkCycle",
    35 => "restart",
    41 => "presolve",
    42 => "improve",
    100 => "feasibilityCheck"
])
