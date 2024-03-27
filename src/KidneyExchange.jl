module KidneyExchange

# Package dependencies
using HiGHS
using DelimitedFiles
using GLPK
using Graphs
using JuMP
using Random
using Requires
using TimerOutputs
using DocStringExtensions

export print_and_check_solution
include("utils.jl")

# data generation functions
export generate_saidman_instance
export generate_sparse_unos_instance
export generate_heterogeneous_instance
export generate_complete_benchmark
include("instance/create_benchmark.jl")
include("instance/matching_instance.jl")

export Mip_model, HPIEF, EXTENDED_EDGE, RELAXED_ARC
export MIP_params
export solve_with_mip
include("models/mip.jl")


export BP_params
export solve_with_BP
export TimerOutput
const ϵ = 1.0e-3
include("branch&price/branch_and_price.jl")

function __init__()
    @require Gurobi = "2e9cd046-0924-5485-92f1-d5272153d98b" @eval using .Gurobi
    @require CPLEX = "a076750e-1247-5638-91d2-ce28b192dca0" @eval using .CPLEX
    @require Clp = "2e9cd046-0924-5485-92f1-d5272153d98b" @eval using .Clp
    @require Cbc = "9961bab8-2fa3-5c5a-9d89-47fab24efd76" @eval using .Cbc
end

end
