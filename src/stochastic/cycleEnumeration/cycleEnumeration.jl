using Cbc
using DelimitedFiles
using GLPK
using Graphs
using Graphs.Experimental
using JuMP
using Random
using Requires
using TimerOutputs


include("../../instance/instance.jl")
include("../../branch&price/typedef.jl")
include("../../branch&price/branch_and_price.jl")
include("../../branch&price/master.jl")
include("../../branch&price/subproblem.jl")
include("../../branch&price/node.jl")
include("../../branch&price/typedef.jl")
include("../../createModel.jl")
include("cycle.jl")



function solve_with_CE(
        filename::String,
        K::Int,
        L::Int,
        bp_params::BP_params = BP_params(),
        timer::TimerOutput = TimerOutput(),
        time_limit::Float64 = 600.0
)
    start_time = time()
    reset_timer!(timer)
    Random.seed!(10)
    if bp_params.verbose
        printstyled("\n********************************************************************************\n Solve $filename with (K,L) = ($K,$L) using branch-and-price and column enumeration\n - time limit is $time_limit seconds\n********************************************************************************\n\n" ; bold = true, color = :magenta)
    end

    # Parsing
    if bp_params.verbose
    printstyled("\n----------------------------------------------------------\n Parse the input file\n----------------------------------------------------------\n\n" ; color = :yellow) end
    instance = @timeit timer "Parser" Instance_stochastic(filename, K, L)

    # Preprocessing
    if bp_params.verbose printstyled("\n----------------------------------------------------------\n Preprocessing: enumerate columns and compute their expected value\n----------------------------------------------------------\n\n" ; color = :yellow) end
    column_pool = @timeit timer "Preprocessing" column_enumeration(instance)

    # Call the branch-and-price algorithm
    if bp_params.verbose printstyled("\n----------------------------------------------------------\n Solve with branch-and-price\n----------------------------------------------------------\n\n" ; color = :yellow) end
    bp_status = @timeit timer "B&P" branch_and_price(instance, subgraphs, bp_params, timer, time_limit - (time() - start_time))
    bp_status.solve_time = TimerOutputs.time(timer["B&P"])/10^9

    # print the number of cycles and chains for each column's length
    print_and_check_solution(bp_status.best_cycles, bp_status.best_chains, instance, bp_params.verbose)

    # Print cpu profiling
    if bp_params.verbose println(timer) end

    return bp_status, Graph_info(instance), Subgraph_info(subgraphs);
    instance = @timeit timer "Parser" Instance(filename, K, L)

    bp_status = @timeit timer "B&P" ce_branch_and_price(instance, cycles, nb_cycles, bp_params, timer, time_limit - (time() - start_time))

    println(timer)

    return bp_status
    
end

