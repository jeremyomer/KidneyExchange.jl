include("hpief.jl")
include("reduced_extended_edge.jl")
include("relaxed_arc.jl")

"""
    solve_with_mip

This is the main function to call to solve the input instance with given bounds on the length of covering cycles and chains and given options

#Input parameters
* `filename::String`: path of the input data files, this should include the name of the files, but not the .dat and .wmd extensions
* `K::Int`: maximum length of cycles
* `L::Int`: Maximum length of chains
* `params::MIP_params`: parameters of the MIP model and solver
* `timer::TimerOutput`: a timer that will provide detail of where the computational time was spent during the branch-and-price
* `time_limit::Float64`: time limit of the algorithm, including parsing and prepreprocessing times

#Output parametes
* `instance::Instance`: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains.
* `subgraphs::Graph_copies`: Description of the graph copies of the extended edge formulation
* `bp_status::BP_status`:  Structure containing every relevant information on the execution of the algorithm (including the optimal solution)
"""
function solve_with_mip(filename::String, K::Int, L::Int, params::MIP_params = MIP_params(), timer::TimerOutput = TimerOutput(), time_limit::Float64 = 600.0)
    start_time = time()
    reset_timer!(timer)
    Random.seed!(10)
    if params.verbose
        printstyled("\n********************************************************************************\n Solve $filename with (K,L) = ($K,$L) using a compact model\n - model type is: $(params.model_type) \n - time limit is $time_limit seconds\n********************************************************************************\n\n" ; bold = true, color = :magenta)
    end

    # Parsing
    if params.verbose
        printstyled("\n----------------------------------------------------------\n Parse the input file\n----------------------------------------------------------\n\n" ; color = :yellow)
    end
    instance = @timeit timer "Parser" Instance(filename, K, L)

    # Preprocessing
    if params.verbose
        printstyled("\n----------------------------------------------------------\n Preprocessing: compute the graph copies\n----------------------------------------------------------\n\n" ; color = :yellow)
    end
    subgraphs = @timeit timer "Preprocessing" preprocess_graph_copies(instance, true, params.reduce_vertices, params.fvs)

    # Call the branch-and-price algorithm
    if params.verbose
        printstyled("\n----------------------------------------------------------\n Solve the compact MIP\n----------------------------------------------------------\n\n" ; color = :yellow)
    end

    val = -1.0
    if params.model_type == HPIEF
        model, κ = @timeit timer "Build MIP" build_hpief_mip(instance, subgraphs, params, time_limit)
        mip_status = @timeit timer "Solve MIP" solve_hpief_mip(model, params, instance, κ)
    elseif params.model_type == EXTENDED_EDGE
        model = @timeit timer "Build MIP" build_reduced_extended_edge_mip(instance, subgraphs, params, time_limit)
        mip_status = @timeit timer "Solve MIP" solve_reduced_extended_edge_mip(model, params, instance, subgraphs)
    elseif params.model_type == RELAXED_ARC
        if L != 0
            error("The relaxed arc formulation is valid only for L=0")
        end
        mip_status = @timeit timer "Solve MIP" relaxed_arc(instance, params, time_limit)
    end

    # print the number of cycles and chains for each column's length
    print_and_check_solution(mip_status.best_cycles, mip_status.best_chains, instance, params.verbose)

    # Print cpu profiling
    if (params.verbose) println(timer) end

    return mip_status, Graph_info(instance), Subgraph_info(subgraphs);
end
