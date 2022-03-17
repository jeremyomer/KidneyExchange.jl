using Cbc
using DelimitedFiles
using GLPK
using Graphs
using JuMP
using Random
using Requires
using TimerOutputs

include("../instance/instance.jl")
include("typedef.jl")

function sto_solve(filename::AbstractString, K::Int, L::Int, bp_params::BP_params, timer::TimerOutput = TimerOutput(), time_limit::Float64 = 600.0)
	sto_solve_bp(filename, K, L, bp_params, timer, time_limit)
end

function sto_solve_bp(filename::AbstractString, K::Int, L::Int, bp_params::BP_params, timer::TimerOutput = TimerOutput(), time_limit::Float64 = 600.0)
	#<!> Add verbose
	start_time = time()
    reset_timer!(timer)
    Random.seed!(10)

    instance = @timeit timer "Parser" Instance(filename, K, L) #<!> Change Instance to reflect the stochastic model

    p_status = @timeit timer "B&P" sto_branch_and_price(instance, bp_params, timer, time_limit - (time() - start_time))
    bp_status.solve_time = TimerOutputs.time(timer["B&P"])/10^9

    println(timer)

    return bp_status

end

function sto_branch_and_price(instance::Instance, bp_params::BP_params, timer::TimerOutput, time_limit::Float64)
	graph = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length
    B = instance.crossmatch_budget
    start_time=time()
    nb_vertices=nv(graph)

    θ_ub = get_expected_upper_bound() #<!> Max number of possible transplants, B ?

    tree = Vector{TreeNode}()

    push!(tree, TreeNode(1, Inf,Vector{Pair{Int,Int}}(),Vector{Pair{Int,Int}}(),Vector{Pair{Int,Int}}(),Vector{Pair{Int,Int}}(), Vector{Int}(), Vector{Int}(), 1, 0))   #the branch and price tree is initialized with the root node

    # initialize the master problem
    initial_time_limit_master_IP = bp_params.time_limit_master_IP
    mastermodel = node_master(instance, column_pool, bp_params, time_limit)
    master_IP = initialize_master_IP(instance, column_pool, bp_params, time_limit)
    bp_params.time_limit_master_IP = initial_time_limit_master_IP

    # initialise branch & price information
    bp_info = BP_info(-Inf,Inf,0)  # LB=-Inf, UB=Inf, nb_col_root=0
    bp_status = BP_status(bp_info, "ON_GOING", -Inf, Inf, Vector{Vector{Int}}(), Vector{Vector{Int}}(), 1, 0.0, 0)

    while length(tree) ≥ 1
        current_node = pop!(tree)
        if current_node.lb ≥ bp_info.UB - ϵ
            continue
        end

        is_fathomed = false

        while(!is_fathomed)

	        # activate the branching constraints of this BP node
	        activate_branching_constraints(mastermodel, current_node, bp_params)

	        master_sol, master_opt = @timeit timer "master" sto_solve_master(mastermodel)

	        deactivate_branching_constraints(mastermodel, current_node, bp_params)

	        if (isnull(master_sol) || master_opt < bp_info.LB)
	        	is_fathomed = true
	        	continue
	        end 

	        fractional_var = check_fractional(master_sol)

	        if !isnull(fractional_var)
	        	branching(fractional_var)
	        	is_fathomed = true
	        	continue
	        end

	        LP_expected_value, dual_var = compute_LP_expected(master_sol)

	        if LP_expected_value ≤ master_sol.θ
	        	add_lshaped_cuts(mastermodel, dual_var)
	        	continue
	        end

	        expected_value = compute_expected(master_sol)

	        obj = master_opt - master_sol.θ + expected_value

	        bp_info.LB = max(bp_info.LB, obj)

	        if expected_value ≥ θ
	        	is_fathomed = true
	        	continue
	        end

	        add_integer_cuts(mastermodel, master_sol, expected_value)
	    end

	    
        
    end
end