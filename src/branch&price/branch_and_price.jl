include("typedef.jl")
include("master.jl")
include("node.jl")
include("subproblem.jl")

function solve(filename::AbstractString, K::Int, L::Int, bp_params::BP_params, timer::TimerOutput = TimerOutput(), time_limit::Float64 = 600.0)
    solve_with_BP(filename, K, L, bp_params, timer, time_limit)
end

"""
    solve_with_BP

This is the main function to call for an execution of the branch-and-price algorithm on input file with given bounds on the length of covering cycles and chains and given options

# Arguments
* `filename::String`: path of the input data files, this should include the name of the files, but not the .dat and .wmd extensions
* `K::Int`: maximum length of cycles
* `L::Int`: maximum length of chains
* `bp_params::BP_params`: solution parameters of the branch-and-price
* `timer::TimerOutput`: a timer that will provide detail of where the computational time was spent during the branch-and-price
* `time_limit::Float64`: time limit of the algorithm, including parsing and prepreprocessing times

#Output parametes
* `instance::Instance`: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains.
* `subgraphs::Graph_copies`: Description of the graph copies of the extended edge formulation
* `bp_status::BP_status`:  Structure containing every relevant information on the execution of the algorithm (including the optimal solution)
"""
function solve_with_BP(filename::String, K::Int, L::Int, bp_params::BP_params = BP_params(), timer::TimerOutput = TimerOutput(), time_limit::Float64 = 600.0)
    start_time = time()
    reset_timer!(timer)
    Random.seed!(10)
    if bp_params.verbose
        printstyled("\n********************************************************************************\n Solve $filename with (K,L) = ($K,$L) using branch-and-price\n - master model uses PIEF = $(bp_params.is_pief) \n - time limit is $time_limit seconds\n********************************************************************************\n\n" ; bold = true, color = :magenta)
    end

    # Parsing
    if bp_params.verbose
    printstyled("\n----------------------------------------------------------\n Parse the input file\n----------------------------------------------------------\n\n" ; color = :yellow) end
    instance = @timeit timer "Parser" Instance(filename, K, L)

    # Preprocessing
    if bp_params.verbose printstyled("\n----------------------------------------------------------\n Preprocessing: compute the graph copies\n----------------------------------------------------------\n\n" ; color = :yellow) end
    subgraphs = @timeit timer "Preprocessing" preprocess_graph_copies(instance, false, bp_params.reduce_vertices, bp_params.fvs)

    # Call the branch-and-price algorithm
    if bp_params.verbose printstyled("\n----------------------------------------------------------\n Solve with branch-and-price\n----------------------------------------------------------\n\n" ; color = :yellow) end
    bp_status = @timeit timer "B&P" branch_and_price(instance, subgraphs, bp_params, timer, time_limit - (time() - start_time))
    bp_status.solve_time = TimerOutputs.time(timer["B&P"])/10^9

    # print the number of cycles and chains for each column's length
    print_and_check_solution(bp_status.best_cycles, bp_status.best_chains, instance, bp_params.verbose)

    # Print cpu profiling
    if bp_params.verbose println(timer) end

    return bp_status, Graph_info(instance), Subgraph_info(subgraphs);
end


"""
    branch_and_price

Core function of the KEP solution with branch-and-price. It requires a parsed instance and the description of the graph copies. The column generation model is that of Riazcos-Alvarez et al (2020), but the many improvements have been added, in particular in the solution of the subproblem.

# Arguments
* `instance::Instance`: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains
* `subgraphs::Graph_copies`: Description of the graph copies of the extended edge formulation
* `bp_params::BP_params`: solution parameters of the branch-and-price
* `timer::TimerOutput`: a timer that will provide detail of where the computational time was spent during the branch-and-price
* `time_limit::Float64`: time limit of the algorithm, including parsing and prepreprocessing times

#Output parametes
* `bp_status::BP_status`:  Structure containing every relevant information on the execution of the algorithm (including the optimal solution)
"""
function branch_and_price(instance::Instance, subgraphs::Graph_copies, bp_params::BP_params, timer::TimerOutput, time_limit::Real)
    # initialization of local variables
    verbose = bp_params.verbose
    graph = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length
    start_time=time()
    nb_vertices=nv(graph)
    nb_subgraph = subgraphs.nb_copies
    column_pool = Vector{Column}()
    if bp_params.is_pief
        initialize_column_pool(instance, column_pool, 2)
        if verbose
            println("Initialize column pool with 2-cycles when using PIEF")
            println("- number of initial columns: $(length(column_pool))\n")
         end
    end
    tree = Vector{TreeNode}()

    push!(tree, TreeNode(1, Inf,Vector{Pair{Int,Int}}(),Vector{Pair{Int,Int}}(),Vector{Pair{Int,Int}}(),Vector{Pair{Int,Int}}(), Vector{Int}(), Vector{Int}(), nv(graph), 0))   #the branch and price tree is initialized with the root node

    # initialize the master problem
    initial_time_limit_master_IP = bp_params.time_limit_master_IP
    mastermodel = node_master(instance, column_pool, bp_params, time_limit)
    master_IP = initialize_master_IP(instance, column_pool, bp_params, time_limit)
    bp_params.time_limit_master_IP = initial_time_limit_master_IP

    # initialise branch & price information
    bp_info = BP_info(-Inf,Inf,0)  # LB=-Inf, UB=Inf, nb_col_root=0
    bp_status = BP_status(bp_info, "ON_GOING", -Inf, Inf, Vector{Vector{Int}}(), Vector{Vector{Int}}(), 1, 0.0, 0, 0, TerminationStatusCode(12))

    while length(tree) >= 1
        current_node = pop!(tree)
        if current_node.ub <= bp_status.bp_info.LB + ϵ
            continue
        end

        # activate the branching constraints of this BP node
        activate_branching_constraints(mastermodel, current_node, bp_params)

         # solve the node relaxation using column generation
        column_flow, pief_flow = @timeit timer "Process_Node" process_node(current_node, instance, mastermodel, subgraphs, bp_status, column_pool, bp_params, master_IP, timer, time_limit - (time() - start_time))

        # update branch & price information
        if bp_status.node_count == 1
            # specific treatment for root node
            bp_info.nb_col_root += length(column_pool)
            if verbose println("After processing root node: LB = $(bp_info.LB), UB = $(current_node.ub)") end

            if bp_info.LB < current_node.ub - ϵ
                if verbose printstyled("- the problem was not solved at root node\n", color=:red) end
                if verbose println("    . deactivate column-disjoint CG and tabu list") end
                # bp_params.is_column_disjoint = false
                # bp_params.is_tabu_list = false
                if bp_params.restart_for_IP
                    if verbose println("    . delete half the columns and restart solving to generate new columns") end
                    deleted_columns = collect(1:length(column_pool))
                    shuffle!(deleted_columns)
                    y = mastermodel[:y]
                    ncols = floor(Int, length(column_pool)/2)
                    for i in deleted_columns[1:ncols]
                        JuMP.set_upper_bound(y[i], 0)
                    end
                    column_flow, pief_flow = @timeit timer "Process_Node" process_node(current_node, instance, mastermodel, subgraphs, bp_status, column_pool, bp_params, master_IP, timer, time_limit - (time() - start_time))
                    for i in deleted_columns
                        JuMP.set_upper_bound(y[i], 1)
                    end
                end
            end
        else
            # deactivate the branching constraints of this BP node for future iterations
            deactivate_branching_constraints(mastermodel, current_node, bp_params)
        end

        # if the relaxation is not infeasible OR not eliminated by bound
        if (!isempty(column_flow) || !isempty(pief_flow)) && current_node.ub > bp_info.LB + ϵ
            branching_done = false
            if (K == 2) && (L == 0)
                total_nb_arcs = 0.0
                for it in column_flow
                    total_nb_arcs += it.second
                end
                for it in pief_flow
                    total_nb_arcs += it.second
                end
                total_nb_arcs = floor(Int, total_nb_arcs+ϵ)
                if total_nb_arcs%2 != 0
                     branch_on_nb_cols(total_nb_arcs, tree, current_node, bp_status.node_count)
                     branching_done = true
                end
            end
            # first, search for a fractional vertex cover to branch on
            is_fractional_vertex = false
            if bp_params.branch_on_vertex && !branching_done
                vertex_to_branch, branching_done = get_branching_vertex(graph, column_flow, pief_flow)

                if  branching_done
                    # branch on the selected fractional vertex
                    branch_on_vertex(vertex_to_branch, mastermodel, tree, current_node, column_pool, bp_status.node_count, bp_params.verbose)
                end
            end

            if !branching_done
                # select a fractional arc to branch on
                arc_to_branch, is_cg_branching = @timeit timer "calc_branch" get_branching_arc(column_flow, pief_flow)

                # branch on the selected fractional arc
                branch_on_arc(arc_to_branch, mastermodel, is_cg_branching, tree, current_node, column_pool, bp_status.node_count, bp_params.verbose)
            end
            bp_status.node_count += 2
        else
            if verbose println("The node is either infeasible or pruned by bound") end
        end

        # Update the global upper bound as the maximum of upper bounds of all open nodes
        bp_info.UB = current_node.ub
        for node in tree
            if node.ub > bp_info.UB
                bp_info.UB = node.ub
            end
        end
        if verbose println("LB = $(bp_info.LB), UB = $(bp_info.UB)") end

        # Stop if the optimality gap is smaller then the tolerance
        if bp_info.UB < bp_info.LB + ϵ
            bp_status.status="OPTIMAL"
            break
        end

        # Stop if the time limit is exceeded
        if time() - start_time >= time_limit
            if verbose println("\e[35m The time limit is exceeded \e[00m") end
            bp_status.status="TIME_LIMIT"
            break
        end
    end

    # Update solution status and gap
    if length(tree) == 0
        bp_status.status = "OPTIMAL"
    end
    if verbose
        printstyled("\n----------------------------------------------------------\n The execution of the branch-and-price is complete\n" ; color = :yellow)
        if bp_status.status == "OPTIMAL"
            printstyled("- the solution is optimal\n" ; color = :yellow)
        elseif bp_status.status == "TIME_LIMIT"
            printstyled("- the time limit is exceeded\n" ; color = :yellow)
        end
    end
    bp_status.objective_value = bp_info.LB
    bp_status.relative_gap = abs(bp_info.UB - bp_info.LB) / abs(bp_info.LB + 10^-10)
    if verbose
        printstyled("- best solution found: value $(bp_status.objective_value) with gap $(100 * bp_status.relative_gap) %\n" ; color = :yellow)
        printstyled("----------------------------------------------------------\n\n" ; color = :yellow)
    end

    return bp_status
end

"""
    get_branching_arc
Find a fractional arc to branch. The fractional arc closest to 0.5 will be
selected to branch, if there is no fractional arc in the solution

# Input parameters
* `column_flow::Dict{Pair{Int,Int}, Float64}` : The dictionary containing the nonzero flow on each arc due to the selection of columns
* `pief_flow::Dict{Pair{Int,Int}, Float64}` : The dictionary containing the flow on each arc from the pief model for chain search (if applicable)

# Output parameters
* `arc_to_branch::Pair{Int,Int}` : The arc to be branched
* `is_cg_branching::Bool`: True if the branching impacts the subproblem of the coluln generation, false if it impacts only the master problem
"""
function get_branching_arc(column_flow::Dict{Pair{Int,Int}, Float64}, pief_flow::Dict{Pair{Int,Int}, Float64})
    arc_to_branch = (0=>0)
    is_cg_branching = false  # true if the branching impacts subproblems of column generation
    val = 0.5  # measure of fractionality
    # first search for a branching on the arcs included in the columns
    for it in column_flow
        if abs(it.second-0.5) < val - ϵ
            val = abs(it.second-0.5)
            arc_to_branch = it.first
            is_cg_branching = true
            if val < ϵ break end
        end
    end

    # if no arc with fractional column flow was found, search for an arc with a fractional chain flow in the solution of the pief model
    if !is_cg_branching
        for it in pief_flow
            if abs(it.second-0.5) < val - ϵ
                val = abs(it.second-0.5)
                arc_to_branch = it.first
                if val < ϵ break end
            end
        end
    end

    if arc_to_branch == (0=>0)
        error("There should always be a fractional arc when entering this function")
    end

    return arc_to_branch, is_cg_branching
end

"""
    get_branching_vertex
Find a fractional vertex cover to branch. The fractional vertex closest to 0.5 will be selected to branch

# Input parameters
* `column_flow::Dict{Pair{Int,Int}, Float64}` : The dictionary containing the nonzero flow on each arc due to the selection of columns
* `pief_flow::Dict{Pair{Int,Int}, Float64}` : The dictionary containing the flow on each arc from the pief model for chain search (if applicable)

# Output parameters
* `vertex_to_branch::Int` : The vertex to be branched on, nothing if none was found
"""
function get_branching_vertex(graph::SimpleDiGraph, column_flow::Dict{Pair{Int,Int}, Float64}, pief_flow::Dict{Pair{Int,Int}, Float64})
    vertex_to_branch = 0
    is_fractional_vertex = false
    val = 0.5
    vertex_cover = zeros(nv(graph))
    # compute the total cover of each vertex
    for it in column_flow
        arc = it.first
        vertex_cover[arc[2]] += it.second
    end
    for it in pief_flow
        arc = it.first
        vertex_cover[arc[2]] += it.second
    end

    # get the most fractional cover
    for v in vertices(graph)
        if abs(vertex_cover[v]-0.5) < val - ϵ
            val = abs(vertex_cover[v]-0.5)
            vertex_to_branch = v
            is_fractional_vertex = true
            if val < ϵ break end
        end
    end

    return vertex_to_branch, is_fractional_vertex
end


"""
    branch_on_arc
Update the branch-and-bound tree with two new nodes by branching on the given arc with the specified branching (only on master problem or both in master and in subproblem)

# Input parameters
* `arc_to_branch::Pair{Int,Int}` : The arc to be branched
* `master::Model`: Master model where the branching constraints are added
* `is_cg_branching::Bool`: True if the branching impacts the subproblem of the coluln generation, false if it impacts only the master problem
* `tree::Vector{TreeNode}`: Branch-and-bound tree
* `current_node::TreeNode`: Branch-and-bound node currently treated
* `column_pool::Vector{Column}`: Pool of all columns in current master problem
* `node_count::Int`: Number of BP nodes enumerated until now
"""
function branch_on_arc(arc_to_branch::Pair{Int,Int}, master::Model,  is_cg_branching::Bool, tree::Vector{TreeNode}, current_node::TreeNode, column_pool::Vector{Column}, node_count::Int, verbose::Bool = true)
    y = master[:y]
    slack = maxter[:slack]
    if is_cg_branching
        if verbose println("Two new nodes are created by branching on variable column_flow[$(arc_to_branch.first), $(arc_to_branch.second)]") end

        # add each branching node in the tree
        node_zero = TreeNode(current_node)
        node_zero.index = node_count + 1
        push!(node_zero.setzero, arc_to_branch)
        push!(tree, node_zero)
        node_one = TreeNode(current_node)
        node_one.index = node_count + 2
        push!(node_one.setone, arc_to_branch)
        push!(tree, node_one)

        # add the corresponding constraints in the master model, but deactivate them by default: they will activated only in relevant BP nodes
        master[:branch_one][arc_to_branch] = @constraint(master, sum(y[c] for c in 1:length(column_pool) if arc_to_branch in (column_pool[c]).arcs) + slack >= 0, base_name = "branch_one[$arc_to_branch]")
        master[:branch_zero][arc_to_branch] = @constraint(master, sum(y[c] for c in 1:length(column_pool) if arc_to_branch in (column_pool[c]).arcs) <= 1, base_name = "branch_zero[$arc_to_branch]")
    else
        if verbose println("Branching on an arc of the master problem: ($(arc_to_branch.first), $(arc_to_branch.second))") end

        # add each branching node in the tree
        node_zero = TreeNode(current_node)
        node_zero.index = node_count + 1
        push!(node_zero.setzero_pief, arc_to_branch)
        push!(tree, node_zero)
        node_one = TreeNode(current_node)
        node_one.index = node_count + 2
        push!(node_one.setone_pief, arc_to_branch)
        push!(tree, node_one)

        # add the corresponding constraints in the master model, but deactivate them by default: they will be activated only in relevant BP nodes
        chain_flow = model[:chain_flow]
        master[:branch_one_pief][arc_to_branch] = @constraint(master, sum(chain_flow[e[1],e[2],k] for k in 1:L) + slack >= 0, base_name = "branch_one_pief[$arc_to_branch]")
        master[:branch_zero_pief][arc_to_branch] = @constraint(master, sum(chain_flow[e[1],e[2],k] for k in 1:L) <= 1, base_name = "branch_zero_pief[$arc_to_branch]")
    end
end

function branch_on_nb_cols(total_nb_arcs::Int, tree::Vector{TreeNode}, current_node::TreeNode, node_count::Int, verbose::Bool = true)
    if verbose println("Two new nodes are created by branching on the total number of arcs to get an even number") end

    # add each branching node in the tree
    node_max = TreeNode(current_node)
    node_max.index = node_count + 1
    node_max.nb_cols_max = (total_nb_arcs - 1)/2
    push!(tree, node_max)
    node_min = TreeNode(current_node)
    node_min.index = node_count + 2
    node_min.nb_cols_min = (total_nb_arcs + 1)/2
    push!(tree, node_min)
end

"""
    branch_on_vertex
Update the branch-and-bound tree with two new nodes by branching on the given vertex.

# Input parameters
* `vertex_to_branch::Int` : The arc to be branched
* `master::Model`: Master model where the branching constraints are added
* `tree::Vector{TreeNode}`: Branch-and-bound tree
* `current_node::TreeNode`: Branch-and-bound node currently treated
* `column_pool::Vector{Column}`: Pool of all columns in current master problem
** `node_count::Int`: Number of BP nodes enumerated until now
"""
function branch_on_vertex(vertex_to_branch::Int, master::Model,  tree::Vector{TreeNode}, current_node::TreeNode, column_pool::Vector{Column}, node_count::Int, verbose::Bool = true)
    y = master[:y]
    if verbose println("Two new nodes are created by branching on vertex cover $vertex_to_branch") end

    # add each branching node in the tree
    node_zero = TreeNode(current_node)
    node_zero.index = node_count + 1
    push!(node_zero.setzero_vertex, vertex_to_branch)
    push!(tree, node_zero)
    node_one = TreeNode(current_node)
    node_one.index = node_count + 2
    push!(node_one.setone_vertex, vertex_to_branch)
    push!(tree, node_one)

    # add the corresponding constraints in the master model, but deactivate them by default: they will activated only in relevant BP nodes
    master[:branch_one_vertex][vertex_to_branch] = @constraint(master, sum(y[c] for c in 1:length(column_pool) if vertex_to_branch in (column_pool[c]).vertices) >= 0, base_name = "branch_one_vertex[$vertex_to_branch]")
    master[:branch_zero_vertex][vertex_to_branch] = @constraint(master, sum(y[c] for c in 1:length(column_pool) if vertex_to_branch in (column_pool[c]).vertices) <= 1, base_name = "branch_zero_vertex[$vertex_to_branch]")
end

"""
    initialize_column_pool
Initialize the pool of columns for the branch-and-price by enumerating all k-cycles up to an input maximum value of k.

# Input parameters
*`graph::SimpleDiGraph`: The KEP graph
*`column_pool::Vector{Column}`: Pool of columns where initial columns will be pushed
*`max_cycle_length::Int`: Maximum length of the cycles to be enumerated at initialization (0 if initialization is deactivated)

# Output parameters
* `column_pool::Vector{Column}` : set of columns initially added to the master problem
"""
function initialize_column_pool(instance::Instance, column_pool::Vector{Column}, max_cycle_length::Int = 0)
    graph = instance.graph
    if max_cycle_length >= 2
        max_nb_cols = 20
        nb_cols = zeros(nv(graph))
        for v in instance.pairs
            col_count = 0
            vtx_list = findall(instance.edge_weight[:,v] .!= 0.0)
            shuffle!(vtx_list)
            for u in vtx_list
                if nb_cols[u] >= max_nb_cols continue end
                if u < v && has_edge(graph, v, u)
                    if nb_cols[v] >= max_nb_cols continue end
                    nb_cols[u] += 1
                    nb_cols[v] += 1
                    path = Vector{Int}(undef, 0)
                    push!(path, u)
                    push!(path, v)
                    if instance.is_vertex_weighted
                        push!(column_pool, Column(path, instance.vertex_weight, true))
                    else
                        push!(column_pool, Column(path, instance.edge_weight, true))
                    end
                    if nb_cols[u] >= max_nb_cols continue end
                end
            end
        end
    end

    return column_pool
end
