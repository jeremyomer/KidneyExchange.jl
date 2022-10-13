"""
    process_node

Columns generation of the node

# Input parameter
- `tree_node::TreeNode`: Branch-and-bound node to solve with column generation

# Output parameter
- `column_flow::Dict{Pair{Int,Int}, Float64}`: value of x_(i,j) of arc (i->j)`
"""
function process_node(tree_node::TreeNode, instance::Instance, mastermodel::Model, subgraphs::Graph_copies, bp_status::BP_status, column_pool::Vector{Column},  bp_params::BP_params, master_IP::Model, timer::TimerOutput, time_limit::Float64)
    # local variables
    graph = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length
    nb_vertices = nv(graph)
    verbose = bp_params.verbose

    # verbose outputs
    if verbose
        printstyled("Processing node $(tree_node.index)\n\n" ; color = :yellow)
        if !isempty(tree_node.setzero_vertex)
            println("- the vertices set to zero are:", tree_node.setzero_vertex)
        end
        if !isempty(tree_node.setone_vertex)
            println("- the vertices set to one are:", tree_node.setone_vertex)
        end
        if !isempty(tree_node.setzero)
            println("- the arcs set to zero are:", tree_node.setzero)
        end
        if !isempty(tree_node.setone)
            println("- the arcs set to one are:", tree_node.setone)
        end
        if (bp_params.is_pief)
            if !isempty(tree_node.setzero_pief)
                println("- the arcs set to zero in pief are:", tree_node.setzero_pief)
            end
            if !isempty(tree_node.setone_pief)
                println("- the arcs set to one in pief are:", tree_node.setone_pief)
            end
        end
        if tree_node.nb_cols_max < nv(graph)
            println("- maximum number of arcs in solution $(2*tree_node.nb_cols_max)")
        end
        if tree_node.nb_cols_min > 0
            println("- minimum number of arcs in solution $(2*tree_node.nb_cols_min)")
        end
    end

    # start by solving the relaxed master problem
    @timeit timer "opt_master" optimize!(mastermodel)

    # initialize local variables specific to the algorithm
    arc_cost = zeros(nb_vertices, nb_vertices)

    to_optimality = true
    tabu_list = Vector{Int}()
    if bp_params.is_tabu_list
        to_optimality = false
    end
    subproblems_for_mip = Vector{Int}()  # list of subproblems that need to be solved optimally if no chain nor cycle was found
    column_flow = Dict{Pair{Int,Int}, Float64}()  # flow on arcs due to the selection of columns (cycles or chains)
    pief_flow = Dict{Pair{Int,Int}, Float64}()  # flow on arcs from the chain flow variables of the pief master model

    # preallocated vectors for Bellman speedup
    pred_for_bellman = Vector{Vector{Int}}(undef, max(K,L))
    for k = 1:(max(K,L))
        pred_for_bellman[k] = zeros(Int, nv(graph))
    end
    d_for_bellman = Vector{Float64}(undef, nv(graph))
    vertex_cost = zeros(nb_vertices)

    # execute the column generation
    nb_iter = 0
    while true
        node_infeasible=false
        nb_iter += 1
        if verbose printstyled("\nIteration $nb_iter:\n") end
        # ==============================================================
        #  Collect values from the solution of the restricted master problem
        # ==============================================================
        if JuMP.termination_status(mastermodel) != MOI.OPTIMAL
            # println(mastermodel)
            return Dict{Pair{Int,Int}, Float64}(), Dict{Pair{Int,Int}, Float64}()
            error("The master problem relaxation has not been solved to optimality")
        end
        master_value = JuMP.objective_value(mastermodel)
        if verbose println("- current master value: ", master_value) end

        # dual values of the vertex-disjoint constraints
        if !JuMP.has_duals(mastermodel)
            # println(mastermodel)
            error("The master problem has no dual solution")
        end
        λ = JuMP.shadow_price.(mastermodel[:capacity])

        # get the dual values of the arc branching constraints
        δ_one = Dict{Pair{Int64, Int64}, Float64}()
        δ_zero = Dict{Pair{Int64, Int64}, Float64}()
        for arc in keys(mastermodel[:branch_one])
            δ_one[arc] = -shadow_price(mastermodel[:branch_one][arc])
            # JuMP has weird definitions of its shadow_price and dual function, pay attention to potential errors
            if δ_one[arc] > 0
                δ_one[arc] = -δ_one[arc]
            end
        end
        for arc in keys(mastermodel[:branch_zero])
            δ_zero[arc] = shadow_price(mastermodel[:branch_zero][arc])
            # JuMP has weird definitions of its shadow_price and dual function, pay attention to potential errors
            if δ_zero[arc] < 0
                δ_zero[arc] = -δ_zero[arc]
            end
        end

        # get the dual values of the vertex branching constraints
        δ_one_vertex = Dict{Int, Float64}()
        δ_zero_vertex = Dict{Int, Float64}()
        for v in keys(mastermodel[:branch_one_vertex])
            δ_one_vertex[v] = -shadow_price(mastermodel[:branch_one_vertex][v])
            # JuMP has weird definitions of its shadow_price and dual function, pay attention to potential errors
            if δ_one_vertex[v] > 0
                δ_one_vertex[v] = -δ_one_vertex[v]
            end
        end
        for v in keys(mastermodel[:branch_zero_vertex])
            δ_zero_vertex[v] = shadow_price(mastermodel[:branch_zero_vertex][v])
            # JuMP has weird definitions of its shadow_price and dual function, pay attention to potential errors
            if δ_zero_vertex[v] < 0
                δ_zero_vertex[v] = -δ_zero_vertex[v]
            end
        end

        # get the duals of the branching constraint on the number of arcs (useful only for K=2, L=0)
        δ_nb_arcs_max = shadow_price(mastermodel[:branch_nb_arcs_max])
        δ_nb_arcs_min = shadow_price(mastermodel[:branch_nb_arcs_min])
        # JuMP has weird definitions of its shadow_price and dual function, pay attention to potential errors
        if δ_nb_arcs_max < 0
            δ_nb_arcs_max = -δ_nb_arcs_max
        end
        if δ_nb_arcs_min > 0
            δ_nb_arcs_min = -δ_nb_arcs_min
        end
        δ_nb_arcs = δ_nb_arcs_max + δ_nb_arcs_min

        # ===============================================================
        #  Check whether the solution is integer or not and recover the corresponding solution if it is integer
        # ===============================================================
        is_integer = true
        column_val = JuMP.value.(mastermodel[:y])
        selected_cycles = Vector{Vector{Int}}()
        selected_chains = Vector{Vector{Int}}()
        for c in 1:length(column_pool)
            if column_val[c] > 1 - ϵ
                if column_pool[c].is_cycle
                    push!(selected_cycles, (column_pool[c]).vertices)
                else
                    push!(selected_chains, (column_pool[c]).vertices)
                end
            elseif column_val[c] > ϵ
                is_integer = false
                break
            end
        end

        selected_arcs_pief = Vector{Pair{Int,Int}}()
        if bp_params.is_pief && is_integer && L >= 1
            # Very odd behaviour of the JuMP.value function here; it takes a VERY long time to get the solution although it is done instantaneously with the corresponding IP, I could not find why, so just deactivated this opportunistic search for interger solution, which is mostly useless anyways
            is_integer = false
            if false
                for u in vertices(graph)
                    for v in outneighbors(graph, u)
                        pief_flow[u=>v] = sum(value(mastermodel[:chain_flow][u,v,k]) for k in 1:L)
                        if pief_flow[u=>v] > 1-ϵ
                            push!(selected_arcs_pief, u=>v)
                        elseif pief_flow[u=>v] > ϵ
                            is_integer = false
                            break
                        end
                    end
                end
                println("done checking integrality")
                if is_integer
                    println("check integer solution")
                    selected_chains = compute_chains_from_pief(mastermodel, instance)
                end
            end
        end

        # ===============================================================
        #  Update BP information if a new incumbent has been found
        # ===============================================================
        if is_integer &&  (master_value > bp_status.bp_info.LB + ϵ)
            bp_status.bp_info.LB = master_value
            bp_status.best_cycles = selected_cycles
            bp_status.best_chains = selected_chains
            if verbose println("\e[32m New incumbent with value $master_value found during the solution of the restricted master \e[00m") end
        end


        # ===============================================================
        #  Solve the subproblem
        # ===============================================================
        # arc costs need to be updated only if not at root node
        max_cost = 0.0
        if !instance.is_vertex_weighted || !isempty(tree_node.setone) || !isempty(tree_node.setzero)
            max_cost = calculate_arc_cost(instance, arc_cost, λ, δ_one, δ_zero, δ_one_vertex, δ_zero_vertex)
            arc_cost_trans = permutedims(arc_cost)
        else
            vertex_cost .= instance.vertex_weight
            vertex_cost .-= λ
            if !isempty(tree_node.setone_vertex) || !isempty(tree_node.setzero_vertex)
                for v in keys(δ_one_vertex)
                    vertex_cost[v] -= δ_one_vertex[v]
                end
                for v in keys(δ_zero_vertex)
                    vertex_cost[v] -= δ_zero_vertex[v]
                end
            end
        end

        # 1. Search only for cycles at first
        cycle_added = false
        nb_cycle_added = 0

        # initialize disjoint column generation
        is_covered = falses(nb_vertices)
        nb_covered = zeros(nb_vertices)

        # shuffle the indices of the subproblems to treat if using column disjoint CG, otherwise, we will always focus on the same subproblems
        inds = collect((instance.nb_altruists+1):length(subgraphs.sources))
        if bp_params.is_column_disjoint
            shuffle!(inds)
        end
        # solve the subproblems
        for l in inds
            vstar = subgraphs.sources[l]  # source of current copy

            # skip this subproblem if the source is already covered or in the tabu list
            if vstar in tabu_list || is_covered[vstar]
                continue
            end

            # initialize the list of vertices of the copy
            is_not_vertex = is_covered  .| .!subgraphs.is_vertex_list[l]

            # in the presence of branching constraints on the arcs covered by columns, there are dual variables on arcs, so a different cycle search must be called
            if !instance.is_vertex_weighted || !isempty(tree_node.setzero) || !isempty(tree_node.setone)
                path = @timeit timer "Bellman-Ford" Bellman_Ford_cycle_search(graph, arc_cost, arc_cost_trans, max_cost, vstar, K, is_not_vertex, pred_for_bellman, d_for_bellman, -δ_nb_arcs)
            else
                path = @timeit timer "Bellman-Ford" Bellman_Ford_cycle_search(graph, vertex_cost, vstar, K, is_not_vertex, pred_for_bellman, d_for_bellman, -δ_nb_arcs)
            end

            # add one column for each positive cycle we just found
            if !isempty(path)
                cycle_added = true
                nb_cycle_added += 1
                if instance.is_vertex_weighted
                    column = Column(path, instance.vertex_weight, true)
                else
                    column = Column(path, instance.edge_weight, true)
                end
                push!(column_pool, column)
                # push!(node_columns, column)
                if bp_params.is_column_disjoint
                    for v in column.vertices
                        nb_covered[v] += 1
                        if nb_covered[v] == bp_params.max_intersecting_columns
                            is_covered[v] = true
                        end
                    end
                end
                # add column to master problem model
                @timeit timer "opt_master" add_column_to_master(column, mastermodel, tree_node)
                @timeit timer "IP_master" add_column_to_master_IP(column, master_IP)
            else
                # make sure that the subproblem will not be solved until we need to prove optimality if it did not produce any new column
                if bp_params.is_tabu_list
                    push!(tabu_list, vstar)
                end
            end
        end
        if verbose println("- nb of cycles added = ", nb_cycle_added) end

        # 2. Search for positive chain if no cycle was added unless we are solving the pief model
        chain_added = false
        if !bp_params.is_pief
            inds = collect(1:instance.nb_altruists)
            if bp_params.is_column_disjoint
                shuffle!(inds)
            end

            # 1. first try and find positive chains with an adaptation of the Bellman algorithm; if L >= K+2, this algorithm cannot guarantee that the optimality of column generation is reached, so we may need to solve some subproblems to optimality using a mip formulation
            is_positive_cycle = false
            empty!(subproblems_for_mip)
            nb_chains_added = 0
            for l in inds
                vstar = subgraphs.sources[l]
                if vstar in tabu_list
                    continue
                end
                is_not_vertex = is_covered  .| .!subgraphs.is_vertex_list[l]

                # in the presence of branching constraints on the arcs covered by columns, there are dual variables on arcs, so a different chain search must be called
                if !instance.is_vertex_weighted || !isempty(tree_node.setone) || !isempty(tree_node.setzero)
                    if to_optimality && !cycle_added
                        positive_path, is_positive_cycle = @timeit timer "Bellman-Ford-chain" Bellman_Ford_chain_search_optimality(graph, arc_cost_trans, max_cost, vstar, L, K, is_not_vertex, pred_for_bellman, d_for_bellman, -δ_nb_arcs)
                    else
                        positive_path = @timeit timer "Bellman-Ford-chain" Bellman_Ford_chain_search(graph, arc_cost_trans, max_cost, vstar, L, K, is_not_vertex, pred_for_bellman, d_for_bellman, -δ_nb_arcs)
                    end
                else
                    if to_optimality && !cycle_added
                        positive_path, is_positive_cycle = @timeit timer "Bellman-Ford-chain" Bellman_Ford_chain_search_optimality(graph, vertex_cost, vstar, L, K, is_not_vertex, pred_for_bellman, d_for_bellman, -δ_nb_arcs)
                    else
                        positive_path = @timeit timer "Bellman-Ford-chain" Bellman_Ford_chain_search(graph, vertex_cost, vstar, L, K, is_not_vertex, pred_for_bellman, d_for_bellman, -δ_nb_arcs)
                    end
                end

                if !isempty(positive_path)
                    # add one column for the positive chain that was just found
                    chain_added = true
                    nb_chains_added += 1
                    if instance.is_vertex_weighted
                        column = Column(positive_path, instance.vertex_weight, false)
                    else
                        column = Column(positive_path, instance.edge_weight, false)
                    end
                    push!(column_pool, column)
                    # push!(node_columns, column)
                    if bp_params.is_column_disjoint
                        for v in column.vertices
                            nb_covered[v] += 1
                            if nb_covered[v] == bp_params.max_intersecting_columns
                                is_covered[v] = true
                            end
                        end
                    end
                    # add column to master problem model
                    @timeit timer "opt_master" add_column_to_master(column, mastermodel, tree_node)
                    @timeit timer "IP_master" add_column_to_master_IP(column, master_IP)
                else
                    if bp_params.is_tabu_list
                        # make sure that the subproblem will not be solved until we need to prove optimality if it did not produce any new column
                        push!(tabu_list, vstar)
                    end
                    if is_positive_cycle
                        # record the subproblem as not solved to optimality if a positive cycle was found and no positive path
                        push!(subproblems_for_mip, l)
                    end
                end
            end
            if verbose println("nb of chains added = $(nb_chains_added)") end
        end

        # 3. if chains or cycles were added, solve the master and step to next iteration
        if chain_added || cycle_added
            @timeit timer "opt_master" optimize!(mastermodel)
            if to_optimality && bp_params.is_tabu_list
                to_optimality = false
            end
            continue
        end

        # 4. if no column was added, but some subproblems were ignored, solve every subproblem from now on (every subproblem is solved at first iteration if no column is generated)
        if !to_optimality && nb_iter >= 2
            if verbose println("- no positive column was added, switch to optimality search") end
            to_optimality = true
            if !isempty(tabu_list)
                empty!(tabu_list)
                continue
            end
        end

        # 5. if every subproblem was solved, but no column added, we may not have reached optimality if L >= K+2 and some positive cycles were found while looking for chains
        if !bp_params.is_pief && (L >= K + 2) && !isempty(subproblems_for_mip)
            if verbose printstyled("- solve $(length(subproblems_for_mip)) subproblems to optimality as MIPs\n" ; color = :red) end
            positive_paths = Vector{Vector{Int}}()
            # initialize the mip for chain search if not already done
            if num_variables(subgraphs.chain_mip) == 0
                # initialize the MIP for chain search if L >= K+2
                @timeit timer "create_chain_mip" subgraphs.chain_mip = create_chain_mip(graph,  L, bp_params.optimizer, time_limit)
            end

            # compute are reduced costs
            max_cost = (instance, arc_cost, λ, δ_one, δ_zero, δ_one_vertex, δ_zero_vertex)

            # solve each subproblem independently
            for l in subproblems_for_mip
                # need to initialize mips first
                is_positive_chain, chain = MIP_chain_search(subgraphs.chain_mip, graph, subgraphs.sources[l], subgraphs.is_vertex_list[l], arc_cost, verbose)
                if (is_positive_chain)
                    push!(positive_paths, chain)
                end
            end
            # add one column for each positive chain we just found
            for path in positive_paths
                chain_added = true
                nb_chains_added += 1
                if instance.is_vertex_weighted
                    column = Column(path, instance.vertex_weight, false)
                else
                    column = Column(path, instance.edge_weight, false)
                end
                push!(column_pool, column)

                # add column to master problem model
                @timeit timer  "opt_master" add_column_to_master(column, mastermodel, tree_node)
                @timeit timer "IP_master" add_column_to_master_IP(column, master_IP)
            end
            if chain_added
                if verbose println("- found at least one positive column") end
                @timeit timer "opt_master" optimize!(mastermodel)
                continue
            end
        end

        if !cycle_added && !chain_added
            if verbose printstyled("\nNode relaxation is solved to optimality\n" ; color = :green) end
            nodeub = master_value
            if nodeub < tree_node.ub
                if abs(nodeub - round(nodeub)) < ϵ
                    nodeub = round(nodeub)
                else
                    nodeub = floor(nodeub)
                end
                tree_node.ub = nodeub
            end
            if verbose print("- node upper bound is $(tree_node.ub), ") end
            if verbose println("tree lower bound is $(bp_status.bp_info.LB)") end

            #check whether artificial column is used
            if tree_node.index >= 2
                if JuMP.value(mastermodel[:slack]) > ϵ
                    node_infeasible = true
                    if verbose printstyled("\n Node is infeasible\n" ; color = :green) end
                    # if the node is infeasible, return an empty vector
                    return Dict{Pair{Int,Int}, Float64}(), Dict{Pair{Int,Int}, Float64}()
                end
            end
            # if not compute the arc flows resulting from the selction of columns
            column_flow = compute_arc_flow(column_val, column_pool)
            break
        end
    end

    # find a feasible solution at current node by solving the master IP
    if bp_params.solve_master_IP &&
        (bp_status.bp_info.LB < tree_node.ub - ϵ) &&
        ( (bp_status.termination_status_last_ip != OPTIMAL) ||
            (length(column_pool) > bp_status.nb_cols_last_ip) )
        if verbose 
            if length(column_pool) > bp_status.nb_cols_last_ip
                println("\nThe number of columns increased: ")
            elseif bp_status.termination_status_last_ip != OPTIMAL
                println("\nLast solution of master IP did not reach optimality: ")
            end
            println(" search for a feasible solution at node $(tree_node.index)")
        end

        @timeit timer "IP_master" solve_master_IP(master_IP, column_pool, instance, bp_status, bp_params)
        bp_status.nb_cols_last_ip = length(column_pool)
        bp_status.node_count_last_ip = bp_status.node_count
        bp_status.termination_status_last_ip = termination_status(master_IP)
        if verbose
            println("Termination status : ", bp_status.termination_status_last_ip)
        end

        if (termination_status != OPTIMAL)
            bp_params.time_limit_master_IP = min(bp_params.time_limit_master_IP + 10.0, 300.0)
            set_time_limit(master_IP, bp_params.time_limit_master_IP, bp_params.optimizer)
        end
    end

    return column_flow, pief_flow
end


"""
    get_feasible_solution

Extract a integer feasible solution from the fractional solution by conserving a set of vertex-disjoint cycles or chains

# Arguments
* `fractional_solution::Array{Array{Float64,1},1}`: The set of value of variables indicating whether or not cycles and chains are selected
* `node_columns::Array{Array{Column,1},1}`: The set of cycles and chains associated with variables

# Output parameters
* `feas_val::Real`:  The objective value of the feasible solution
* `columns::Array{Array{Int,1},1}`: The set of selected cycles and chains of the feasible solution
"""
function get_feasible_solution(fractional_solution::Vector{Float64}, node_columns::Vector{Column})
    columns = Vector{Vector{Int}}()
    column_weights = Vector{Float64}() # weight of the corresponding cycle
    sol_vals = Vector{Float64}() # the fractional value of selected columns

    # Add progressively vertex-disjoint cycles and chains into the set of selected columns
    ind_sorted = sortperm(fractional_solution, rev = true)
    for c in ind_sorted
        if fractional_solution[c] > ϵ
            column = node_columns[c]
            column_vertices = column.vertices
            column_weight = column.weight
            add = true
            # check if there exits columns in the set of selected columns that have common vertices with the cycle to be added
            redundant_cycle_indices = check_used_vertices(column_vertices, columns)
            # if not, add the cycle, if yes, check if it is worthy to add the cycle and delete the columns that have common vertices
            if !isempty(redundant_cycle_indices)
                # if the cycle contributed more in terms of the objective value then remove the columns that have common vertices and add the cycle
                if fractional_solution[c] * column_weight > sum(sol_vals[i] * column_weights[i] for i in redundant_cycle_indices)
                    deleteat!(columns, redundant_cycle_indices)
                    deleteat!(column_weights, redundant_cycle_indices)
                    deleteat!(sol_vals, redundant_cycle_indices)
                else
                    add = false
                end
            end
            if add
                push!(columns, column_vertices)
                push!(column_weights, column_weight)
                push!(sol_vals, fractional_solution[c])
            end
        else
            break
        end
    end
    return sum(column_weights), columns
end

"""
    check_used_vertices

Gives a set of indices of cycle in cycles who have common vertices with cycle

# Arguments
* `cycle::Array{Int,1}`: The vertices of cycle
* `cycles::Array{Array{Int,1},1}`: The set of cycles

#Output parametes
* `redundant_cycle_indices::Array{Int,1}`: the set of indices of cycle in cycles who have common vertices with cycle

"""
function check_used_vertices(cycle::Vector{Int}, cycles::Vector{Vector{Int}})
    redundant_cycle_indices = Vector{Int}()
    vertices_to_be_checked = copy(cycle)
    for c in 1:length(cycles)
        checked_vertices_id = Vector{Int}()
        for i in 1:length(vertices_to_be_checked)
            if vertices_to_be_checked[i] in cycles[c]
                push!(redundant_cycle_indices, c)
                push!(checked_vertices_id, i)
            end
        end
        deleteat!(vertices_to_be_checked, checked_vertices_id)
    end
    return sort(unique(redundant_cycle_indices))
end

"""
add_column_to_master

# Arguments
* `column::Column`: Column to add to the master model
* `mastermodel::Model`: JuMP model for current master problem
* `treenode::TreeNode`: Information on current tree node
# Return values: None
"""

function add_column_to_master(column::Column, mastermodel::Model, tree_node::TreeNode)
    # add the new variable
    y = mastermodel[:y]
    push!(y, @variable(mastermodel, lower_bound=0))
    set_name(y[end], "y_$(length(y))")

    # set the coefficient of the variable in the capacity constraint and objective function
    capacity = mastermodel[:capacity]
    for v in column.vertices
        set_normalized_coefficient(capacity[v], y[end], 1)
    end
    set_objective_coefficient(mastermodel, y[end], column.weight)

    # set the coefficient of the new column in every branching constraint (including those that are not active at this node)
    set_normalized_coefficient(mastermodel[:branch_nb_arcs_max], y[end], 1)
    set_normalized_coefficient(mastermodel[:branch_nb_arcs_min], y[end], 1)
    branch_one = mastermodel[:branch_one]
    for arc in keys(branch_one)
        if arc in column.arcs
            set_normalized_coefficient(branch_one[arc], y[end], 1)
        end
    end
    branch_zero = mastermodel[:branch_zero]
    for arc in keys(branch_zero)
        if arc in column.arcs
            set_normalized_coefficient(branch_zero[arc], y[end], 1)
        end
    end
    branch_one_vertex = mastermodel[:branch_one_vertex]
    for v in keys(branch_one_vertex)
        if v in column.vertices
            set_normalized_coefficient(branch_one_vertex[v], y[end], 1)
        end
    end
    branch_zero_vertex = mastermodel[:branch_zero_vertex]
    for v in keys(branch_zero_vertex)
        if v in column.vertices
            set_normalized_coefficient(branch_zero_vertex[v], y[end], 1)
        end
    end
end

function add_column_to_master_IP(column::Column, master_IP::Model)
    # add the new variable
    y = master_IP[:y]
    push!(y, @variable(master_IP, binary = true))
    set_name(y[end], "y_$(length(y))")

    # set the coefficient of the variable in the capacity constraint and objective function
    capacity = master_IP[:capacity]
    for v in column.vertices
        set_normalized_coefficient(capacity[v], y[end], 1)
    end
    set_objective_coefficient(master_IP, y[end], column.weight)
end

"""
    compute_arc_flow

Calculates the relaxed value of the decision variable x_(i,j) of arc (i->j) in A.
x_(i,j) whose value >0 is stored in x::Dict{Pair{Int,Int}, Float64}

# Arguments
* `mastersol::Vector{Float64}`: solution value of the master problem, ie: y[c] for c in cycles
* `node_columns::Vector{Column}`: The corresponding cycles of current node

# Return values
* `x:: Dict{Pair{Int,Int}, Float64}`: value of x_(i,j) of arc (i->j)
"""
function compute_arc_flow(mastersol::Vector{Float64}, node_columns::Vector{Column})
    x = Dict{Pair{Int,Int}, Float64}()
    for c in 1:length(mastersol)
        if mastersol[c] > ϵ
            column = node_columns[c]
            for arc in column.arcs
                if haskey(x, arc)
                    x[arc] = x[arc] + mastersol[c]
                else
                    push!(x, (arc => mastersol[c]))
                end
            end
        end
    end
    return x
end

function compute_chains_from_pief(master::Model, instance::Instance)
    graph = instance.graph
    L = instance.max_chain_length
    selected_chains = Vector{Vector{Int}}()
    for u in instance.altruists
        for v in outneighbors(graph, u)
            chain = Vector{Int}()
            if JuMP.value(master[:chain_flow][u,v,1]) > 1 - ϵ
                push!(chain, u)
                push!(chain, v)
                k = 2
                while k <= L
                    is_end_of_chain = true
                    for w in outneighbors(graph, v)
                        if JuMP.value(master[:chain_flow][v,w,k]) > 1 - ϵ
                            push!(chain, w)
                            v = w
                            k += 1
                            is_end_of_chain = false
                            break
                        end
                    end
                    if is_end_of_chain
                        break
                    end
                end
                push!(selected_chains, chain)
                break
            end
        end
    end
    return selected_chains
end

function solve_master_IP(master_IP::Model, column_pool::Vector{Column}, instance::Instance, bp_status::BP_status, bp_params::BP_params)
    L = instance.max_chain_length
    if bp_params.verbose println("- number of columns in master IP: $(length(master_IP[:y]))\n") end

    # solve the integer master problem
    optimize!(master_IP)

    # update the lower bound and the best solution if an improving one was found
    if has_values(master_IP)
        obj_val = floor(JuMP.objective_value(master_IP) + ϵ)
        if obj_val > bp_status.bp_info.LB + ϵ
            bp_status.bp_info.LB = round(obj_val)
            is_integer = true
            column_val = JuMP.value.(master_IP[:y])
            selected_cycles = Vector{Vector{Int}}()
            selected_chains = Vector{Vector{Int}}()
            for c in 1:length(column_pool)
                if column_val[c] > 1 - ϵ
                    if column_pool[c].is_cycle
                        push!(selected_cycles, (column_pool[c]).vertices)
                    else
                        push!(selected_chains, (column_pool[c]).vertices)
                    end
                end
            end
            bp_status.best_cycles = selected_cycles
            if bp_params.is_pief && instance.max_chain_length >= 1
                bp_status.best_chains = compute_chains_from_pief(master_IP, instance)
            else
                bp_status.best_chains = selected_chains
            end
            if bp_params.verbose printstyled("\nNew incumbent found with value $obj_val found by solving the IP with every column of the pool\n" ; color = :green) end
        end
    end
end
