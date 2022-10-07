"""
    node_master

Initialization of the restricted master problem

# Arguments
* `instance::Instance`: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains.
* `column_pool::Vector{Column}`: the set of initial columns of the master
* `bp_params::BP_params`: parameters of the branch-and-price
* `time_limit::Float64`: time limit for each solution of the master relaxation

# Return values
* `master::Model`: the model of the restricted master problem
"""
function node_master(instance::Instance, column_pool::Vector{Column}, bp_params::BP_params = BP_params(), time_limit::Float64 = 10000.0)
    # Local variables
    graph = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length

    # enumerate all the triangles of the graph
    triangles = Vector{Tuple{Int,Int,Int}}()
    if K == 2
        for u in instance.pairs
            outlist = outneighbors(graph, u)
            for i in eachindex(outlist) 
                v = outlist[i]
                for j in i+1:length(outlist)
                    w = outlist[j]
                    if !has_edge(graph, v, u) || !has_edge(graph, w, u) continue end
                    if has_edge(graph, v, w) 
                        if has_edge(graph, w, v)
                            if u < v && u < w
                                push!(triangles, (u,v,w))
                            end
                        else
                            push!(triangles, (u,v,w))
                        end
                    elseif has_edge(graph, w, v) 
                        println("not really a triangle")
                        push!(triangles, (u,v,w))
                    end
                end
            end
        end
    end
    end
    error("nb of triangles = $(length(triangles))")


    # Initialize the JuMP model
    master = create_model(time_limit, bp_params.optimizer, false, false)

    # Decision variables
    # - column variables
    @variable(master, y[c in 1:length(column_pool)] >= 0)
    @variable(master, nb_cycles[k in 2:K] >= 0)
    @variable(master, nb_chains[l in 1:L] >= 0)

    # - artificial variable to handle infeasibility after branching
    @variable(master, slack >=0)

    # - in the position-indexed model, we need arc variables with position indexes for the chains
    if bp_params.is_pief
        l_min = 2*ones(nv(graph))
        l_max = L*ones(nv(graph))
        l_min[instance.altruists] .= 1
        l_max[instance.altruists] .= 1
        is_altruist = falses(nv(graph))
        for u in instance.altruists
             is_altruist[u] = true
        end
        altruist_neighbors = Vector{Vector{Int}}(undef, nv(graph))
        pair_neighbors = Vector{Vector{Int}}(undef, nv(graph))
        for v in vertices(graph)
            altruist_neighbors[v] = []
            pair_neighbors[v] = []
        end
        for v in instance.pairs
            for u in inneighbors(graph, v)
                if is_altruist[u]
                    push!(altruist_neighbors[v], u)
                else
                    push!(pair_neighbors[v], u)
                end
            end
        end
        @variable(master, chain_flow[u in vertices(graph), v in outneighbors(graph, u), k in l_min[u]:l_max[u]] >= 0)
    end

    # Constraints
    # - vertex disjoint constraints
    if !bp_params.is_pief || L == 0
        @constraint(master, capacity[v in  vertices(graph)], sum(y[c] for c in 1:length(column_pool) if v in (column_pool[c]).vertices) <= 1)
    else
        @constraint(master, capacity[v in  vertices(graph)], sum(y[c] for c in 1:length(column_pool) if v in (column_pool[c]).vertices) + sum(chain_flow[u,v,k] for u in inneighbors(graph,v), k in l_min[u]:l_max[u]) <= 1)
    end

    # flow constraints of the position-indexed model
    if bp_params.is_pief && L >= 1
        # position-indexed flow conservation constraints
        if L >= 2
            @constraint(master, [u in instance.pairs], sum(chain_flow[v,u,1] for v in altruist_neighbors[u]) - sum(chain_flow[u,v,2] for v in outneighbors(graph, u)) >= 0)
        end
        @constraint(master, [u in instance.pairs, k in 2:L-1], sum(chain_flow[v,u,k] for v in pair_neighbors[u]) - sum(chain_flow[u,v,k+1] for v in outneighbors(graph, u)) >= 0)
        # altruists can start no more than one chain
        @constraint(master, [u in instance.altruists], sum(chain_flow[u,v,1] for v in outneighbors(graph, u)) <= 1)
    end

    # matching polytope cut: for each triangle, the number of arcs selected on the triangle is no more than two
    if K == 2
        # @constraint(master, matching_cuts[(u,v,w) in triangles])
    end

    # - branching constraints on the arcs covered by the columns
    master[:branch_one] = Dict{Pair{Int64, Int64}, ConstraintRef}()
    master[:branch_zero] = Dict{Pair{Int64, Int64}, ConstraintRef}()
    master[:branch_one_vertex] = Dict{Int, ConstraintRef}()
    master[:branch_zero_vertex] = Dict{Int, ConstraintRef}()

    # - branching constraints on the arcs covered by the chain variables of the master when using pief model
    if bp_params.is_pief
        master[:branch_one_pief] = Dict{Pair{Int64, Int64}, ConstraintRef}()
        master[:branch_zero_pief] = Dict{Pair{Int64, Int64}, ConstraintRef}()
    end

    # branching constraint on the total number of cycles and chains depending on their length
    @constraint(master, nb_cycles_per_size[k in 2:K], sum(y[c] for c in 1:length(column_pool) if column_pool[c].is_cycle && (column_pool[c].length == k)) == nb_cycles[k])
    @constraint(master, nb_chains_per_size[l in 1:L], sum(y[c] for c in 1:length(column_pool) if !column_pool[c].is_cycle && (column_pool[c].length == l)) == nb_chains[l])
    @constraint(master, branch_nb_cycles_max[k in 2:K], nb_cycles[k] <= nv(graph))
    @constraint(master, branch_nb_cycles_min[k in 2:K], nb_cycles[k] + slack >= 0)
    @constraint(master, branch_nb_chains_max[l in 1:L], nb_chains[l] <= nv(graph))
    @constraint(master, branch_nb_chains_min[l in 1:L], nb_chains[l] + slack >= 0)

    # objective
    W = instance.edge_weight
    max_cost = maximum(W)
    if bp_params.is_pief && L >= 1
        @objective(master, Max, sum(column_pool[c].weight * y[c] for c in 1:length(column_pool)) + sum(instance.edge_weight[u,v] * chain_flow[u,v,1] for  u in instance.altruists for v in outneighbors(graph, u)) + sum(instance.edge_weight[u,v] * chain_flow[u,v,k] for u in instance.pairs for v in outneighbors(graph, u) for k  in 2:L) - nv(graph) * max_cost * slack)
    else
        @objective(master, Max, sum(column_pool[c].weight * y[c] for c in 1:length(column_pool)) - nv(graph) * max_cost * slack)
    end

    return master
end

"""
    activate_branching_constraints

Activate the all the branching constraints corresponding to a given node of the branch-and-price enumeration tree.
"""
function activate_branching_constraints(master::Model, tree_node::TreeNode, bp_params::BP_params, instance::Instance)
    K = instance.max_cycle_length
    L = instance.max_chain_length

    # - branching constraints on the arcs covered by the columns
    branch_one = master[:branch_one]
    branch_zero = master[:branch_zero]
    for e in tree_node.setone
        set_normalized_rhs(branch_one[e], 1)
    end
    for e in tree_node.setzero
        set_normalized_rhs(branch_zero[e], 0)
    end

    # - branching constraints on the vertices covered by the columns
    branch_one_vertex = master[:branch_one_vertex]
    branch_zero_vertex = master[:branch_zero_vertex]
    for v in tree_node.setone_vertex
        set_normalized_rhs(branch_one_vertex[v], 1)
    end
    for v in tree_node.setzero_vertex
        set_normalized_rhs(branch_zero_vertex[v], 0)
    end

    # - branching constraints on the arcs covered by the chain variables of the master when using pief model
    if bp_params.is_pief
        branch_one_pief = master[:branch_one_pief]
        branch_zero_pief = master[:branch_zero_pief]
        for e in tree_node.setone_pief
            set_normalized_rhs(branch_one_pief[e], 1)
        end
        for e in tree_node.setzero_pief
            set_normalized_rhs(branch_zero_pief[e], 0)
        end
    end

    # - branching constraints on the number of cycles and chains of each size
    for k in 2:K
        set_normalized_rhs(master[:branch_nb_cycles_max][k], tree_node.nb_cycles_max[k-1])
        set_normalized_rhs(master[:branch_nb_cycles_min][k], tree_node.nb_cycles_min[k-1])
    end
    for l in 1:L
        set_normalized_rhs(master[:branch_nb_chains_max ][l], tree_node.nb_chains_max[l])
        set_normalized_rhs(master[:branch_nb_chains_min][l], tree_node.nb_chains_min[l])
    end
end

"""
    deactivate_branching_constraints

Deactivate the all the branching constraints corresponding to a given node of the branch-and-price enumeration tree.
"""
function deactivate_branching_constraints(master::Model, tree_node::TreeNode, bp_params::BP_params, instance::Instance)
    K = instance.max_cycle_length
    L = instance.max_chain_length
    g = instance.graph
    # - branching constraints on the arcs covered by the columns
    branch_one = master[:branch_one]
    branch_zero = master[:branch_zero]
    for e in tree_node.setone
        set_normalized_rhs(branch_one[e], 0)
    end
    for e in tree_node.setzero
        set_normalized_rhs(branch_zero[e], 1)
    end

    # - branching constraints on the vertices covered by the columns
    branch_one_vertex = master[:branch_one_vertex]
    branch_zero_vertex = master[:branch_zero_vertex]
    for v in tree_node.setone_vertex
        set_normalized_rhs(branch_one_vertex[v], 0)
    end
    for v in tree_node.setzero_vertex
        set_normalized_rhs(branch_zero_vertex[v], 1)
    end

    # - branching constraints on the arcs covered by the chain variables of the master when using pief model
    if bp_params.is_pief
        branch_one_pief = master[:branch_one_pief]
        branch_zero_pief = master[:branch_zero_pief]
        for e in tree_node.setone_pief
            set_normalized_rhs(branch_one_pief[e], 0)
        end
        for e in tree_node.setzero_pief
            set_normalized_rhs(branch_zero_pief[e], 1)
        end
    end

    # - branching constraints on the number of cycles and chains of each size
    for k in 2:K
        set_normalized_rhs(master[:branch_nb_cycles_max][k], nv(g))
        set_normalized_rhs(master[:branch_nb_cycles_min][k], 0)
    end
    for l in 1:L
        set_normalized_rhs(master[:branch_nb_chains_max ][l], nv(g))
        set_normalized_rhs(master[:branch_nb_chains_min][l], 0)
    end
end

"""
    initialize_master_IP

Initialization of the master problem with integer columns

# Arguments
* `instance::Instance`: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains.
* `column_pool::Vector{Column}`: the set of initial columns of the master
* `bp_params::BP_params`: parameters of the branch-and-price
* `time_limit::Float64`: time limit for each solution of the master relaxation

# Return values
* `master::Model`: the model of the restricted master problem
"""
function initialize_master_IP(instance::Instance, column_pool::Vector{Column}, bp_params::BP_params = BP_params(), time_limit::Float64 = 10000.0)
    # Local variables
    graph = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length

    # Set time limit properly depending on the master model
    if bp_params.time_limit_master_IP < 0.005 * instance.nb_pairs
        bp_params.time_limit_master_IP = 0.005 * instance.nb_pairs
    end
    if bp_params.is_pief
        bp_params.time_limit_master_IP = time_limit/2.0
    end

    # Initialize the JuMP model
    master = create_model(bp_params.time_limit_master_IP, bp_params.optimizer, true, bp_params.verbose)

    # Decision variables
    # - column variables
    @variable(master, y[c in 1:length(column_pool)], Bin)
    @variable(master, nb_cycles[k in 2:K] >= 0, Int)
    @variable(master, nb_chains[l in 1:L] >= 0, Int)
    # - in the position-indexed model, we need arc variables with position indexes for the chains
    if bp_params.is_pief
        l_min = 2*ones(nv(graph))
        l_max = L*ones(nv(graph))
        l_min[instance.altruists] .= 1
        l_max[instance.altruists] .= 1
        is_altruist = falses(nv(graph))
        for u in instance.altruists
             is_altruist[u] = true
        end
        altruist_neighbors = Vector{Vector{Int}}(undef, nv(graph))
        pair_neighbors = Vector{Vector{Int}}(undef, nv(graph))
        for v in vertices(graph)
            altruist_neighbors[v] = []
            pair_neighbors[v] = []
        end
        for v in instance.pairs
            for u in inneighbors(graph, v)
                if is_altruist[u]
                    push!(altruist_neighbors[v], u)
                else
                    push!(pair_neighbors[v], u)
                end
            end
        end
        @variable(master, chain_flow[u in vertices(graph), v in outneighbors(graph, u), k in l_min[u]:l_max[u]], Bin)
    end

    # Constraints
    # - vertex disjoint constraints
    if !bp_params.is_pief || L == 0
        @constraint(master, capacity[v in  vertices(graph)], sum(y[c] for c in 1:length(column_pool) if v in (column_pool[c]).vertices) <= 1)
    else
        @constraint(master, capacity[v in  vertices(graph)], sum(y[c] for c in 1:length(column_pool) if v in (column_pool[c]).vertices) + sum(chain_flow[u,v,k] for u in inneighbors(graph,v), k in l_min[u]:l_max[u]) <= 1)
    end

    # flow constraints of the position-indexed model
    if bp_params.is_pief && L >= 1
        # position-indexed flow conservation constraints
        if L >= 2
            @constraint(master, [u in instance.pairs], sum(chain_flow[v,u,1] for v in altruist_neighbors[u]) - sum(chain_flow[u,v,2] for v in outneighbors(graph, u)) >= 0)
        end
        @constraint(master, [u in instance.pairs, k in 2:L-1], sum(chain_flow[v,u,k] for v in pair_neighbors[u]) - sum(chain_flow[u,v,k+1] for v in outneighbors(graph, u)) >= 0)
        # altruists can start no more than one chain
        @constraint(master, [u in instance.altruists], sum(chain_flow[u,v,1] for v in outneighbors(graph, u)) <= 1)
    end

    # branching constraint on the total number of cycles and chains depending on their length
    @constraint(master, nb_cycles_per_size[k in 2:K], sum(y[c] for c in 1:length(column_pool) if column_pool[c].is_cycle && (column_pool[c].length == k)) == nb_cycles[k])
    @constraint(master, nb_chains_per_size[l in 1:L], sum(y[c] for c in 1:length(column_pool) if !column_pool[c].is_cycle && (column_pool[c].length == l)) == nb_chains[l])
    @constraint(master, branch_nb_cycles_max[k in 2:K], nb_cycles[k] <= nv(graph))
    @constraint(master, branch_nb_cycles_min[k in 2:K], nb_cycles[k] >= 0)
    @constraint(master, branch_nb_chains_max[l in 1:L], nb_chains[l] <= nv(graph))
    @constraint(master, branch_nb_chains_min[l in 1:L], nb_chains[l] >= 0)

    #objective
    W = instance.edge_weight
    nb_altruists = instance.nb_altruists
    nb_pairs = instance.nb_pairs
    if bp_params.is_pief && L >= 1
        # break symmetry in the model by privileging cycles over chains in the objective function (without loosing the proof of optimality of the solution
        max_cost = maximum(W)
        @objective(master, Max,  (1.0 + 1/(nb_pairs*max_cost)) * sum(column_pool[c].weight * y[c] for c in 1:length(column_pool)) + sum(W[u,v] * chain_flow[u,v,1] for u in instance.altruists for v in outneighbors(graph, u)) + sum(W[u,v] * chain_flow[u,v,k] for u in instance.pairs for v in outneighbors(graph, u) for k  in 2:L))
    else
        @objective(master, Max, sum(column_pool[c].weight * y[c] for c in 1:length(column_pool)))
    end

    return master
end



"""
    activate_branching_constraints

Activate some branching constraints in the master IP.
"""
function activate_branching_constraints_IP(master_IP::Model, tree_node::TreeNode, instance::Instance)
    K = instance.max_cycle_length
    L = instance.max_chain_length

    # - branching constraints on the number of cycles and chains of each size
    for k in 2:K
        set_normalized_rhs(master_IP[:branch_nb_cycles_max][k], tree_node.nb_cycles_max[k-1])
        set_normalized_rhs(master_IP[:branch_nb_cycles_min][k], tree_node.nb_cycles_min[k-1])
    end
    for l in 1:L
        set_normalized_rhs(master_IP[:branch_nb_chains_max ][l], tree_node.nb_chains_max[l])
        set_normalized_rhs(master_IP[:branch_nb_chains_min][l], tree_node.nb_chains_min[l])
    end
end

"""
    deactivate_branching_constraints

Deactivate some branching constraints in the master IP.
"""
function deactivate_branching_constraints_IP(master_IP::Model, instance::Instance)
    K = instance.max_cycle_length
    L = instance.max_chain_length
    g = instance.graph

    # - branching constraints on the number of cycles and chains of each size
    for k in 2:K
        set_normalized_rhs(master_IP[:branch_nb_cycles_max][k], nv(g))
        set_normalized_rhs(master_IP[:branch_nb_cycles_min][k], 0)
    end
    for l in 1:L
        set_normalized_rhs(master_IP[:branch_nb_chains_max ][l], nv(g))
        set_normalized_rhs(master_IP[:branch_nb_chains_min][l], 0)
    end
end
