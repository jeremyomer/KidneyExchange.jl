"""
    calculate_arc_cost

Calculates in parallel the cost of arcs of the original graph

#Input parameters
* `original_graph::MyDirectedGraph` : The original graph
* `arc_cost::SharedMatrix{Float64}` : The cost of arcs to be calculated
* `λ:: Vector{Float64}` : The dual values of the vertex disjoint constraints
* `δ_arc::Array{Array{Int,1},1}`: The arcs to be branched, if arc (i,j) is branched, δ_arc[i] will contain j
* `δ_val::Array{Array{Float64,1},1}`: dual value of the branching constraints associated with arcs to be branched
"""
# JO: in my view, this function is not useful, we could directly use the vectors of dual variables in the subproblem, and it would make the code easier to read
function calculate_arc_cost(instance::Instance, arc_cost::Matrix{Float64}, λ::Vector{Float64}, δ_one::Dict{Pair{Int64, Int64}, Float64}, δ_zero::Dict{Pair{Int64, Int64}, Float64}, δ_one_vertex::Dict{Int, Float64}, δ_zero_vertex::Dict{Int, Float64})
    # initialize all arc costs with arc weights and retrieve destination dual cost for all vertices; also retrieve source dual for arcs outgoing from altruists
    for u in vertices(instance.graph)
        for v in outneighbors(instance.graph, u)
            # add a small perturbation to guarantee that the only zero costs are where there is no edge
            arc_cost[u,v] = instance.edge_weight[u,v] - λ[v] + rand()*1.0e-6
        end
    end
    for u in instance.altruists
        for v in outneighbors(instance.graph, u)
            arc_cost[u,v] -= λ[u]
        end
    end

    # finally, retrieve the duals of all branching constraints
    for arc in keys(δ_one)
        arc_cost[arc[1],arc[2]] -= δ_one[arc]
    end
    for arc in keys(δ_zero)
        arc_cost[arc[1],arc[2]] -= δ_zero[arc]
    end
    for v in keys(δ_one_vertex)
        for u in inneighbors(instance.graph, v)
            arc_cost[u,v] -= δ_one_vertex[v]
        end
    end
    for v in keys(δ_zero_vertex)
        for u in inneighbors(instance.graph, v)
            arc_cost[u,v] -= δ_zero_vertex[v]
        end
    end
    return maximum(arc_cost)
end

"""
    traverse_preds

Inner function for Bellman_Ford_cycle_search to return path of exactly (n-1) of length ending at v from the table of predecessor

#Input parameters
* `v::Int` : The ending vertex
* `pred::Array{Array{Tuple{Int,Int},1},1}`: table of predecessor containing tuples of predecessor and legnth of path to arrive the vertex
* `n::Int`: the length of the path

#Output Parameters
* `c=Array{Int,1}`: path of exactly (n-1) of length ending at v
"""

function traverse_preds(v::Int, pred::Vector{Vector{Int}}, n::Int)
    c = Vector{Int}()
    push!(c, v)
    current = v
    for i = n:-1:1
        u = pred[i][current]
        pushfirst!(c, u)
        current = u
    end
    return c
end

"""
    Bellman_Ford_cycle_search

Bellman-Ford style search for one positive cost cycle

#Input parameters
* `graph::SimpleDiGraph` : The directed graph with cost on each arc
* `arc_cost::Matrix{Float64}`:
* `source::Int` : The local vertex index from which starts the search
* `K::Int` : The maximal length of cycles
* pred[k][v] contains u if u is a predecessor of v in a path of length k from source

#Output Parameters
* `cycle::Vector{Int}`: the positive cycle found, [] if none
"""
function Bellman_Ford_cycle_search(graph::SimpleDiGraph, vertex_cost::Vector{Float64}, source::Int, K::Int, is_covered::BitVector, pred::Vector{Vector{Int}}, d:: Vector{Float64})
    if K <= 1 return [] end

    # initialize distances
    tmp_d = Vector{Float64}(undef, nv(graph))
    for u in 1:nv(graph)
        d[u] = -Inf
        pred[1][u]=0
    end
    d[source] = vertex_cost[source]  # we can count the return arc cost right away

    # Bellman-Ford search
    # a. treat the first iteration differently as we know it propagates the costs along arcs outgoing from the source
    is_to_source = falses(nv(graph))
    for v in inneighbors(graph, source)
        is_to_source[v] = true
    end
    @inbounds for v in outneighbors(graph, source)
        if is_covered[v]
            continue
        end
        # always update the distance to neighbors since they are not reached yet
        d[v] = d[source]
        pred[1][v] = source
        # always try and take a return arc to the starting vertex in order to stop the algorithm prematurely
        if is_to_source[v] && (d[v] + vertex_cost[v] > ϵ)
            return [source;v]
        end
    end
    if K == 2 return [] end

    # b. following iterations
    max_cost = maximum(vertex_cost)
    is_ignored = BitArray(undef, nv(graph))
    cycle = []
    for k in 2:K-2
        is_ignored .= is_covered
        pred[k] .= zeros(Int, nv(graph))
        nb_edges_left = K-k
        vtx_list = findall(pred[k-1] .!= 0)
        tmp_d[vtx_list] .= d[vtx_list] .+ vertex_cost[vtx_list]
        sort_ind = sortperm(tmp_d[vtx_list] ; rev=true)
        vtx_list = vtx_list[sort_ind]
        @inbounds for u in vtx_list
            # the cost of an arc cannot be larger than max_cost, so we can only hope to get a nb_edges_left*max_cost extra cost with the remaining arcs of the cycle
            if tmp_d[u] + nb_edges_left * max_cost < ϵ
                 break
            end
            d_with_tol = tmp_d[u] - ϵ
            @inbounds for v in outneighbors(graph,u)
                if is_ignored[v]
                    continue
                end
                is_ignored[v] = true
                # update if the distance through u is larger than current distance
                if d[v] < d_with_tol
                    pred[k][v] = u
                    d[v] = tmp_d[u]
                    # always try and take a return arc to the starting vertex when the distance to a vertex is updated in order to stop the algorithm prematurely
                    if is_to_source[v] && (d[v] + vertex_cost[v] > ϵ)
                        cycle = traverse_preds(v, pred, k)
                        # subcycles need to be checked only for cycles with more than 4 edges, we never update a distance if it uses a subcycle
                        if k >= 3
                            for i in 1:length(cycle)-1
                                for j in i+1:length(cycle)
                                    if cycle[i] == cycle[j]
                                        cycle = cycle[i:j-1]
                                    end
                                end
                            end
                        end
                        return cycle
                    end
                end
            end
        end
    end

    # c. last iteration
    #  ignore the vertices that are not predecessor to the source
    is_ignored .= .!is_to_source .| is_covered

    # no need for tmp array here, we will not overwrite d in the inside loop
    vtx_list = findall(pred[K-2] .!= 0)
    d[vtx_list] .+= vertex_cost[vtx_list]
    sort_ind = sortperm(d[vtx_list] ; rev=true)
    vtx_list = vtx_list[sort_ind]
    @inbounds for u in vtx_list
        # the cost of an arc cannot be larger than max_cost, so we can only hope to get an extra max_cost with the remaining arc of the cycle
        if d[u] + max_cost < ϵ
             break
        end
        @inbounds for v in outneighbors(graph,u)
            if is_ignored[v]
                continue
            end
            is_ignored[v] = true
            # we know that there is a return arc to the source since we kept only the in-neighbors of the source
            if d[u] + vertex_cost[v] > ϵ
                pred[K-1][v] = u
                cycle = traverse_preds(v, pred, K-1)
                # subcycles need to be checked only for cycles with more than 4 edges
                if K >= 4
                    for i in 1:length(cycle)-1
                        for j in i+1:length(cycle)
                            if cycle[i] == cycle[j]
                                cycle = cycle[i:j-1]
                            end
                        end
                    end
                end
                return cycle
            end
        end
    end

    return []
end

function Bellman_Ford_cycle_search(graph::SimpleDiGraph, arc_cost::Array{Float64,2}, arc_cost_trans::Matrix{Float64}, max_cost::Float64, source::Int, K::Int, is_covered::BitVector, pred::Vector{Vector{Int}}, d:: Vector{Float64})
    if K <= 1 return [] end
    # initialize distances
    for u in 1:nv(graph)
        d[u] = -Inf
    end
    pred[1] .= zeros(Int, nv(graph))
    d[source] = 0.0

    # Bellman-Ford search
    # a. treat the first iteration differently as we know it propagates the costs along arcs outgoing from the source
    # the dual cost of a cycle is the sum of dual costs of the vertices of the cycle
    for v in outneighbors(graph, source)
        if is_covered[v]
            continue
        end
        # always update the distance to neighbors since they are not reached yet
        d[v] = d[source] + arc_cost_trans[v, source]
        pred[1][v] = source
        # always try and take a return arc to the starting vertex in order to stop the algorithm prematurely (arc_cost is exactly zero only if there is no edge)
        if arc_cost[v, source] != 0.0 && (d[v] + arc_cost[v,source] > ϵ)
            return [source;v]
        end
    end
    if K == 2 return [] end

    # b. following iterations
    cycle = []
    tmp_d = Vector{Float64}(undef, nv(graph))
    for k in 2:K-1
        tmp_d .= d
        pred[k] .= zeros(Int, nv(graph))
        nb_edges_left = K-k+1  # the return cost to starting vertex has not been counted yet
        for u in shuffle!(findall(pred[k-1] .!= 0))
            # the cost of an arc cannot be larger than max_cost, so we can only hope to get a nb_edges_left*max_cost extra cost with the remaining arcs of the cycle
            if tmp_d[u] + nb_edges_left * max_cost < ϵ
                 continue
            end
            if k >= 3
                cycle = traverse_preds(u, pred, k-1)
            end
            for v in outneighbors(graph,u)
                if is_covered[v]
                    continue
                end
                # update if the distance through u is larger than current distance
                cost_to_v = tmp_d[u] + arc_cost_trans[v,u]
                if cost_to_v > d[v] + ϵ
                    if k >= 3
                        if v ∈ cycle
                            continue
                        end
                    end
                    d[v] = cost_to_v
                    pred[k][v] = u
                    # always try and take a return arc to the starting vertex when the distance to a vertex is updated in order to stop the algorithm prematurely
                    if arc_cost[v, source] != 0.0 && (d[v] + arc_cost[v,source] > ϵ)
                        # at this stage there cannot be any subcycle
                        return traverse_preds(v, pred, k)
                    end
                end
            end
        end
    end
    return []
end


function Bellman_Ford_chain_search(graph::SimpleDiGraph,    vertex_cost::Vector{Float64}, source::Int, L::Int, K::Int, is_covered::BitVector, pred::Vector{Vector{Int}}, d:: Vector{Float64}, verbose::Bool = true)
    if L == 0   return []    end
    # initialize distances to vertices
    tmp_d = Vector{Float64}(undef, nv(graph))
    for i in 1:nv(graph)
        d[i] = -Inf
        pred[1][i] = 0
    end
    d[source] = vertex_cost[source]  # count the dual cost of the altruist source vertex

    # Bellman-Ford search
    # a. treat the first iteration differently as we know it propagates the costs along arcs outgoing from the source
    for v in outneighbors(graph, source)
        if is_covered[v]
            continue
        end
        # always update the distance to neighbors since they are not reached yet
        d[v] = d[source]
        pred[1][v] = source
        # check if the chain improves the incumbent
        if d[v] + vertex_cost[v] > ϵ
            return [source;v]
        end
    end
    if L == 1 return []  end

    # b. in the next iterations, we need to be cautious with positive subcycles
    max_cost = maximum(vertex_cost)
    is_ignored = falses(nv(graph))
    for l in 2:L
        nb_edges_left = L-l+1
        is_ignored .= is_covered
        pred[l] .= zeros(Int, nv(graph))
        vtx_list = findall(pred[l-1] .!= 0)
        for u in vtx_list
            tmp_d[u] = d[u] + vertex_cost[u]
        end
        sort_ind = sortperm(tmp_d[vtx_list]; rev=true)
        vtx_list = vtx_list[sort_ind]
        for u in vtx_list
            # the cost of an arc cannot be larger than max_cost, so we can only hope to get a  nb_edges_left*max_cost extra cost with the remaining arcs of the cycle
            if tmp_d[u] + nb_edges_left * max_cost < ϵ
                 break
            end
            if l >= 3
                chain_to_u = traverse_preds(u, pred, l-1)
            end
            for v in outneighbors(graph,u)
                if is_ignored[v]
                    continue
                end
                is_ignored[v] = true
                # update if the distance through u is larger than current distance
                if tmp_d[u] > d[v] + ϵ
                    # do not record and propagate a chain if not elementary
                    if l >= 3
                        if v ∉ chain_to_u[2:end-1]
                            d[v] = tmp_d[u]
                            pred[l][v] = u
                        else
                            # consider the vertex in other chains since this one was wrongly considered
                            is_ignored[v] = false
                            continue
                        end
                    else
                        d[v] = tmp_d[u]
                        pred[l][v] = u
                    end
                    # check if the chain improves the incumbent
                    if d[v] + vertex_cost[v] > ϵ
                        return traverse_preds(v, pred, l)
                    end
                end
            end
        end
    end
    return []
end

function Bellman_Ford_chain_search_optimality(graph::SimpleDiGraph,    vertex_cost::Vector{Float64}, source::Int, L::Int, K::Int, is_covered::BitVector, pred::Vector{Vector{Int}}, d:: Vector{Float64}, verbose::Bool = true)
    if L == 0   return [], false    end
    # initialize distances to vertices
    tmp_d = Vector{Float64}(undef, nv(graph))
    for i in 1:nv(graph)
        d[i] = -Inf
        pred[1][i] = 0
    end
    d[source] = vertex_cost[source]  # count the dual cost of the altruist source vertex
    is_positive_cycle = false  # true if a positive cycle was found at some point

    # Bellman-Ford search
    # a. treat the first iteration differently as we know it propagates the costs along arcs outgoing from the source
    for v in outneighbors(graph, source)
        if is_covered[v]
            continue
        end
        # always update the distance to neighbors since they are not reached yet
        d[v] = d[source]
        pred[1][v] = source
        # check if the chain improves the incumbent
        if d[v] + vertex_cost[v] > ϵ
            return [source;v], false
        end
    end
    if L == 1 return [], false  end

    # b. in the next iterations, we need to be cautious with positive subcycles
    max_cost = maximum(vertex_cost)
    is_ignored = falses(nv(graph))
    for l in 2:L
        nb_edges_left = L-l+1
        is_ignored .= is_covered
        pred[l] .= zeros(Int, nv(graph))
        vtx_list = findall(pred[l-1] .!= 0)
        for u in vtx_list
            tmp_d[u] = d[u] + vertex_cost[u]
        end
        sort_ind = sortperm(tmp_d[vtx_list];rev=true)
        vtx_list = vtx_list[sort_ind]

        for u in vtx_list
            # the cost of an arc cannot be larger than max_cost, so we can only hope to get a  nb_edges_left*max_cost extra cost with the remaining arcs of the cycle
            if tmp_d[u] + nb_edges_left*max_cost < ϵ
                break
            end
            for v in outneighbors(graph,u)
                if is_ignored[v]
                    continue
                end
                is_ignored[v] = true
                # update if the distance through u is larger than current distance
                if tmp_d[u] > d[v] + ϵ
                    d[v] = tmp_d[u]
                    pred[l][v] = u
                    # check if the chain improves the incumbent
                    if d[v] + vertex_cost[v] > ϵ
                        if l >= K+2
                            chain_to_v = traverse_preds(v, pred, l)
                            if allunique(chain_to_v)
                                return chain_to_v, false
                            else
                                is_positive_cycle = true
                            end
                        else
                            return traverse_preds(v, pred, l), false
                        end
                    end
                end
            end
        end
    end
    return [], is_positive_cycle
end

function Bellman_Ford_chain_search(graph::SimpleDiGraph,    arc_cost_trans::Matrix{Float64}, max_cost::Float64, source::Int, L::Int, K::Int, is_covered::BitVector, pred::Vector{Vector{Int}}, d:: Vector{Float64}, verbose::Bool = true)
    if L == 0
        return []
    end
    # initialize distances to vertices
    tmp_d = Vector{Float64}(undef, nv(graph))
    for i in 1:nv(graph)
        d[i] = -Inf
        pred[1][i] = 0
    end
    d[source] = 0.0

    # Bellman-Ford search
    # a. treat the first iteration differently as we know it propagates the costs along arcs outgoing from the source
    # the dual cost of a cycle is the sum of dual costs of the vertices of the cycle
    for v in outneighbors(graph, source)
        if is_covered[v]
            continue
        end
        # always update the distance to neighbors since they are not reached yet
        d[v] = d[source] + arc_cost_trans[v,source]
        pred[1][v] = source
        # check if the chain improves the incumbent
        if d[v] > ϵ
            return [source;v]
        end
    end

    # b. in the last iterations, we need to be cautious with positive subcycles
    chain_to_u = []
    for l in 2:L
        nb_edges_left = L-l+1
        pred[l] .= zeros(Int, nv(graph))
        tmp_d .= d
        for u in shuffle!(findall(pred[l-1] .!= 0))
            # the cost of an arc cannot be larger than max_cost, so we can only hope to get a nb_edges_left*max_cost extra cost with the remaining arcs of the cycle
            if tmp_d[u] + nb_edges_left * max_cost < ϵ
                continue
            end
            if l >= 3
                chain_to_u = traverse_preds(u, pred, l-1)
            end
            for v in outneighbors(graph,u)
                if is_covered[v]
                    continue
                end
                # update if the distance through u is larger than current distance
                cost_to_v = tmp_d[u] + arc_cost_trans[v,u]
                if cost_to_v > d[v] + ϵ
                    if l >= 3
                        if v ∈ chain_to_u[2:end-1]
                            # do not record and propagate a chain if not elementary
                            continue
                        end
                    end
                    d[v] = cost_to_v
                    pred[l][v] = u

                    # check if the chain improves the incumbent
                    if d[v] > ϵ
                        return traverse_preds(v, pred, l)
                    end
                end
            end
        end
    end

    return []
end

function Bellman_Ford_chain_search_optimality(graph::SimpleDiGraph,    arc_cost_trans::Matrix{Float64}, max_cost::Float64, source::Int, L::Int, K::Int, is_covered::BitVector, pred::Vector{Vector{Int}}, d:: Vector{Float64}, verbose::Bool = true)
    if L == 0
        return [], false
    end
    # initialize distances to vertices
    tmp_d = Vector{Float64}(undef, nv(graph))
    for i in 1:nv(graph)
        d[i] = -Inf
        pred[1][i] = 0
    end
    d[source] = 0.0
    is_positive_cycle = false  # true if a positive cycle was found at some point

    # Bellman-Ford search
    # a. treat the first iteration differently as we know it propagates the costs along arcs outgoing from the source
    # the dual cost of a cycle is the sum of dual costs of the vertices of the cycle
    for v in outneighbors(graph, source)
        if is_covered[v]
            continue
        end
        # always update the distance to neighbors since they are not reached yet
        d[v] = d[source] + arc_cost_trans[v, source]
        pred[1][v] = source
        # check if the chain improves the incumbent
        if d[v] > ϵ
            return [source;v], false
        end
    end

    # b. in the last iterations, we need to be cautious with positive subcycles
    for l in 2:L
        nb_edges_left = L-l+1 # the return cost to starting vertex has not been counted yet
        pred[l] .= zeros(Int, nv(graph))
        tmp_d .= d
        for u in findall(pred[l-1] .!= 0)  # useless to shuffle in optimality search
            # the cost of an arc cannot be larger than max_cost, so we can only hope to get a nb_edges_left*max_cost extra cost with the remaining arcs of the cycle
            if tmp_d[u] + nb_edges_left * max_cost < ϵ
                continue
            end
            for v in outneighbors(graph,u)
                if is_covered[v]
                    continue
                end
                # update if the distance through u is larger than current distance
                cost_to_v = tmp_d[u] + arc_cost_trans[v,u]
                if cost_to_v > d[v] + ϵ
                    # propagate any chain when trying to prove optimality
                    d[v] = cost_to_v
                    pred[l][v] = u
                    # check if the chain improves the incumbent
                    if d[v] > ϵ
                        # return the positive chain if does not contain any subcycle
                        if l >= K+2
                            chain_to_v = traverse_preds(v, pred, l)
                            if allunique(chain_to_v)
                                return chain_to_v, false
                            else
                                is_positive_cycle = true
                            end
                        else
                            return traverse_preds(v, pred, l), false
                        end
                    end
                end
            end
        end
    end
    return [], is_positive_cycle
end


"""
    create_chain_mip

Initialize the MIP model with cycle constraint generation for the optimal search of positive cost chains. Only one model is created for every copy to save a great amount of initialization time and memory. The model will then need to be modified for each graph copy to keep only the vertices of the graph and select the right source vertex

#Input parameters
* `graph::SimpleDiGraph` : The directed graph with cost on each arc
* `L::Int`: The maximal length of chains
* `optimizer::String`: Name of the MIP sover
* `time_limit::Real`: Time limit of the solver
* `gurobi_env`: Gurobi environment if Guroib is used (avoids many messages from the solver)
ln("   . solution found: ", arcs, ", objective value: ", objective_value(mip))
                println("   . cost of other vertices = ", [vertex_cost[e[2]] for e in arcs])
                println("   . the mip search did not find any positive cost chain")
            end
#Output Parameters
* `mip::Model`: Initial JuMP model for the search of a positive chain
"""
function create_chain_mip(graph::SimpleDiGraph, L::Int, optimizer::String, time_limit::Float64)
    # set of arcs of the graph copy represented as pairs of vertices
    E = [e.src=>e.dst for e in edges(graph)]

    # add a sink node with one incoming arc from each vertex
    sink = nv(graph) + 1
    for v in vertices(graph)
        push!(E, v=>sink)
    end

    # create the JuMP model and add one variable per arc
    mip = create_model(time_limit, optimizer, true, false)
    @variable(mip, is_arc[e in E], Bin)

    # add the flow constraints at each vertex
    @constraint(mip, ct_flow_conservation[v in vertices(graph)], sum(is_arc[v=>w] for w in outneighbors(graph, v)) + is_arc[v=>sink] - sum(is_arc[u=>v] for u in inneighbors(graph, v)) == 0)

    # add constraints that will be used to ignore the vertices that are not in the considered graph copy
    @constraint(mip, ct_flow_max[v in vertices(graph)], sum(is_arc[u=>v] for u in inneighbors(graph, v)) <= 1)

    # flow constraint at sink
    @constraint(mip, sum(is_arc[v=>sink] for v in vertices(graph)) == 1)

    # maximum length of the chain
    @constraint(mip, sum(is_arc[e] for e in E) <= L + 1)

    return mip
end

"""
    MIP_chain_search

IP model with subtour elimination constraints. The constraints are the  generalized cutset inequalities (GCS) and they are added dynamically in a row generation algorithm.
Refer for instance to the following reference for a presentation of the GCS.
Taccari, Leonardo. « Integer Programming Formulations for the Elementary Shortest Path Problem ». European Journal of Operational Research 252, nᵒ 1 (2016).

#Input parameters
* `mip::Model`: The JuMP model initialized with the flow conservation constraints and the bound on the length of the chain
* `graph::SimpleDiGraph` : The directed graph with cost on each arc
* `source::Int`: Index of the source vertex
* `is_vertex::BitVector`: For each vertex, indicates if it is in the considered subgraph
* `arc_cost::Vector{Float64}`: Matrix of reduced costs of every arc

#Output Parameters
* `is_positivie_chain::Bool`: True if a positive chain was found
* `chain::Vector{Int}:` Positive chain that was found
"""
function MIP_chain_search(mip::Model, graph::SimpleDiGraph, source::Int, is_vertex::BitVector, arc_cost::Matrix{Float64}, verbose::Bool = true)
    sink = nv(graph) + 1
    if verbose
        println("- solve mip chain search")
        println("   . $source is the souce and $sink is the artificial sink" )
    end

    # set the objective of the MIP
    is_arc = mip[:is_arc]
    E = [e.src=>e.dst for e in edges(graph) if (is_vertex[e.src] && is_vertex[e.dst])]
    @objective(mip, Max, sum(arc_cost[e[1],e[2]] * is_arc[e[1]=>e[2]] for e in E))

    # modify the constraints of the MIP to consider only the vertices in the graph copy and set the right source
    set_normalized_rhs(mip[:ct_flow_conservation][source], 1)
    set_normalized_rhs(mip[:ct_flow_max][source], 0)
    for v in vertices(graph)
        if !is_vertex[v]
            set_normalized_rhs(mip[:ct_flow_max][v], 0)
        end
    end

    is_positive_chain = false
    chain = Vector{Int}()
    while true
        optimize!(mip)
        if (termination_status(mip) != MOI.OPTIMAL)
            error("The mip chain search model was not solved to optimality")
        end

        # get the graph induced by the support of the solution
        is_arc_val = value.(is_arc)
        arcs = E[findall([is_arc_val[e] for e in E] .>= 1 - ϵ)]
        vertex_list = [e[1] for e in arcs]  # list of covered vertices
        successor = Dict{Int,Int}()  # successor of each covered vertex
        for e in arcs
            successor[e[1]] = e[2]
        end
        if verbose println("   . solution found: ", arcs, ", objective value: ", objective_value(mip)) end
        if verbose println("   . arc from source to sink: $(is_arc_val[source=>sink])") end

        # stop if the optimal value is non-positive
        if objective_value(mip) <= ϵ
            if verbose println("   . the mip search did not find any positive cost chain") end
            break
        end

        # build the positive chain or detect a positive cycle
        # - initialize the chain with arc leaving the source
        chain = Vector{Int}()
        subtour = copy(vertex_list)
        ind_source = findfirst(vertex_list .== source)
        if ind_source != nothing
            push!(chain, source)
            push!(chain, successor[source])

            # - build the path iteratively from the source; no cycle can be found because of the flow constraints
            while haskey(successor, chain[end])
                 push!(chain, successor[chain[end]])
            end

             # at this stage, we have found a chain starting from source and ending at sink: check its cost
             cost =  sum(arc_cost[chain[i],chain[i+1]] for i in 1:length(chain)-2)
             if cost > ϵ
                 is_positive_chain = true
                 break
             end
         else
             if verbose println("    . source goes directly to sink in solution") end
         end

         # if the cost is non-positive, there is a positive cycle: delete the chain from the solution and add a cut to forbid what's left
         subtour = Vector{Int}()
         is_in_subtour = trues(length(vertex_list))
         for v in chain
             is_in_subtour[v] = false
         end
         subtour = vertex_list[findall(is_in_subtour .== true)]
         is_in_subtour = falses(nv(graph))
         for v in subtour
             is_in_subtour[v] = true
         end

         if verbose println("add constraints to eliminate subtour: ", subtour) end
         delta_subtour = [u=>v for u in subtour for v in outneighbors(graph, u) if !is_in_subtour[v]]
         delta_subtour = [delta_subtour ; [u=>sink for u in subtour]]
         for u in subtour
             @constraint(mip, sum(is_arc[e] for e in delta_subtour) >= sum(is_arc[u=>v] for v in outneighbors(graph, u)))
         end
    end

    # reset the constraints of the MIP
    set_normalized_rhs(mip[:ct_flow_conservation][source], 0)
    set_normalized_rhs(mip[:ct_flow_max][source], 1)
    for v in vertices(graph)
        if !is_vertex[v]
            set_normalized_rhs(mip[:ct_flow_max][v], 1)
        end
    end

    return is_positive_chain, chain
end
