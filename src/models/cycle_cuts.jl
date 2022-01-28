"""
    create_chain_mip

Initialize the MIP model with cycle constraint generation for the optimal search of positive cost chains. Only one model is created for every copy to save a great amount of initialization time and memory. The model will then need to be modified for each graph copy to keep only the vertices of the graph and select the right source vertex

# Arguments
* `graph::SimpleDiGraph` : The directed graph with cost on each arc
* `L::Int`: The maximal length of chains
* `optimizer::String`: Name of the MIP sover
* `time_limit::Real`: Time limit of the solver

# Return values
* `model::Model`: Initial JuMP model for the search of a positive chain
"""
function buid_cycle_cuts(instance::Instance, subgraphs::Graph_copies, params::MIP_params, maxtime::Float64 = 600)
    g = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length
    # set of arcs of the graph as pairs of vertices
    E = [e.src=>e.dst for e in edges(g)]

    # add a sink node with one incoming arc from each vertex
    sink = nv(g) + 1
    for v in vertices(g)
        push!(E, v=>sink)
    end


    graph_vertices = vertices(g)
    nb_vertices = nv(g)
    P = instance.pairs
    A = instance.altruists
    nb_subgraphs = subgraphs.nb_copies
    nb_altruists = instance.nb_altruists


    # extract node existance information from subgraphs
    V = Vector{Vector{Int}}(undef, nb_subgraphs)
    is_vertex = subgraphs.is_vertex_list
    for l in 1:nb_subgraphs
        V[l] = Vector{Int}()
        for v in 1:nb_vertices
            if is_vertex[l][v] == true
                push!(V[l], v)
            end
        end
    end
    pair_graph_inds = nb_altruists+1:nb_subgraphs

    # create the JuMP model
    model = create_model(maxtime, params.optimizer, true, params.verbose)

    # one set of cycle flow variables per subgraph
    @variable(model, cycle_flow[l in pair_graph_inds, i in V[l], j in outneighbors(g, i); is_vertex[l][j]], Bin)

    # create one set of chain flow variables for each altruist vertex
    @variable(model, chain_flow[a in A, u in 1:nb_vertices, v in outneighbors(g, u)], Bin)

    # the chain flow constraints at each vertex
    @constraint(model, ct_chain_flow[u in A, v in P], sum(chain_flow[u,v,w] for w in outneighbors(g, v)) - sum(chain_flow[u,w,v] for w in inneighbors(g, v)) <= 0)
    @constraint(model, ct_flow_max[u in A], sum(chain_flow[u,u,v] for v in outneighbors(g,u)) <= 1)
    @constraint(model, [u in A, v in A; u!= v], sum(chain_flow[u,v,w] for w in outneighbors(g,v)) == 0)

    # cycle flow conservation constraints for each subgraph
    @constraint(model, ct_cycle_flow[l in pair_graph_inds, u in V[l]], sum(cycle_flow[l,v,u] for v in inneighbors(g, u) if is_vertex[l][v]) == sum(cycle_flow[l,u,v] for v in outneighbors(g,u) if is_vertex[l][v]))

    # cycle and chain maximum lengths constraints
    @constraint(model, ct_cycle_length[l in pair_graph_inds], sum(cycle_flow[l,u,v] for u in V[l], v in outneighbors(g, u) if is_vertex[l][v]) <= K)
    @constraint(model, ct_chain_length[a in A], sum(chain_flow[a,u,v] for u in vertices(g), v in outneighbors(g, u)) <= L)

    # vertex disjoint constraints for pairs
    @constraint(model, ct_vertex_disjoint[u in P], sum(cycle_flow[l,u,v] for l in pair_graph_inds, v in outneighbors(g,u) if  is_vertex[l][u] && is_vertex[l][v]) + sum(chain_flow[a,v,u] for v in inneighbors(g,u), a in A) <= 1)

    # the cycle flow going out of each vertex of each copy cannot be larger than that going out of the source
    S = subgraphs.sources
    @constraint(model, [l in pair_graph_inds, u in V[l]], sum(cycle_flow[l,u,v] for v in outneighbors(g,u) if is_vertex[l][v]) <= sum(cycle_flow[l,S[l],v] for v in outneighbors(g, S[l]) if is_vertex[l][v]))
    @constraint(model, [a in A, v in P], sum(chain_flow[a,u,v] for u in inneighbors(g,v)) <= sum(chain_flow[a,a,u] for u in outneighbors(g, a)))

    # objective function maximizes the weight of selected cycles
    # - break symmetry by giving a preference to copies with smallest index
    weight = instance.edge_weight
    max_weight = maximum(weight)
    if L == 0
        @objective(model, Max, sum((1.0 + l/(max_weight*nb_vertices^2)) * weight[u,v] * cycle_flow[l,u,v] for l in pair_graph_inds, u in V[l], v in outneighbors(g,u) if is_vertex[l][v]))
    else
        @objective(model, Max, sum((1.0 + l/(max_weight*nb_vertices^2)) * weight[u,v] * cycle_flow[l,u,v] for l in pair_graph_inds, u in V[l], v in outneighbors(g,u) if is_vertex[l][v]) + sum(weight[u,v] * chain_flow[a,u,v] for a in A, u in vertices(g), v in outneighbors(g,u)))
    end

    optimize!(model)

    # Build the MIP status structure
    mip_status = Solution_status()

    # a. first get the solution status
    if termination_status(model) == MOI.OPTIMAL
        mip_status.status = "OPTIMAL"
    elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
        mip_status.status = "TIME_LIMIT"
    else
        mip_status.status = "NOT_FEASIBLE"
        printstyled("- no solution was found\n" ; color = :red)
        return mip_status
    end

    # b. get the main characteristics of the solution
    mip_status.objective_value = floor(objective_value(model))
    mip_status.relative_gap = mip_status.relative_gap = round(10^4*relative_gap(model))/10^4
    # mip_status.node_count = node_count(model)  # the function has errors in last JuMP version
    mip_status.solve_time = solve_time(model)

    return mip_status
end

"""
    MIP_chain_search

MIP model with cycle constraint generation for the optimal search of positive cost chains

# Arguments
* `graph::SimpleDiGraph` : The directed graph with cost on each arc
* `source::Int`: Index of the source vertex
* `is_vertex::BitVector`: For each vertex, indicates if it is in the considered subgraph
* `vertex_cost::Vector{Float64}`: Dual costs of the vertices

# Return values
* `is_positivie_chain::Bool`: True if a positive chain was found
* `chain::Vector{Int}:` Positive chain that was found
"""
# function MIP_chain_search(model::Model, graph::SimpleDiGraph, source::Int, is_vertex::BitVector, vertex_cost::Vector{Float64}, verbose::Bool = true)
#     if verbose
#         println("- solve model chain search")
#         println("   . cost of source vertex = ", vertex_cost[source])
#     end
#     # set the objective of the MIP
#     E = [e.src=>e.dst for e in edges(graph) if (is_vertex[e.src] && is_vertex[e.dst])]
#     chain_flow = model[:chain_flow]
#     @objective(model, Max, -vertex_cost[source] + sum((1 - vertex_cost[e[2]]) * chain_flow[e[1]=>e[2]] for e in E))
#
#     # modify the constraints of the MIP to consider only the vertices in the graph copy and set the right source
#     set_normalized_rhs(model[:ct_flow_conservation][source], 1)
#     set_normalized_rhs(model[:ct_flow_max][source], 0)
#     for v in vertices(graph)
#         if !is_vertex[v]
#             set_normalized_rhs(model[:ct_flow_max][v], 0)
#         end
#     end
#
#     is_positive_chain = false
#     chain = Vector{Int}()
#     while true
#         optimize!(model)
#         if (termination_status(model) != MOI.OPTIMAL)
#             error("The model chain search model was not solved to optimality")
#         end
#         is_arc_val = value.(chain_flow)
#         arcs = E[findall([is_arc_val[e] for e in E] .>= 1 - ϵ)]
#         if objective_value(model) <= ϵ
#             if verbose
#                 println("   . solution found: ", arcs, ", objective value: ", objective_value(model))
#                 println("   . cost of other vertices = ", [vertex_cost[e[2]] for e in arcs])
#                 println("   . the model search did not find any positive cost chain")
#             end
#             break
#         end
#
#         # build the positive chain or detect a positivie cycle
#         # - initialize the chain with arc leaving the source
#         chain = Vector{Int}()
#         chain_arcs = Vector{Pair{Int,Int}}()
#         is_in_chain = falses(length(arcs))
#         for i in 1:length(arcs)
#             e = arcs[i]
#             if e[1] == source
#                 push!(chain, e[1])
#                 push!(chain, e[2])
#                 push!(chain_arcs, e)
#                 is_in_chain[i] = true
#                 break
#             end
#         end
#         if isempty(chain)
#             error("The source should be covered by the selected arcs")
#         end
#
#         # - build the path iteratively from the source; no cycle can be found because of the flow constraints
#         ind_next = findfirst([e[1] for e in arcs] .== chain[end])
#         if verbose println("arcs = ", arcs) end
#         while ind_next != nothing
#              next_arc = arcs[ind_next]
#              is_in_chain[ind_next] = true
#              push!(chain_arcs, next_arc)
#              push!(chain, next_arc[2])
#              ind_next = findfirst([e[1] for e in arcs] .== chain[end])
#          end
#
#          # at this stage, we have found a chain: check its cost
#          cost =  length(chain)-1 - sum(vertex_cost[v] for v in chain)
#          if cost > ϵ
#              is_positive_chain = true
#              break
#          end
#
#          # if the cost is non-positive, there is a positive cycle: delete the chain from the solution and add a cut to forbid what's left
#          deleteat!(arcs, is_in_chain)
#          if verbose println("add a cut to forbid cycle(s) : ", arcs) end
#          @constraint(model, sum(chain_flow[e] for e in arcs) <= length(arcs)-1)
#     end
#
#     # reset the constraints of the MIP
#     set_normalized_rhs(model[:ct_flow_conservation][source], 0)
#     set_normalized_rhs(model[:ct_flow_max][source], 1)
#     for v in vertices(graph)
#         if !is_vertex[v]
#             set_normalized_rhs(model[:ct_flow_max][v], 1)
#         end
#     end
#
#     return is_positive_chain, chain
# end
