"""
$(SIGNATURES)

Initialize the MIP model with cycle constraint generation for the optimal search of positive cost chains. Only one model is created for every copy to save a great amount of initialization time and memory. The model will then need to be modified for each graph copy to keep only the vertices of the graph and select the right source vertex

# Arguments
* `graph::SimpleDiGraph` : The directed graph with cost on each arc
* `L::Int`: The maximal length of chains
* `optimizer::String`: Name of the MIP sover
* `time_limit::Real`: Time limit of the solver

# Return values
* `model::Model`: Initial JuMP model for the search of a positive chain
"""
function buid_cycle_cuts(
    instance::Instance,
    subgraphs::Graph_copies,
    params::MIP_params,
    maxtime::Float64 = 600,
)
    g = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length
    # set of arcs of the graph as pairs of vertices
    E = [e.src => e.dst for e in edges(g)]

    # add a sink node with one incoming arc from each vertex
    sink = nv(g) + 1
    for v in vertices(g)
        push!(E, v => sink)
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
    for l = 1:nb_subgraphs
        V[l] = Vector{Int}()
        for v = 1:nb_vertices
            if is_vertex[l][v] == true
                push!(V[l], v)
            end
        end
    end
    pair_graph_inds = nb_altruists+1:nb_subgraphs

    # create the JuMP model
    model = create_model(maxtime, params.optimizer, true, params.verbose)

    # one set of cycle flow variables per subgraph
    @variable(
        model,
        cycle_flow[
            l in pair_graph_inds,
            i in V[l],
            j in outneighbors(g, i);
            is_vertex[l][j],
        ],
        Bin
    )

    # create one set of chain flow variables for each altruist vertex
    @variable(model, chain_flow[a in A, u in 1:nb_vertices, v in outneighbors(g, u)], Bin)

    # the chain flow constraints at each vertex
    @constraint(
        model,
        ct_chain_flow[u in A, v in P],
        sum(chain_flow[u, v, w] for w in outneighbors(g, v)) -
        sum(chain_flow[u, w, v] for w in inneighbors(g, v)) <= 0
    )
    @constraint(
        model,
        ct_flow_max[u in A],
        sum(chain_flow[u, u, v] for v in outneighbors(g, u)) <= 1
    )
    @constraint(
        model,
        [u in A, v in A; u != v],
        sum(chain_flow[u, v, w] for w in outneighbors(g, v)) == 0
    )

    # cycle flow conservation constraints for each subgraph
    @constraint(
        model,
        ct_cycle_flow[l in pair_graph_inds, u in V[l]],
        sum(cycle_flow[l, v, u] for v in inneighbors(g, u) if is_vertex[l][v]) ==
        sum(cycle_flow[l, u, v] for v in outneighbors(g, u) if is_vertex[l][v])
    )

    # cycle and chain maximum lengths constraints
    @constraint(
        model,
        ct_cycle_length[l in pair_graph_inds],
        sum(
            cycle_flow[l, u, v] for u in V[l], v in outneighbors(g, u) if is_vertex[l][v]
        ) <= K
    )
    @constraint(
        model,
        ct_chain_length[a in A],
        sum(chain_flow[a, u, v] for u in vertices(g), v in outneighbors(g, u)) <= L
    )

    # vertex disjoint constraints for pairs
    @constraint(
        model,
        ct_vertex_disjoint[u in P],
        sum(
            cycle_flow[l, u, v] for l in pair_graph_inds, v in outneighbors(g, u) if
            is_vertex[l][u] && is_vertex[l][v]
        ) + sum(chain_flow[a, v, u] for v in inneighbors(g, u), a in A) <= 1
    )

    # the cycle flow going out of each vertex of each copy cannot be larger than that going out of the source
    S = subgraphs.sources
    @constraint(
        model,
        [l in pair_graph_inds, u in V[l]],
        sum(cycle_flow[l, u, v] for v in outneighbors(g, u) if is_vertex[l][v]) <=
        sum(cycle_flow[l, S[l], v] for v in outneighbors(g, S[l]) if is_vertex[l][v])
    )
    @constraint(
        model,
        [a in A, v in P],
        sum(chain_flow[a, u, v] for u in inneighbors(g, v)) <=
        sum(chain_flow[a, a, u] for u in outneighbors(g, a))
    )

    # objective function maximizes the weight of selected cycles
    # - break symmetry by giving a preference to copies with smallest index
    weight = instance.edge_weight
    max_weight = maximum(weight)
    if L == 0
        @objective(
            model,
            Max,
            sum(
                (1.0 + l / (max_weight * nb_vertices^2)) *
                weight[u, v] *
                cycle_flow[l, u, v] for
                l in pair_graph_inds, u in V[l], v in outneighbors(g, u) if is_vertex[l][v]
            )
        )
    else
        @objective(
            model,
            Max,
            sum(
                (1.0 + l / (max_weight * nb_vertices^2)) *
                weight[u, v] *
                cycle_flow[l, u, v] for
                l in pair_graph_inds, u in V[l], v in outneighbors(g, u) if is_vertex[l][v]
            ) + sum(
                weight[u, v] * chain_flow[a, u, v] for a in A, u in vertices(g),
                v in outneighbors(g, u)
            )
        )
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
        printstyled("- no solution was found\n"; color = :red)
        return mip_status
    end

    # b. get the main characteristics of the solution
    mip_status.objective_value = floor(objective_value(model))
    if params.optimizer != "HiGHS"
        mip_status.relative_gap = round(10^4 * relative_gap(model)) / 10^4
    end
    # mip_status.node_count = node_count(model)  # the function has errors in last JuMP version
    mip_status.solve_time = solve_time(model)

    return mip_status
end
