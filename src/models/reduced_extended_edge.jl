"""
    build_reduced_extended_edge_mip(instance, subgraphs, params; maxtime = 600)

A compact MIP formulation originally proposed by [Constantino2013](@cite) for the cycles-only variant of the problem.
Here it is adapted to the graph copies based on FVS and chains are considered with position-indexed variables.
"""
function build_reduced_extended_edge_mip(instance::Instance, subgraphs::Graph_copies, params::MIP_params, maxtime::Float64 = 600)
    if params.verbose println("- build the JuMP model") end
    # extract relevant information from the instance and preprocess
    graph = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length
    graph_edges = collect(edges(graph))
    graph_vertices = vertices(graph)
    nb_vertices = nv(graph)
    nb_pairs = instance.nb_pairs
    nb_altruists = instance.nb_altruists
    nb_edges = ne(graph)
    nb_subgraph = subgraphs.nb_copies

    # extract node existence information from subgraphs
    subgraph_vertices = Array{Vector{Int},1}(undef, nb_subgraph)
    for l in 1:nb_subgraph
        subgraph_vertices[l] = Vector{Int}()
        for _v in 1:nb_vertices
            if subgraphs.is_vertex_list[l][_v] == true
                push!(subgraph_vertices[l], _v)
            end
        end
    end
    pair_graph_inds = nb_altruists+1:nb_subgraph

    # extract edge existence information from subgraphs
    subgraph_edges = Vector{Vector{Tuple}}(undef, nb_subgraph)
    for l in 1:nb_subgraph
        subgraph_edges[l] = Vector{Tuple}()
    end
    is_arc = subgraphs.is_arc_list
    arc_ind = Dict{Pair{Int,Int}, Int}()
    for e in 1:nb_edges
        arc_ind[graph_edges[e].src=>graph_edges[e].dst] = e
    end
    for l in 1:nb_subgraph
        for e in 1:nb_edges
            if is_arc[l][e] == true
                push!(subgraph_edges[l], (graph_edges[e].src,graph_edges[e].dst))
            end
        end
    end

    model = create_model(maxtime, params.optimizer, true, params.verbose)

    # create variables x^l_(i,j)
    @variable(model, x[l in pair_graph_inds, (u,v) in subgraph_edges[l]], Bin)

    # create flow variables for altruist subgraphs
    @variable(model, y[u in vertices(graph), v in outneighbors(graph, u), k in 1:L], Bin)

    # create flow conservation constraints for altruist subgraphs
    @constraint(model, chain_flow[u in instance.pairs, k in 1:L-1], sum(y[v,u,k] for v in inneighbors(graph, u)) - sum(y[u,v,k+1] for v in outneighbors(graph, u)) >= 0)
    if L >= 1
        @constraint(model,[u in instance.pairs, v in outneighbors(graph,u)], y[u,v,1] == 0)
        @constraint(model, altruist_flow[u in instance.altruists], sum(y[u,v,1] for v in outneighbors(graph, u)) <= 1)
    end
    @constraint(model, [u in instance.altruists, v in outneighbors(graph,u), k in 2:L], y[u,v,k] == 0)

    # create vertex disjoint constraints for pairs
    @constraint(model, vertex_disjoint[u in instance.pairs], sum(x[l,(u,v)] for l in pair_graph_inds, v in outneighbors(graph,u) if is_arc[l][arc_ind[u=>v]]) + sum(y[v,u,k] for v in inneighbors(graph,u), k in 1:L) <= 1)

    # create flow conservation constraints for each subgraph
    @constraint(model, cycle_flow[l in pair_graph_inds, u in subgraph_vertices[l]], sum(x[l,(v,u)] for v in inneighbors(graph, u) if is_arc[l][arc_ind[v=>u]]) == sum(x[l,(u,v)] for v in outneighbors(graph,u) if is_arc[l][arc_ind[u=>v]]))

    # create cycle_limit constraints
    @constraint(model, cycle_length[l in pair_graph_inds], sum(x[l,e] for e in subgraph_edges[l]) <= K)

    # the flow going out of each vertex of each copy cannot be larger than that going out of the source
    S = subgraphs.sources
    @constraint(model, symmetry_break[l in pair_graph_inds, u in subgraph_vertices[l]], sum(x[l,(u,v)] for v in outneighbors(graph,u) if is_arc[l][arc_ind[u=>v]]) <= sum(x[l,(S[l],v)] for v in outneighbors(graph, S[l]) if is_arc[l][arc_ind[S[l]=>v]]))

    # objective function maximizes the weight of selected cycles
    # - break symmetry by giving a preference to copies with smallest index
    weight = instance.edge_weight
    if L == 0
        if params.symmetry_break
            max_weight = maximum(weight)
            @objective(model, Max, sum((1.0 + l/(max_weight*nb_vertices^2)) * weight[e[1],e[2]] * x[l,e] for l in pair_graph_inds, e in subgraph_edges[l]))
        else
            @objective(model, Max, sum(weight[e[1],e[2]] * x[l,e] for l in pair_graph_inds, e in subgraph_edges[l]))
        end
    else
        if params.symmetry_break
            max_weight = maximum(weight)
            @objective(model, Max, sum((1.0 + l/(max_weight*nb_vertices^2)) * weight[e[1],e[2]] * x[l,e] for l in pair_graph_inds, e in subgraph_edges[l]) + sum(weight[u,v] * y[u,v,k] for u in instance.pairs, v in outneighbors(graph, u), k  in 2:L) + sum(weight[u,v] * y[u,v,1] for u in instance.altruists, v in outneighbors(graph, u)))
        else
            @objective(model, Max, sum(weight[e[1],e[2]] * x[l,e] for l in pair_graph_inds, e in subgraph_edges[l]) + sum(weight[u,v] * y[u,v,k] for u in instance.pairs, v in outneighbors(graph, u), k  in 2:L) + sum(weight[u,v] * y[u,v,1] for u in instance.altruists, v in outneighbors(graph, u)))
        end
    end

    return model
end

function solve_reduced_extended_edge_mip(model::Model, params::MIP_params, instance::Instance, subgraphs::Graph_copies)
    if params.verbose println("- solve the model with $(params.optimizer)") end
    # Local variables
    g = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length

    # Solve the model
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
    if params.optimizer != "HiGHS"
        mip_status.relative_gap = round(10^4*relative_gap(model))/10^4
    end
    # mip_status.node_count = node_count(model)  # the function has errors in last JuMP version
    mip_status.solve_time = solve_time(model)

    # c. get the solution in a readable format
    # get the flow corresponding to chains
    nb_edges = ne(g)
    nb_subgraph = subgraphs.nb_copies
    graph_edges = collect(edges(g))
    subgraph_edges = Vector{Vector{Tuple}}(undef, nb_subgraph)
    for l in 1:nb_subgraph
        subgraph_edges[l] = Vector{Tuple}()
        for e in 1:nb_edges
            arc_key = (graph_edges[e].src,graph_edges[e].dst)
            if subgraphs.is_arc_list[l][e] == true
                push!(subgraph_edges[l], arc_key)
            end
        end
    end
    pair_graph_inds = instance.nb_altruists+1:nb_subgraph

    flow_cycle = Dict{Tuple{Int,Int}, Float64}()
    x = value.(model[:x])
    for e in edges(g)
        flow_cycle[(e.src,e.dst)] = 0.0
    end
    for l in pair_graph_inds
        for e in subgraph_edges[l]
            flow_cycle[e] += x[l,e]
        end
    end

    # get the cycles
    is_covered = falses(nv(g))
    for u in instance.pairs
        if is_covered[u] continue end
        next = u
        cycle = []
        while next != 0
            cur = next
            next = 0
            for v in outneighbors(g, cur)
                if is_covered[v] continue end
                if flow_cycle[(cur,v)] > 1-ϵ
                    next = v
                    push!(cycle, v)
                    is_covered[v] = true
                    break
                end
            end
        end
        if !isempty(cycle)
            is_covered[u] = true
            if cycle[end] != u
                error("not a cycle")
            end
            push!(mip_status.best_cycles, cycle)
        end
    end

    # get the flow corresponding to chains
    if L >= 1
        flow_chain = Dict{Pair{Int,Int}, Float64}()
        for u in vertices(g)
            for v in outneighbors(g, u)
                flow_chain[u=>v] = sum(value(model[:y][u,v,k]) for k in 1:L)
            end
        end

        # get the chains
        for u in instance.altruists
            next = u
            chain = []
            while next != 0
                cur = next
                next = 0
                for v in outneighbors(g, cur)
                    if flow_chain[cur=>v] > 1-ϵ
                        next = v
                        push!(chain, v)
                        break
                    end
                end
            end
            if !isempty(chain)
                pushfirst!(chain, u)
                push!(mip_status.best_chains, chain)
            end
        end
    end


    # Output a completion message if in verbose
    if params.verbose
        printstyled("\n----------------------------------------------------------\n The solution of the EXTENDED EDGE model is complete\n" ; color = :yellow)
        if termination_status(model) == MOI.OPTIMAL
            printstyled("- the solution is optimal\n" ; color = :yellow)
            printstyled("- best solution found: value $(mip_status.objective_value)" ; color = :yellow)
            printstyled("----------------------------------------------------------\n\n" ; color = :yellow)
        elseif termination_status(model) == MOI.TIME_LIMIT
            printstyled("- the time limit is exceeded\n" ; color = :yellow)
            if params.optimizer != "HiGHS"
                printstyled("- best solution found: value $(mip_status.objective_value) with gap $(mip_status.relative_gap) %\n" ; color = :yellow)
            else
                printstyled("- best solution found: value $(mip_status.objective_value)\n" ; color = :yellow)
            end
            printstyled("----------------------------------------------------------\n\n" ; color = :yellow)
        end
    end

    return mip_status
end
