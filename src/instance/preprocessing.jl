"""
$(SIGNATURES)

Get the data of each graph copy used in the subproblems of the column generation. The graph copies are charaterized by a source vertex and the set of vertices (optionally set of arcs) kept in the copy.
Each altruist donor is the source of one graph copy. Preprocessing the vertices that can be reached with a given maximum number of arcs does not bring any improvement to the subproblem solution, so we do nothing at this stage.

For the other copies, either each donor/patient pair is the source of one copy, or we compute a feedback vertex set (FVS) that heuristically minimizes the cardinality of the set.  Each vertex of the FVS will then be the source of one copy in the column generation subproblems. In both cases, we then compute distances from and to the source in each subgraph. These will be used in the Bellman-Ford search at each solution of the subproblems.

# Input parameters
* `instance::Instance` : The structure describing the characteristics of an instance including the compatibility graph
* `reduce_arcs::Bool` : A boolean parameter indicating whether or not arcs that are not useful in a copy should be eliminated from it (default value is false)
* `reduce_vertices::Bool` : A boolean parameter indicating whether or not vertices that are not useful in copy should be eliminated from it (default value is true)
* `fvs::Bool` : True if the copies correspond to an FVS, false if there is one copy per donor/patient pair

# Output parameters
* `subgraphs::Graph_copies` : The structure describing the preprocessed data of the graph copies
"""
function preprocess_graph_copies(
    instance::Instance,
    reduce_arcs::Bool = false,
    reduce_vertices::Bool = true,
    fvs::Bool = true,
)
    K = instance.max_cycle_length
    L = instance.max_chain_length
    graph = instance.graph
    nb_vertices = nv(graph)
    sources = Vector{Int}()
    vertex_in_subgraph_list = Vector{BitVector}()
    arc_in_subgraph_list = Vector{BitVector}()
    d_to_vstar_list = Vector{Vector{Int}}()
    d_from_vstar_list = Vector{Vector{Int}}()
    nb_copies = 0
    chain_mip = Model()


    # if there are altruist donors, each one of them will have its own subproblem;
    is_deleted = falses(nv(graph))
    for vstar in instance.altruists
        nb_copies += 1
        push!(sources, vstar)
        push!(d_to_vstar_list, zeros(Int, nv(graph)))
        push!(d_from_vstar_list, zeros(Int, nv(graph)))
        push!(vertex_in_subgraph_list, trues(nv(graph)))
        # do not include the vertices that are farther than L from the vertex, this is useful only for MIP formulations of the chain subproblem
        if reduce_vertices
            bfs(graph, vstar, L, .!is_deleted, d_from_vstar_list[end])
            for v in instance.pairs
                if d_from_vstar_list[end][v] > L
                    vertex_in_subgraph_list[end][v] = false
                end
            end
        end
        for v in instance.altruists
            vertex_in_subgraph_list[end][v] = false
        end
        vertex_in_subgraph_list[end][vstar] = true

        # remove the altruist for what comes next, since it is in no cycle
        is_deleted[vstar] = true
    end

    if fvs
        # for the pairs an FVS is computed to get the source vertices, the score used in the heuristic requires the degrees of the vertices
        in_degrees = zeros(Int, nv(graph))
        out_degrees = zeros(Int, nv(graph))
        fvs_scores = zeros(Int, nv(graph))
        for v in instance.pairs
            in_degrees[v] = indegree(graph, v)
            out_degrees[v] = outdegree(graph, v)
        end

        # do not count the deleted vertices in the degrees
        for u in instance.altruists
            for v in outneighbors(graph, u)
                in_degrees[v] -= 1
            end
        end
    end

    # iteratively compute the vertices of the FVS until there is no more cycle with size less than K
    while true
        # select the vertex with max score to be a feedback vertex: the score is the total degree of the vertex

        if fvs
            # recursively delete the vertices with indegree or outdegree equal to zero, since they cannot be a part of any cycle
            has_deleted = true
            deletion_list = Vector{Int}()
            for v in findall(is_deleted .== false)
                if min(in_degrees[v], out_degrees[v]) <= 0
                    push!(deletion_list, v)
                end
            end
            while !isempty(deletion_list)
                v = pop!(deletion_list)
                if is_deleted[v]
                    continue
                end
                is_deleted[v] = true
                for u in inneighbors(graph, v)
                    if !is_deleted[u]
                        out_degrees[u] -= 1
                        if out_degrees[u] == 0
                            push!(deletion_list, u)
                        end
                    end
                end
                for w in outneighbors(graph, v)
                    if !is_deleted[w]
                        in_degrees[w] -= 1
                        if in_degrees[w] == 0
                            push!(deletion_list, w)
                        end
                    end
                end
            end
        end
        vertices = findall(is_deleted .== false)
        if isempty(vertices)
            break
        end

        # get the source of next graph copy
        if fvs
            # get the vertex with largest total degree
            fvs_scores = in_degrees .+ out_degrees
            vstar = vertices[argmax(fvs_scores[vertices])]
            if fvs_scores[vstar] == 0
                break
            end
        else
            vstar = vertices[1]
        end

        # compute the distance from and to the reference vertex
        push!(sources, vstar)
        push!(d_to_vstar_list, zeros(Int, nv(graph)))
        push!(d_from_vstar_list, zeros(Int, nv(graph)))
        push!(vertex_in_subgraph_list, trues(nv(graph)))
        vertex_in_subgraph_list[end] .= .!is_deleted
        if reduce_vertices
            bfs(graph, vstar, K - 1, .!is_deleted, d_from_vstar_list[end])
            bfs_reverse(graph, vstar, K - 1, .!is_deleted, d_to_vstar_list[end])
            for v in vertices
                if d_from_vstar_list[end][v] + d_to_vstar_list[end][v] > K
                    vertex_in_subgraph_list[end][v] = false
                end
            end
        end

        # delete current source vertex in next copies
        is_deleted[vstar] = true
        if fvs
            for v in inneighbors(graph, vstar)
                out_degrees[v] -= 1
            end
            for v in outneighbors(graph, vstar)
                in_degrees[v] -= 1
            end
        end
    end

    # reduce arcs to get smaller models: this is time-consuming and useful only if there is a real need to get smaller models, so this should be used only when a compact MIP model is solved
    if reduce_arcs
        alledges = collect(edges(graph))
        for s = 1:instance.nb_altruists
            arc_in_subgraph = falses(ne(graph))
            for i = 1:ne(graph)
                u = alledges[i].src
                # keep an arc in the copy only if it can belong to a chain with length at most L originated from the source vertex
                if d_from_vstar_list[s][u] + 1 <= L
                    arc_in_subgraph[i] = true
                end
            end
            push!(arc_in_subgraph_list, arc_in_subgraph)
        end
        for s = instance.nb_altruists+1:length(sources)
            arc_in_subgraph = falses(ne(graph))
            for i = 1:ne(graph)
                u = alledges[i].src
                v = alledges[i].dst
                # keep an arc in the copy only if it can belong to a cycle with length at most K going through the source vertex
                if d_from_vstar_list[s][u] + d_to_vstar_list[s][v] + 1 <= K
                    arc_in_subgraph[i] = true
                end
            end
            push!(arc_in_subgraph_list, arc_in_subgraph)
        end
        return Graph_copies(
            sources,
            vertex_in_subgraph_list,
            d_to_vstar_list,
            d_from_vstar_list,
            chain_mip,
            arc_in_subgraph_list,
        )
    end

    if reduce_vertices
        return Graph_copies(
            sources,
            vertex_in_subgraph_list,
            d_to_vstar_list,
            d_from_vstar_list,
        )
    end
    return Graph_copies(sources, vertex_in_subgraph_list)
end

"""
$(SIGNATURES)

Breadth-first search in the graph where arcs are reversed
# Input parameters
* `g::AbstractGraph{T}` : The graph
* `source::Int` : The source vertex
* `K::Int`: The maximum number of arcs
* `vertex_in_subgraph::Array{Bool}`: Table indicating for each vertex if it is in the subgraph under consideration

# Output parameters
* `dists::Vector{Float64}`: shortest distance from source to vertices
"""
function bfs_reverse(
    g::SimpleDiGraph,
    source::Int,
    K::Int,
    vertex_in_subgraph::BitVector = trues(nv(g)),
    d_to_vstar::Vector{Int} = zeros(Int, nv(g)),
)
    n = nv(g)
    for s = 1:n
        d_to_vstar[s] = n
    end
    d_to_vstar[source] = 0
    to_treat = trues(n)

    to_treat .= vertex_in_subgraph
    to_treat[source] = false
    for v in inneighbors(g, source)
        if to_treat[v]
            d_to_vstar[v] = 1
            to_treat[v] = false
        end
    end
    d = 2
    vertices = findall(to_treat)
    while !isempty(vertices) && d <= K
        for v in vertices
            for u in outneighbors(g, v)
                if d_to_vstar[u] == d - 1
                    d_to_vstar[v] = d
                    to_treat[v] = false
                    break
                end
            end
        end
        if d <= K - 1
            vertices = findall(to_treat)
        end
        d += 1
    end
    return nothing
end

"""
$(SIGNATURES)

Breadth-first search to calculate the shortest path in terms of number of arcs
from source to other vertices of graph g. The search stops when it is at level K
even if there is vertex to be visited. Because we are trying to find cycles of
at most K length, vertices of level >=K are never be considered in cycles.

# Input parameters
* `g::AbstractGraph{T}` : The graph
* `source::Int` : The source vertex
* `K::Int`: The maximum number of arcs
* `vertex_in_subgraph::Array{Bool}`: Table indicating for each vertex if it is in the subgraph under consideration

# Output parameters
* `dists::Vector{Float64}`: shortest distance from source to vertices
"""
function bfs(
    g::SimpleDiGraph,
    source::Int,
    K::Int,
    vertex_in_subgraph::BitVector = trues(nv(g)),
    d_from_vstar::Vector{Int} = zeros(Int, nv(g)),
)
    n = nv(g)
    for s = 1:n
        d_from_vstar[s] = n
    end
    d_from_vstar[source] = 0
    if K == 0
        return nothing
    end

    to_treat = trues(n)
    to_treat .= vertex_in_subgraph
    to_treat[source] = false
    for v in outneighbors(g, source)
        if to_treat[v]
            d_from_vstar[v] = 1
            to_treat[v] = false
        end
    end
    if K == 1
        return nothing
    end

    d = 2
    vertices = findall(to_treat)
    while !isempty(vertices) && d <= K
        for v in vertices
            for u in inneighbors(g, v)
                if d_from_vstar[u] == d - 1
                    d_from_vstar[v] = d
                    to_treat[v] = false
                    break
                end
            end
        end
        if d <= K - 1
            vertices = findall(to_treat)
        end
        d += 1
    end
    return nothing
end
