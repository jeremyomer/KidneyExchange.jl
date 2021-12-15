@enum BloodType A B AB O

struct Instance
    graph::SimpleDiGraph
    vertex_weight::Vector{Float64}
    edge_weight::Matrix{Float64}
    pairs::Vector{Int}
    altruists::Vector{Int}
    nb_pairs::Int
    nb_altruists::Int
    max_cycle_length::Int
    max_chain_length::Int
    is_vertex_weighted::Bool  # true if all weights are actually on the vertices

    # default constructor: assumes weights on vertices; the vertices with zero weight are the altruists
    function Instance(g::SimpleDiGraph, vertex_weight::Array{Float64}, K::Int, L::Int = 0)
        P = [v for v in vertices(g) if vertex_weight[v] != 0.0]
        A = [v for v in vertices(g) if vertex_weight[v] == 0.0]
        edge_weight = zeros(nv(g), nv(g))
        for e in edges(g)
            edge_weight[e.src,e.dst] = vertex_weight[e.dst]
        end

        return new(g, vertex_weight, edge_weight, P, A, length(P), length(A), K, L, true)
    end

    # Parse instance from file
    function Instance(filename::String, K::Int, L::Int = 0)
        inst = string(filename)
        data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-2)], "/"), "data")
        wmd_file = joinpath(data_folder, join([inst, ".wmd"]))
        dat_file = joinpath(data_folder, join([inst, ".dat"]))

        g, edge_weight, is_altruist = read_kep_file(wmd_file, dat_file)
        P = [v for v in vertices(g) if !is_altruist[v]]
        A = [v for v in vertices(g) if is_altruist[v]]
        vertex_weight = zeros(nv(g))
        is_vertex_weighted = true
        for v in P
            vertex_weight[v] = 1.0
            if indegree(g, v) >= 1
                vertex_weight[v] = edge_weight[inneighbors(g,v)[1],v]
                for u in inneighbors(g,v)
                    if edge_weight[u,v] != vertex_weight[v]
                        is_vertex_weighted = false
                        break
                    end
                end
                if !is_vertex_weighted  break   end
            end
        end
        if !is_vertex_weighted
            println("the instance is not vertex weighted!")
        end
        return new(g, vertex_weight, edge_weight, P, A, length(P), length(A), K, L, is_vertex_weighted)
    end
end

mutable struct SubgraphsData
    sources::Vector{Int}
    is_vertex_list::Vector{BitVector}
    d_to_vstar_list::Vector{Vector{Int}}
    d_from_vstar_list::Vector{Vector{Int}}
    nb_copies::Int
    is_arc_list::Vector{BitVector}
    chain_mip::Model

    function SubgraphsData(_sources, _isvertex_list, _d_to::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, 0), _d_from::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, 0), _chain_mip = Model(), _isarc_list::Vector{BitVector}=Vector{BitVector}(undef,0))
        return new(_sources, _isvertex_list, _d_to, _d_from, length(_sources), _isarc_list, _chain_mip)
    end
end
