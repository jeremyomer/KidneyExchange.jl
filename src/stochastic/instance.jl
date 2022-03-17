struct ce_Instance
    graph::SimpleDiGraph
    vertex_weight::Vector{Float64}
    edge_weight::Matrix{Float64}
    pairs::Vector{Int}
    altruists::Vector{Int}
    nb_pairs::Int
    nb_altruists::Int
    max_cycle_length::Int
    max_chain_length::Int
    pa::Float64
    pv::Float64
    is_vertex_weighted::Bool

    function ce_Instance(g::SimpleDiGraph, vertex_weight::Array{Float64}, K::Int, L::Int = 0, pa::Float64, pv::Float64)
        P = [v for v in vertices(g) if vertex_weight[v] != 0.0]
        A = [v for v in vertices(g) if vertex_weight[v] == 0.0]
        edge_weight = zeros(nv(g), nv(g))
        for e in edges(g)
            edge_weight[e.src,e.dst] = vertex_weight[e.dst]
        end

        return new(g, vertex_weight, edge_weight, P, A, length(P), length(A), K, L, pa, pv, true)
    end

    function ce_Instance(filename::String, K::Int, L::Int = 0, pa::Float64, pv::Float64)
        inst = string(filename)
		data_folder = "data" # joinpath(pkgdir(KidneyExchange), "data")
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
        return new(g, vertex_weight, edge_weight, P, A, length(P), length(A), K, L, pa, pv, is_vertex_weighted)
    end
    
end