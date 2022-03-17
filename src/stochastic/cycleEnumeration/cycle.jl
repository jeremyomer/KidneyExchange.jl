try
    mutable struct Cycle
        index::Int
        type::Int
        expectation::Float64
        nodes::Vector{Int}
    
        function Cycle(_index::Int, instance::Instance, _nodes::Vector{Int}, isCycle::Bool)
            subgraph = induced_subgraph(instance.graph, _nodes)[1]
            _type = get_type(subgraph, isCycle)
            _expectation = get_expectation(instance.pv, instance.pa, _type)
            return new(_index, _type, _expectation, _nodes)
        end
    end
catch 
end

function get_type(
        subgraph::SimpleDiGraph,
        isCycle::Bool
)
    #Works for K = 2 and K = 3, not for K = 4
    if nv(subgraph) ≤ 3 && isCycle
        return ne(subgraph) - nv(subgraph) + 1
    end

    return 0
end

#Need to consider K = 4 and different weights and probabilities for vertices and edges
function get_expectation(
        pv::Float64, #probability of vertice failure
        pa::Float64, #probability of edge failure
        type::Int
)
    if type == 0
        return 0
    end
    if type == 1
        return 2*((1-pv)*(1-pa))^2
    end
    
    return 3*((1-pv)*(1-pa))^3 + 2*(type-2)((1-pv)*(1-pa))^2

end

function cycle_recursion(
        instance::Instance,
        cycles::Vector{Vector{Cycle}},
        cycle_index::Int,
        node::Int,
        path::Vector{Int},
        path_length::Int
)
    for next_node in neighbors(instance.graph, node)
        if next_node < path[1]
            continue
        end
        if next_node == path[1]
            cycle_index = cycle_index + 1
            cycle = Cycle(cycle_index, instance, path, true)
            push!(cycles[cycle.type], cycle)
            continue
        end
        if path_length < K
            new_path = [path; next_node]
            cycle_index = cycle_recursion(instance, cycles, cycle_index, next_node, new_path, path_length+1)
        end
    end

    return cycle_index
end

function cycle_enumeration(instance::Instance)
    cycles = Vector{Vector{Cycle}}(undef, 5) 
    for i in 1:5
        cycle[i] = Vector{Cycle}[]
    end
    cycle_index = 0
        for node in vertices(instance.graph)
            path = [node]
            path_length = 1
            cycle_index = cycle_recursion(instance, cycles, cycle_index, node, path, path_length)
        end

    return cycles, cycle_index
end

