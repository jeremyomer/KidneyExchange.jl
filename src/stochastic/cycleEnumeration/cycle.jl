# try
#     mutable struct Cycle
#         index::Int
#         type::Int
#         expectation::Float64
#         nodes::Vector{Int}
    
#         function Cycle(_index::Int, instance::Instance, _nodes::Vector{Int}, isCycle::Bool)
#             subgraph = induced_subgraph(instance.graph, _nodes)[1]
#             _type = get_type(subgraph, isCycle)
#             _expectation = get_expectation(instance.pv, instance.pa, _type)
#             return new(_index, _type, _expectation, _nodes)
#         end
#     end
# catch 
# end

# function get_type(
#         subgraph::SimpleDiGraph,
#         isCycle::Bool
# )
#     #Works for K = 2 and K = 3, noclut for K = 4
#     if nv(subgraph) ≤ 3 && isCycle
#         return ne(subgraph) - nv(subgraph) + 1
#     end

#     if nv(subgraph) == 4
        
#     end

#     return 0
# end

# function get_expectation(
#         pv::Float64, #probability of vertice failure
#         pa::Float64, #probability of edge failure
#         type::Int
# )
#     if type == 0
#         return 0
#     end
#     if type == 1
#         return 2*((1-pv)*(1-pa))^2
#     end
    
#     return 3*((1-pv)*(1-pa))^3 + 2*(type-2)((1-pv)*(1-pa))^2

# end

function get_type(graph::SimpleDiGraph, column_list::Vector{Column})
    type_dict = Dict{SimpleDiGraph, Any}()
    isomorph_list = SimpleDiGraph[]
    println(typeof(isomorph_list))
    type_index = 1
    sort!(column_list, by = c -> length(c.vertices))
    
    for (index, column) in enumerate(column_list)

        if index%100000 == 0
             println(index, " cycles got a type !")
        end
        if length(column.vertices) == 2
            column.type = type_index
            continue
        end

        subgraph = induced_subgraph(graph, column.vertices)[1]
        in_search = true
        for isomorph in isomorph_list
            if in_search && has_isomorph(subgraph, isomorph)
                column.type = type_dict[isomorph]
                in_search = false
            end
        end

        if in_search
            push!(isomorph_list, subgraph)
            type_index = type_index + 1
            type_dict[subgraph] = type_index
            column.type = type_index
        end
    end
end

function get_expectation()
    
end

function process_c_list(instance::Instance, column_list::Vector{Column})
    get_type(instance::Instance, column_list)
    get_expectation(instance::Instance, column_list)
end

function path_recursion(
        graph::SimpleDiGraph,
        column_list::Vector{Column},
        node::Int,
        path::Vector{Int},
        path_length::Int
)
    for next_node in neighbors(graph, node)
        if next_node < path[1] || next_node in path[2:path_length]
            continue
        end
        if next_node == path[1]
            column = Column(path, true)
            push!(column_list, column)
            continue
        end
        if path_length < 4
            new_path = [path; next_node]
            path_recursion(graph, column_list, next_node, new_path, path_length+1)
        end
    end

    return
end

function column_enumeration(graph::SimpleDiGraph)
    to = TimerOutput()
    reset_timer!(to)

    #column_list initialization
    column_list = Column[]
    
    #graph = instance.graph
    @timeit to "Enumeration" for node in vertices(graph)
        path = [node]
        path_length = 1
        path_recursion(graph, column_list, node, path, path_length)
    end
    println(length(column_list), "columns enumerated.")
    #Get types and expectation
    
    @timeit to "type" process_c_list(graph, column_list)
    println(to)
    return column_list
end

