using Graphs
using Graphs.Experimental



function getAdjacency(n::Int64)
    BaseGraph = zeros((n,n))

    for i in 1:n-1
        BaseGraph[i,i+1] = 1
    end

    BaseGraph[n,1] = 1

    Adjacencies = Matrix{}[BaseGraph]

    for i in 1:2^8-1
        array = digits(i, base = 2, pad = 8)
        graph = copy(BaseGraph)
        graph[1,3] = array[1]
        graph[1,4] = array[2]
        graph[2,1] = array[3]
        graph[2,4] = array[4]
        graph[3,1] = array[5]
        graph[3,2] = array[6]
        graph[4,2] = array[7]
        graph[4,3] = array[8]
        push!(Adjacencies, graph)
    end

    return SimpleDiGraph.(Adjacencies)
end

function get_isomorph_dict()
    Graphs = getAdjacency()
    iso_pool = Vector{SimpleDiGraph}[]
    for g in Graphs
        pass = false 
        for i in iso_pool
            if !pass && has_isomorph(g, i[1])
                push!(i, g)
                pass = true
            end
        end
        if !pass
            push!(iso_pool, [g])
        end
    end
    sort!(adj_iso4, by = g -> ne(g[1]))
    return iso_pool
end

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