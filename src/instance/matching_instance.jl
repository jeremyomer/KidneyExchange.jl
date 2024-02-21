using Graphs

struct PrefLibInstance

    file_path :: String
    file_name :: String
    data_type :: String
    modification_type :: String
    relates_to :: String
    related_files :: String
    title :: String
    description :: String
    publication_date :: String
    modification_date :: String
    num_alternatives :: Int
    alternatives_name :: Dict{Int, String}
    num_voters :: Int
    lines :: Vector{String}

    function PrefLibInstance(filepath)
    
        file_path = filepath
        file_name = splitdir(filepath)[1]
        data_type = splitext(filepath)[2]


        modification_type = ""
        relates_to = ""
        related_files = ""
        title = ""
        description = ""
        publication_date = ""
        modification_date = ""
        num_alternatives = 0
        alternatives_name = Dict()
        num_voters = 0
        lines = readlines(filepath)

        for line in lines

            if startswith(line, "# FILE NAME")
                file_name = split(line,":") |> last |> strip
            elseif startswith(line, "# TITLE")
                title = split(line,":") |> last |> strip
            elseif startswith(line, "# DESCRIPTION")
                description = split(line,":") |> last |> strip
            elseif startswith(line, "# DATA TYPE")
                data_type = split(line,":") |> last |> strip
            elseif startswith(line, "# MODIFICATION TYPE")
                modification_type = split(line,":") |> last |> strip
            elseif startswith(line, "# RELATES TO")
                relates_to = split(line,":") |> last |> strip
            elseif startswith(line, "# RELATED FILES")
                related_files = split(line,":") |> last |> strip
            elseif startswith(line, "# PUBLICATION DATE")
                publication_date = split(line,":") |> last |> strip
            elseif startswith(line, "# MODIFICATION DATE")
                modification_date = split(line,":") |> last |> strip
            elseif startswith(line, "# NUMBER ALTERNATIVES")
                num_alternatives = parse(Int, last(split(line, ":")))
            elseif startswith(line, "# NUMBER VOTERS")
               num_voters = parse(Int, last(split(line, ":")))
            elseif startswith(line, "# ALTERNATIVE NAME")
               m = match(r"# ALTERNATIVE NAME (\d+): (.*)", line)
               alt = parse(Int, m.captures[1])
               alt_name = m.captures[2]
               alternatives_name[alt] = alt_name
           end
        end
        new( file_path, file_name, data_type, modification_type, relates_to,
             related_files, title, description, publication_date, 
             modification_date, num_alternatives, alternatives_name,
             num_voters, lines )

    end
end


struct WeightedDiGraph

    node_mapping :: Dict{Int, Set{Int}}
    weights :: Dict{Tuple{Int,Int}, Float64}

end

neighbours(graph :: WeightedDiGraph, node) = self.node_mapping[node]

outgoing_edges(self :: WeightedDiGraph, node) = Set((node, n, self.weights[(node, n)]) for n in self.node_mapping[node])

function add_node!(self :: WeightedDiGraph, node)
    if node ∉ keys(self.node_mapping)
        self.node_mapping[node] = Set{Int}()
    end
end

function add_edge!(self :: WeightedDiGraph, node1, node2, weight)
    add_node!(self, node1)
    add_node!(self, node2)
    push!(self.node_mapping[node1], node2)
    self.weights[(node1, node2)] = weight
end

function edges(graph :: WeightedDiGraph)

    s = Edge{Int64}[]
    for n1 in keys(graph.node_mapping)
        for n2 in graph.node_mapping[n1]
            push!(s, Edge(n1, n2))
        end
    end

    return s

end

nodes(self :: WeightedDiGraph) = keys(self.node_mapping)

struct MatchingInstance

     instance :: PrefLibInstance
     graph :: WeightedDiGraph

     function MatchingInstance( filepath )

         instance = PrefLibInstance(filepath)
         node_mapping = Dict()
         weights = Dict()
         graph = WeightedDiGraph(node_mapping, weights)
         num_edges = 0

         i = 0
         for line in strip.(instance.lines)
             if startswith(line, "#")
                 if startswith(line, "# NUMBER EDGES")
                     @show num_edges = parse( Int, last(split(line, ":")))
                 end
                 i += 1
             else
                 break
             end
         end

         num_voters = instance.num_alternatives

         for line in strip.(instance.lines[i+1:end])
             vertex1, vertex2, weight = split(line, ",")
             add_edge!(graph, parse(Int, vertex1), parse(Int, vertex2), parse(Float64, weight))
         end

         num_edges = sum(length(edge_set) for edge_set in values(graph.node_mapping))

         new(instance, graph)

     end

end

function read_wmd_file(filepath)

    instance = MatchingInstance(filepath)
    graph = SimpleDiGraph(edges(instance.graph))

    nb_vertices = maximum(max.(keys(instance.graph.weights)...))
    edge_weights = zeros(nb_vertices, nb_vertices)
    for ((n1, n2), w) in instance.graph.weights
        edge_weights[n1, n2] = w
    end

    return graph, edge_weights 

end

read_wmd_file(joinpath("00036-00000001.wmd"))
