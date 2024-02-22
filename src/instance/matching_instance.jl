using Graphs, SimpleWeightedGraphs

struct PrefLibInstance

    file_path::String
    file_name::String
    data_type::String
    modification_type::String
    relates_to::String
    related_files::String
    title::String
    description::String
    publication_date::String
    modification_date::String
    num_alternatives::Int
    alternatives_name::Dict{Int,String}
    num_voters::Int
    lines::Vector{String}

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
        lines = strip.(readlines(filepath))

        for line in lines

            if startswith(line, "# FILE NAME")
                file_name = split(line, ":") |> last |> strip
            elseif startswith(line, "# TITLE")
                title = split(line, ":") |> last |> strip
            elseif startswith(line, "# DESCRIPTION")
                description = split(line, ":") |> last |> strip
            elseif startswith(line, "# DATA TYPE")
                data_type = split(line, ":") |> last |> strip
            elseif startswith(line, "# MODIFICATION TYPE")
                modification_type = split(line, ":") |> last |> strip
            elseif startswith(line, "# RELATES TO")
                relates_to = split(line, ":") |> last |> strip
            elseif startswith(line, "# RELATED FILES")
                related_files = split(line, ":") |> last |> strip
            elseif startswith(line, "# PUBLICATION DATE")
                publication_date = split(line, ":") |> last |> strip
            elseif startswith(line, "# MODIFICATION DATE")
                modification_date = split(line, ":") |> last |> strip
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
        new(
            file_path,
            file_name,
            data_type,
            modification_type,
            relates_to,
            related_files,
            title,
            description,
            publication_date,
            modification_date,
            num_alternatives,
            alternatives_name,
            num_voters,
            lines,
        )

    end
end


struct MatchingInstance

    preflib_instance::PrefLibInstance
    graph::SimpleWeightedDiGraph{Int,Float64}

    function MatchingInstance(filepath)

        preflib_instance = PrefLibInstance(filepath)

        i = 0
        for line in preflib_instance.lines
            if startswith(line, "#")
                if startswith(line, "# NUMBER EDGES")
                    num_edges = parse(Int, last(split(line, ":")))
                end
                i += 1
            else
                break
            end
        end

        sources = Int[]
        destinations = Int[]
        weights = Float64[]

        for line in preflib_instance.lines[i+1:end]
            vertex1, vertex2, weight = split(line, ",")
            push!(sources, parse(Int, vertex1))
            push!(destinations, parse(Int, vertex2))
            push!(weights, parse(Float64, weight))
        end

        graph = SimpleWeightedDiGraph( sources, destinations, weights)

        new(preflib_instance, graph)

    end

end

export read_wmd_file

function read_wmd_file(filepath)

    @assert last(splitext(filepath)) == ".wmd" "not a wmd file"
    instance = MatchingInstance(filepath)
    graph = SimpleDiGraph([Edge(e.src, e.dst) for e in edges(instance.graph)])

    weights = collect(instance.graph.weights)

    related_files = instance.preflib_instance.related_files
    if isfile(related_files)
        lines, header = readdlm(related_files, '\n'; header = true)
        is_altruist = parse.(Bool, last.(lines)) |> vec
    end

    return graph, weights, is_altruist

end
