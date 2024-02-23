import Dates
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

        graph = SimpleWeightedDiGraph(sources, destinations, weights)

        new(preflib_instance, graph)

    end

end

export read_wmd_file

"""
$(SIGNATURES)

Contruct a KEP graph from a `.wmd` file. The format of these files
can be found on the [PrefLib website](https://www.preflib.org/dataset/00036). 
The extra `.dat` file mentionned in `.wmd` file metadata must be provided. It contains
individual information on the patient and donor of each pair such
as blood type. The .wmd file describes the edges of the KEP graph.
For example, in the first instance of the benchmark, the .dat looks
like this :

```
	Pair,Patient,Donor,Wife-P?,%Pra,Out-Deg,Altruist
	1,A,B,0,0.05,2,0
	2,O,A,0,0.05,4,0
	3,A,B,0,0.05,2,0
	4,O,A,1,0.5875,3,0
	5,B,AB,0,0.05,0,0
	6,B,A,0,0.45,3,0
	7,O,A,0,0.45,4,0
	8,B,A,0,0.45,3,0
	9,O,A,0,0.45,4,0
	10,O,O,1,0.2875,11,0
	11,A,AB,0,0.05,0,0
	12,A,A,1,0.5875,3,0
	13,O,O,1,0.925,10,0
	14,O,A,0,0.45,4,0
	15,O,A,0,0.05,3,0
	16,O,B,1,0.2875,3,0
```

And a preview of the .wmd file (including only a subset of the arcs) looks like this: 

```
	# FILE NAME: 00036-00000001.wmd
	# TITLE: Kidney Matching - 16 with 0
	# DESCRIPTION:
	# DATA TYPE: wmd
	# MODIFICATION TYPE: synthetic
	# RELATES TO:
	# RELATED FILES: 00036-00000001.dat
	# PUBLICATION DATE: 2013-08-17
	# MODIFICATION DATE: 2022-09-16
	# NUMBER ALTERNATIVES: 16
	# NUMBER EDGES: 59
	# ALTERNATIVE NAME 1: Pair 1
	# ALTERNATIVE NAME 2: Pair 2
	# ALTERNATIVE NAME 3: Pair 3
	# ALTERNATIVE NAME 4: Pair 4
	...
	1, 5, 1.0
	1, 6, 1.0
	2, 1, 1.0
	2, 3, 1.0
	2, 11, 1.0
	2, 12, 1.0
	3, 5, 1.0
	3, 8, 1.0
	4, 1, 1.0
	4, 3, 1.0
	...
```

# Parameters

* `wmd_file::String` : Absolute path of the `.wmd` file.
"""
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

"""
$(SIGNATURES)

Write a `.wmd` and a `.dat` files to store the input KEP graph in the same form as in Preflib.

"""
function write_preflib_file(
    kep_graph::SimpleDiGraph,
    edge_weight::Matrix{Float64},
    donorBT::Vector{Blood_type},
    patientBT::Vector{Blood_type},
    wifeP::BitArray,
    patientPRA::Vector{Float64},
    is_altruist::BitArray,
    file_name::String
)

    dat_filename = file_name * ".dat"
    wmd_filename = file_name * ".wmd"
    io_wmd = open(wmd_filename, "w")

    num_vertices = nv(kep_graph)
    num_edges = ne(kep_graph)

    title = "KEP"
    description = "KEP"
    num_alternatives = length(is_altruist)

    println(io_wmd, "# FILE NAME: $wmd_filename")
    println(io_wmd, "# TITLE: $title")
    println(io_wmd, "# DESCRIPTION:$description")
    println(io_wmd, "# DATA TYPE: wmd")
    println(io_wmd, "# MODIFICATION TYPE: synthetic")
    println(io_wmd, "# RELATES TO:")
    println(io_wmd, "# RELATED FILES: $dat_filename")
    println(io_wmd, "# PUBLICATION DATE: $(Dates.today())")
    println(io_wmd, "# MODIFICATION DATE: $(Dates.today())")
    println(io_wmd, "# NUMBER ALTERNATIVES: $(num_alternatives)")
    println(io_wmd, "# NUMBER EDGES: $num_edges")

    P = findall(is_altruist .== false)
    A = findall(is_altruist .== true)

    for (i,v) in enumerate(P)
        println(io_wmd, "# ALTERNATIVE NAME $v : Pair $v")
    end
    for v in A
        println(io_wmd, "# ALTERNATIVE NAME $v : Altruist $v")
    end

    for e in edges(kep_graph)
        println(io_wmd, "$(e.src),$(e.dst), $(edge_weight[e.src,e.dst])")
    end

    close(io_wmd)

    io_dat = open(dat_filename, "w")
    println(io_dat, "Pair,Patient,Donor,Wife-P?,%Pra,Out-Deg,Altruist")

    for v in P
        println(io_dat, "$v,$(patientBT[v]),$(donorBT[v]),$(Int(wifeP[v])),$(patientPRA[v]),$(outdegree(kep_graph,v)),0")
    end
    for v in A
     	println(io_dat, "$v,O,$(donorBT[v]),0,0.0,$(outdegree(kep_graph,v)),1")
    end
    close(io_dat)

end
