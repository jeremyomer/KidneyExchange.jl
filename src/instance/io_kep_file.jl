export read_kep_file

"""
$(SIGNATURES)

Contruct a KEP graph from a `.wmd` and a `.dat` input files. The format of these files can be found on the PrefLib website (see [here](https://www.preflib.org/dataset/00036) for more details).
The `.dat` file provides individual information on the patient and donor of each pair such as blood type. The .wmd file describes the edges of the KEP graph. For example, in the first instance of the benchmark, the .dat looks like this :

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
* `dat_file::String` : Absolute path of the `.dat` file.
"""
function read_kep_file(wmd_file::AbstractString, dat_file::AbstractString)

    wmd_file_name = first(splitext(basename(wmd_file)))
    dat_file_name = first(splitext(basename(dat_file)))

    @assert wmd_file_name == dat_file_name ".wmd and .dat files do not correspond to the same dataset."

	@assert isfile(abspath(wmd_file)) "$(abspath(wmd_file)): file not found."
	@assert isfile(abspath(dat_file)) ".dat file not found."

	# Get the number of vertices and edges from the first line of wmd file
    wmd_io = open(wmd_file, "r")
	splitted_line = split(readline(wmd_io), ',')
	nb_vertices = parse(Int, splitted_line[1])
	nb_edges = parse(Int, splitted_line[2])

    # Extract meta information from the .dat file
    file = readdlm(dat_file, '\n')
	is_altruist = falses(nb_vertices)
	ind = 1
    for line in file[2:end]
        splitted_line = split(line, ',')
        if Bool(parse(Int, splitted_line[7]))
			is_altruist[ind] = true
		end
		ind += 1
    end

    # Extract the graph structure from the .wmd file
    # skip next nb_vertices lines, which are redundant with the data contained in the .dat file
    for i in 1:nb_vertices
        readline(wmd_io)
    end

    # read the set of edges
	in_list = Vector{Vector{Int}}(undef, nb_vertices)
	out_list = Vector{Vector{Int}}(undef, nb_vertices)
	for u in 1:nb_vertices
		in_list[u] = Vector{Int}()
		out_list[u] = Vector{Int}()
	end
	ne = 0
	edge_weight = zeros(nb_vertices, nb_vertices)
    while !eof(wmd_io)
        splitted_line = split(readline(wmd_io), ',')
        # /!\ Pairs are numbered from 0 in the second part of the file
        src = parse(Int, splitted_line[1]) + 1
        dst = parse(Int, splitted_line[2]) + 1
        weight = parse(Float64, splitted_line[3])

        # do not add an edge that has an altruist as destination or that has a zero weight
		if !is_altruist[dst] && weight > 0.0
			ne += 1
			push!(out_list[src], dst)
			push!(in_list[dst], src)
			edge_weight[src, dst] = weight
		end
    end
	for u in 1:nb_vertices
		sort!(out_list[u])
		sort!(in_list[u])
	end
    return SimpleDiGraph(ne, out_list, in_list), edge_weight, is_altruist
end

"""
$(SIGNATURES)

Write a `.wmd` and a `.dat` files to store the input KEP graph in the same form as in Preflib. Store the files in ./data directory.

# Parameters
* `kep_graph::SimpleDiGraph`: KEP graph to write
* `edge_weight::Matrix{Float64}`: weight of each eadge of the KEP graph
* `donorBT::Vector{Blood_type}`: blood type of the donor of each pair
* `patientBT::Vector{Blood_type}`: blood type of the patient of each pair
* `wifeP::BitArray`: for each pair true iff the donor and patient are married
* `patientPRA::Vector{Float64}`: PRA of the patient of each pair
* `is_altruist::BitArray`: for each pair, true if the donor is an altruist donor
* `file_name::String`: Name of the files (before the extensions)
"""
function write_kep_file(kep_graph::SimpleDiGraph, edge_weight::Matrix{Float64}, donorBT::Vector{Blood_type}, patientBT::Vector{Blood_type}, wifeP::BitArray, patientPRA::Vector{Float64}, is_altruist::BitArray, file_name::String)
	data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-2)], "/"), "data/")
	io_dat = open(data_folder * file_name * ".dat", "w")
	io_wmd = open(data_folder * file_name * ".wmd", "w")

	# Initialize the files
	write(io_wmd, "$(nv(kep_graph)),$(ne(kep_graph))\n")
	write(io_dat, "Pair,Patient,Donor,Wife-P?,%Pra,Out-Deg,Altruist\n")

	# Write the data relative to vertices
	P = findall(is_altruist .== false)
	A = findall(is_altruist .== true)
	println("Write $file_name")
	println("- write vertices")
	for v in P
		write(io_wmd, "$v,Pair $v\n")
		write(io_dat, "$v,$(patientBT[v]),$(donorBT[v]),$(Int(wifeP[v])),$(patientPRA[v]),$(outdegree(kep_graph,v)),0\n")
	end
	for v in A
		write(io_wmd, "$v,Altruist $v\n")
		write(io_dat, "$v,O,$(donorBT[v]),0,0.0,$(outdegree(kep_graph,v)),1\n")
	end
	println("- write edges")
	for e in edges(kep_graph)
		write(io_wmd, "$(e.src - 1),$(e.dst - 1), $(edge_weight[e.src,e.dst])\n")
	end
	println("$file_name is written")
	close(io_dat)
	close(io_wmd)
end
