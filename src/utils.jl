using Cbc
using DelimitedFiles
using GLPK
using Graphs
using JuMP
using Random
using TimerOutputs

include("instance/instance.jl")
include("instance/io_kep_file.jl")
include("instance/preprocessing.jl")
include("typedef.jl")
include("createModel.jl")

"""
    print_and_check_solution

Function that prints the main characteristics of a solution and check that it satisfies all teh constraints

# Input parameters
* `cycles::Vector{Vector{Int}}`: list of cycles in the solution to display
* `chains::Vector{Vector{Int}}`: list of chains in the solution to display
* `instance::Instance`: KEP instance whose solution it is
* `verbose::Bool`: true if the characteristics of the solution are printed, false if only the verification need be done

# Output parameters: None
"""
function print_and_check_solution(cycles::Vector{Vector{Int}}, chains::Vector{Vector{Int}}, instance::Instance, verbose::Bool = true)
    cycles_lengths = [length(c) for c in cycles]
    chains_lengths = [length(c) for c in chains]

    # Display the number of chains and cycles for each path length up to the maximum allowed
    nb_vertices = 0
    max_cycles_lengths = 0
    if !isempty(cycles_lengths)
        max_cycles_lengths = maximum(cycles_lengths)
    end
    if verbose println("Numbers of cycles per cycle length") end
    for k in 2:max_cycles_lengths
        nb_cycles = length(findall(cycles_lengths .== k))
        if verbose  && nb_cycles >= 1
            println("- k = $k: $nb_cycles cycles")
         end
        nb_vertices += nb_cycles * k
    end
    max_chains_lengths = 0
    if !isempty(chains_lengths)
        max_chains_lengths = maximum(chains_lengths)
    end
    if verbose println("In total, $nb_vertices pairs are covered by cycles\n") end
    if instance.nb_altruists > 0
        if verbose println("Numbers of chains per chain length") end
        nb_vertices = 0
        for l in 2:max_chains_lengths
            nb_chains = length(findall(chains_lengths .== l))
            if verbose && nb_chains >= 1
                println("- l = $(l-1): $nb_chains chains")
            end
            nb_vertices += nb_chains * (l-1)
        end
        if verbose println("In total, $nb_vertices pairs are covered by chains\n") end
    end

    # check that the cycles and chains of the solution use only edges of the KEP graph and compute the objective value of the solution at the same time
    g = instance.graph
    vertex_coverage = zeros(Int, nv(g))
    obj_val = 0.0
    for c in cycles
        vertex_coverage[c[1]] += 1
        for i in 2:length(c)
            if !has_edge(g, c[i-1], c[i])
                error("An edge of the solution is not in the KEP graph")
            end
            obj_val += instance.edge_weight[c[i-1], c[i]]
            vertex_coverage[c[i]] += 1
        end
        if !has_edge(g, c[end], c[1])
            error("An edge of the solution is not in the KEP graph")
        end
        obj_val += instance.edge_weight[c[end], c[1]]
    end
    for c in chains
        vertex_coverage[c[1]] += 1
        for i in 2:length(c)
            if !has_edge(g, c[i-1], c[i])
                error("An edge of the solution is not in the KEP graph")
            end
            obj_val += instance.edge_weight[c[i-1], c[i]]
            vertex_coverage[c[i]] += 1
        end
    end
    if verbose println("The computed cost of the solution is $obj_val") end

    # check that every vertex is covered at most once by the solution
    for i in 1:nv(g)
        if vertex_coverage[i] >= 2
            printstyled("\n vertex coverage: $vertex_coverage\n" ; color = :red)
            for c in chains
                if i ∈ c
                    printstyled("- by chain: $c\n" ; color = :red)
                end
            end
            for c in cycles
                if i ∈ c
                    printstyled("- by cycle: $c\n" ; color = :red)
                end
            end
            error("A vertex is covered more than once in the solution")
        end
    end
end
