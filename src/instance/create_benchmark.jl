include("pool_generator.jl")
include("saidman_pool_generator.jl")
include("heterogeneous_pool_generator.jl")

"""
    generate_saidman_instance

Generate a KEP dataset and write in two text files with extensions .dat and .wmd. This generator is usually refererred to as the "Saidman Generator".

# Arguments
* `nb_pairs::Int`: number of pairs of incompatible donor and patient
* `nb_altruists::Int`: number of altruist donors
* `index::Int`: index of the dataset, which will appear in the filename

# Reference

[Saidman2014](@cite)

"""
function generate_saidman_instance(nb_pairs::Int, nb_altruists::Int, index::Int)
    kep_graph, edge_weights, donorBT, patientBT, wifeP, patientPRA, is_altruist = generate_saidman_kep_graph(nb_pairs, nb_altruists);
    write_kep_file(kep_graph, edge_weights, donorBT, patientBT, wifeP, patientPRA, is_altruist, "saidman/saidman_$(nb_pairs)_$(nb_altruists)_$(index)");
    kep_graph = 0;
    return nothing;
end

"""
    generate_sparse_unos_instance

Generate a KEP dataset and write in two text files with extensions .dat and .wmd. This generator is obtained by applying the Saidman generator with small modifications in the probabilities of compatibility to roughtly mimic the UNOS pool as of April 15, 2013.  The corresponding data was taken from the KPD Work
Group Data Analysis - CMR - June 2013 report.
Seel also [Dickerson2012](@cite) for the motivation and description of this generator.

# Arguments
* `nb_pairs::Int`: number of pairs of incompatible donor and patient
* `nb_altruists::Int`: number of altruist donors
* `index::Int`: index of the dataset, which will appear in the filename
 """
function generate_sparse_unos_instance(nb_pairs::Int, nb_altruists::Int, index::Int)
    kep_graph, edge_weights, donorBT, patientBT, wifeP, patientPRA, is_altruist = generate_sparse_unos_kep_graph(nb_pairs, nb_altruists);
    write_kep_file(kep_graph, edge_weights, donorBT, patientBT, wifeP, patientPRA, is_altruist, "sparse/sparse_$(nb_pairs)_$(nb_altruists)_$(index)");
    kep_graph = 0;
    return nothing;
end

"""
    generate_heterogeneous_instance(nb_pairs, nb_altruists, index)

Generate a KEP dataset and write in two text files with extensions .dat and .wmd; the graph generator is based on the following article: Kidney Exchange in Dynamic Sparse Heterogeneous Pools: [Ashlagi2013](@cite).

# Arguments
* `nb_pairs::Int`: number of pairs of incompatible donor and patient
* `nb_altruists::Int`: number of altruist donors
* `index::Int`: index of the dataset, which will appear in the filename
"""
function generate_heterogeneous_instance(nb_pairs::Int, nb_altruists::Int, index::Int)
    kep_graph, edge_weights, donorBT, patientBT, wifeP, patientPRA, is_altruist = generate_heterogeneous_kep_graph(nb_pairs, nb_altruists);
    write_kep_file(kep_graph, edge_weights, donorBT, patientBT, wifeP, patientPRA, is_altruist, "heterogeneous/heterogeneous_$(nb_pairs)_$(nb_altruists)_$(index)");
    kep_graph = 0;
    return nothing;
end

"""
	generate_abraham_benchmark()

Generate a benchmark using the same method as that used in [Abraham2007](@cite)

"""
function generate_abraham_benchmark()
    Random.seed!(30112021)
    list_nb_pairs = [100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    for nb_pairs in list_nb_pairs
        for k in 1:10
            generate_saidman_instance(nb_pairs, 0, k)
        end
    end
end

"""
	generate_complete_benchmark()

Generate all the instances that are used in addition to the PrefLib to assess the solution methods of this package in the article describing the package.

*add citation later*

The instances are stored in the subdirectories `heterogeneous`, `sparse` and `saidman` of the `data` directory. Beware that if using the package with Pkg.add(), the `data` directory is in the directory where the package is stored.
"""
function generate_complete_benchmark()
    Random.seed!(30112021)
    list_nb_pairs = [128,256,512,1024,2048,6000,1000]
    for nb_pairs in list_nb_pairs
        for k in 1:10
            generate_heterogeneous_instance(nb_pairs, round(Int, 0.10*nb_pairs), k)
            generate_sparse_unos_instance(nb_pairs, round(Int, 0.10*nb_pairs), k)
        end
    end
    list_nb_pairs = [6000,1000]
    for nb_pairs in list_nb_pairs
        for k in 1:10
            generate_saidman_instance(nb_pairs, round(Int, 0.10*nb_pairs), k)
        end
    end
end
