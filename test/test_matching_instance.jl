import Downloads
using KidneyExchange
using Printf
using Test

@testset "Read and write preflib files" begin

    dataset = @sprintf "%08d" 1

    prefix = "00036-" * dataset
    dat_file = prefix * ".dat"
    wmd_file = prefix * ".wmd"

    Downloads.download("https://www.preflib.org/static/data/kidney/" * dat_file, dat_file)
    Downloads.download("https://www.preflib.org/static/data/kidney/" * wmd_file, wmd_file)

    KidneyExchange.read_wmd_file(wmd_file)

    bp_status, graph_info = solve_with_BP(prefix, 3, 3)

    @test bp_status.best_cycles ≈ [[3, 8], [1, 6]]
    @test bp_status.objective_value ≈ 4.0

    nb_pairs = 500
    nb_altruists = 25

    kep_graph, edge_weights, donorBT, patientBT, wifeP, patientPRA, is_altruist =
        KidneyExchange.generate_sparse_unos_kep_graph(nb_pairs, nb_altruists)

    filepath = @sprintf "sparse%05d%08d" nb_pairs nb_altruists
    title = "Sparse unos graph with $nb_pairs pairs and $nb_altruists altruists"
    description = "Saidman generator with small modifications in the probabilities of compatibility to roughtly mimic the UNOS pool"
    KidneyExchange.write_wmd_file(
        kep_graph,
        edge_weights,
        is_altruist,
        filepath,
        title,
        description,
    )
    KidneyExchange.write_dat_file(
        kep_graph,
        donorBT,
        patientBT,
        wifeP,
        patientPRA,
        is_altruist,
        filepath,
    )

    new_dataset = generate_heterogeneous_instance(nb_pairs, nb_altruists, index = 1)

    bp_status, graph_info = solve_with_BP(new_dataset, 3, 3)

end
