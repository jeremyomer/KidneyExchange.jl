using KidneyExchange
using Test

include("test_matching_instance.jl")

const KEP = KidneyExchange

@testset "KidneyExchange.jl" begin
    tests_K3_L0 = Dict(
        joinpath("tests_cycles", "MD-00001-00000001") => 4,
        joinpath("tests_cycles", "MD-00001-00000120") => 83,
        joinpath("tests_cycles", "heterogeneous_128_0_1") => 85,
        joinpath("tests_cycles", "sparse_128_0_1") => 28,
    )
    # "tests_cycles/MD-00001-00000160" => 159, "tests_cycles/heterogeneous_256_0_1" => 170, "tests_cycles/sparse_256_0_1" => 100

    for filename in keys(tests_K3_L0)
        printstyled("\nRun every algorithm on $filename with K=3, L=0\n"; bold = true)        # test branch-and-price algorithms
        stats = KidneyExchange.solve_with_BP(filename, 3, 0, BP_params(false))
        println("\t", @test stats[1].objective_value == tests_K3_L0[filename])

        stats = KidneyExchange.solve_with_BP(filename, 3, 0, BP_params(true))
        println("\t", @test stats[1].objective_value == tests_K3_L0[filename])

        # test compact formulations
        stats =
            KidneyExchange.solve_with_mip(filename, 3, 0, MIP_params(KidneyExchange.HPIEF))
        println("\t", @test stats[1].objective_value == tests_K3_L0[filename])

        stats = KidneyExchange.solve_with_mip(
            filename,
            3,
            0,
            MIP_params(KidneyExchange.EXTENDED_EDGE),
        )
        println("\t", @test stats[1].objective_value == tests_K3_L0[filename])

        stats = KidneyExchange.solve_with_mip(
            filename,
            3,
            0,
            MIP_params(KidneyExchange.RELAXED_ARC),
        )
        println("\t", @test stats[1].objective_value >= tests_K3_L0[filename])
    end

    tests_K4_L0 = Dict(
        joinpath("tests_cycles", "MD-00001-00000001") => 4,
        joinpath("tests_cycles", "MD-00001-00000120") => 86,
        joinpath("tests_cycles", "heterogeneous_128_0_1") => 90,
        joinpath("tests_cycles", "sparse_128_0_1") => 36,
    )
    # "tests_cycles/MD-00001-00000160"=>159, "tests_cycles/heterogeneous_256_0_1" => 174, "tests_cycles/sparse_256_0_1" => 127
    for filename in keys(tests_K4_L0)
        printstyled("\nRun every algorithm on $filename with K=4, L=0\n"; bold = true)        # test branch-and-price algorithms
        stats = KidneyExchange.solve_with_BP(filename, 4, 0, BP_params(false))
        println("\t", @test stats[1].objective_value == tests_K4_L0[filename])

        stats = KidneyExchange.solve_with_BP(filename, 4, 0, BP_params(true))
        println("\t", @test stats[1].objective_value == tests_K4_L0[filename])

        # test compact formulations
        stats =
            KidneyExchange.solve_with_mip(filename, 4, 0, MIP_params(KidneyExchange.HPIEF))
        println("\t", @test stats[1].objective_value == tests_K4_L0[filename])

        stats = KidneyExchange.solve_with_mip(
            filename,
            4,
            0,
            MIP_params(KidneyExchange.EXTENDED_EDGE),
        )
        println("\t", @test stats[1].objective_value == tests_K4_L0[filename])

        stats = KidneyExchange.solve_with_mip(
            filename,
            4,
            0,
            MIP_params(KidneyExchange.RELAXED_ARC),
        )
        println("\t", @test stats[1].objective_value >= tests_K4_L0[filename])
    end

    tests_K3_L6 = Dict(
        joinpath("tests_chains", "MD-00001-00000015") => 16,
        joinpath("tests_chains", "MD-00001-00000127") => 82,
        joinpath("tests_chains", "heterogeneous_128_19_1") => 102,
        joinpath("tests_chains", "sparse_128_19_1") => 51,
    )

    for filename in keys(tests_K3_L6)
        printstyled("\nRun every algorithm on $filename with K=3, L=6\n"; bold = true)
        # test branch-and-price algorithms
        stats = KidneyExchange.solve_with_BP(filename, 3, 6, BP_params(false))
        println("\t", @test stats[1].objective_value == tests_K3_L6[filename])

        stats = KidneyExchange.solve_with_BP(filename, 3, 6, BP_params(true))
        println("\t", @test stats[1].objective_value == tests_K3_L6[filename])

        # test compact formulations
        stats =
            KidneyExchange.solve_with_mip(filename, 3, 6, MIP_params(KidneyExchange.HPIEF))
        println("\t", @test stats[1].objective_value == tests_K3_L6[filename])

        stats = KidneyExchange.solve_with_mip(
            filename,
            3,
            6,
            MIP_params(KidneyExchange.EXTENDED_EDGE),
        )
        println("\t", @test stats[1].objective_value == tests_K3_L6[filename])
    end
end
