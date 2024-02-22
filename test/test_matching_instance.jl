import Downloads
using Printf

@testset "Read and write preflib files" begin

    dataset = @sprintf "%08d" 1
    
    prefix = "00036-" * dataset 
    dat_file = prefix * ".dat"
    wmd_file = prefix * ".wmd"
    
    Downloads.download("https://www.preflib.org/static/data/kidney/" * dat_file, dat_file)
    Downloads.download("https://www.preflib.org/static/data/kidney/" * wmd_file, wmd_file)
    
    KidneyExchange.read_wmd_file( wmd_file )
    
    bp_status, graph_info = solve_with_BP( prefix, 3, 3 )
    
    @test bp_status.best_cycles ≈ [[3, 8], [1, 6]]
    @test bp_status.objective_value ≈ 4.0

end
