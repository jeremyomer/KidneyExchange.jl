include("src/KEP.jl")
using Main.KEP
using ArgParse
using Gurobi
using TimerOutputs
using CPLEX

argparse_settings = ArgParseSettings()
@add_arg_table! argparse_settings begin
    "filename"
        help = "file name of the instance to run"
        arg_type = String
        required = true
    "k"
        help = "maximum number of arcs per cycle"
        arg_type = Int
        required = true
    "l"
        help = "maximum number of arcs per chain"
        arg_type = Int
        required = true
end
parsed_args = parse_args(ARGS, argparse_settings)

filename = parsed_args["filename"]
cycle_limit = parsed_args["k"]
chain_limit = parsed_args["l"]

bp_params = BP_params()
bp_params.optimizer = "CPLEX"
timer = TimerOutput()
max_time = 600.0

solve_with_BP("heterogeneous/heterogeneous_128_0_1", 3, 4, bp_params, timer, max_time)

bp_status, graph_info, subgraph_info = solve_with_BP(filename, cycle_limit, chain_limit, bp_params, timer, max_time)

t_total =  TimerOutputs.time(timer["Parser"]) + TimerOutputs.time(timer["Preprocessing"]) + TimerOutputs.time(timer["B&P"])

t_parser = 100*TimerOutputs.time(timer["Parser"])/t_total
t_preprocess = 100*TimerOutputs.time(timer["Preprocessing"])/t_total
t_node_master = 100*TimerOutputs.time(timer["B&P"]["Process_Node"]["opt_master"]) / t_total
n_call_node_master = TimerOutputs.ncalls(timer["B&P"]["Process_Node"]["opt_master"])
t_bellman = 100*TimerOutputs.time(timer["B&P"]["Process_Node"]["Bellman-Ford"])/t_total
ncall_bellman = TimerOutputs.ncalls(timer["B&P"]["Process_Node"]["Bellman-Ford"])
if haskey(timer["B&P"]["Process_Node"], "Bellman-Ford-chain")
    t_bellman += 100*TimerOutputs.time(timer["B&P"]["Process_Node"]["Bellman-Ford-chain"])/t_total
    ncall_bellman += TimerOutputs.ncalls(timer["B&P"]["Process_Node"]["Bellman-Ford-chain"])
end

t_mip_master = 0.0
ncall_mip_master = 0
if haskey(timer["B&P"]["Process_Node"], "IP_master")
    t_mip_master = 100*TimerOutputs.time(timer["B&P"]["Process_Node"]["IP_master"]) / t_total
    ncall_mip_master = TimerOutputs.ncalls(timer["B&P"]["Process_Node"]["IP_master"])
end

# store collected data
stockfilename = "$(filename)_$(cycle_limit)_$(chain_limit)"
graph_log_file_name = string("/BP_profiling/graph_info_", stockfilename, ".csv")
graph_log = open(string(dirname(@__FILE__),graph_log_file_name), "a")

bp_log_file_name = string("/BP_profiling/bp_info_", stockfilename, ".csv")
bp_log = open(string(dirname(@__FILE__),bp_log_file_name), "a")

println(graph_log, "$filename; $(graph_info.nb_pairs); $(graph_info.nb_altruists); $(round(100*graph_info.density)/100); $(subgraph_info.nb_copies); $(round(Int, subgraph_info.nb_vertices_per_copy))")

println(bp_log, "$filename; $(graph_info.nb_pairs); $(graph_info.nb_altruists); $cycle_limit; $chain_limit; BP; $(bp_status.objective_value); $(round(t_total/10^8)/10); $(100*bp_status.relative_gap); $(bp_status.bp_info.nb_col_root); $(bp_status.node_count); $(round(10*t_node_master)/10); $n_call_node_master; $(round(10*t_bellman)/10); $ncall_bellman; $(round(10*t_mip_master)/10) %; $ncall_mip_master; $(round(10*t_parser)/10) %; $(round(10*t_preprocess)/10) %")

close(graph_log)
close(bp_log)
