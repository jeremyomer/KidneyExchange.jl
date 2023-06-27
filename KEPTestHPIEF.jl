include("src/KidneyExchange.jl")
using Main.KidneyExchange
using Gurobi
using ArgParse
using TimerOutputs

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

filename = parsed_args["filename"] #
cycle_limit = parsed_args["k"]
chain_limit = parsed_args["l"]


mip_params = MIP_params()
mip_params.optimizer = "Gurobi"
mip_params.model_type = HPIEF
timer = TimerOutput()
max_time = 600.0

solve_with_mip("heterogeneous/heterogeneous_128_0_1", 3, 4, mip_params, timer, max_time)

status, graph_info, subgraph_info = solve_with_mip(filename, cycle_limit, chain_limit, mip_params, timer, max_time);

t_build_MIP = TimerOutputs.time(timer["Build MIP"])/10^9
t_solve_MIP = TimerOutputs.time(timer["Solve MIP"])/10^9
t_total = TimerOutputs.time(timer["Parser"]) + TimerOutputs.time(timer["Preprocessing"]) + TimerOutputs.time(timer["Build MIP"]) + TimerOutputs.time(timer["Solve MIP"])

stockfilename = "$(filename)_$(cycle_limit)_$(chain_limit)"
mip_log_file_name = string("/BP_profiling/mip_info_", stockfilename, ".csv")
mip_log = open(string(dirname(@__FILE__),mip_log_file_name), "a")
println(mip_log, "$filename; $(graph_info.nb_pairs);  $(graph_info.nb_altruists); $(cycle_limit); $(chain_limit); $(mip_params.model_type); $(status.objective_value); $(round(t_total/10^8)/10); $(status.relative_gap); $(round(10*t_solve_MIP)/10); $(status.node_count); $(round(10*t_build_MIP)/10);")
close(mip_log)
