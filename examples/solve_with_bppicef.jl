# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# +
using KidneyExchange
using ArgParse
using HiGHS
using TimerOutputs
using Downloads

filename = "00036-00000001"
cycle_limit = 3
chain_limit = 2

bp_params = BP_params()
bp_params.optimizer = "HiGHS"
bp_params.is_pief = true
timer = TimerOutput()
max_time = 7200.0

dat_file = filename * ".dat"
wmd_file = filename * ".wmd"

Downloads.download("https://www.preflib.org/static/data/kidney/" * dat_file, dat_file)
Downloads.download("https://www.preflib.org/static/data/kidney/" * wmd_file, wmd_file)


bp_status, graph_info, subgraph_info = solve_with_BP(filename, cycle_limit, chain_limit, bp_params, timer, max_time)



# -


