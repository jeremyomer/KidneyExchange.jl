# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.10.1
#     language: julia
#     name: julia-1.10
# ---

# # Kidney Data (00036)
#
#
# https://www.preflib.org/dataset/00036

# +
using KidneyExchange
using Printf
using Downloads
using Gurobi

bp_params = BP_params()
bp_params.optimizer = "Gurobi"
timer = TimerOutput()
max_time = 7200.0

dataset = @sprintf "%08d" 2  #  (1 - 310) instance index of synthetic kidney donor pools

filename = "00036-" * dataset
dat_file = filename * ".dat"
wmd_file = filename * ".wmd"

Downloads.download("https://www.preflib.org/static/data/kidney/" * dat_file, dat_file)
Downloads.download("https://www.preflib.org/static/data/kidney/" * wmd_file, wmd_file)

cycle_limit, chain_limit =  3, 2

@time bp_status, graph_info, subgraph_info = solve_with_BP(filename, cycle_limit, chain_limit, bp_params, timer, max_time);

# +
bp_params.is_pief = true

@time bp_status, graph_info, subgraph_info = solve_with_BP(filename, cycle_limit, chain_limit, bp_params, timer, max_time);
# -


