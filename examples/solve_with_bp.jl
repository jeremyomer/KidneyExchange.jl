# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.10.2
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
using HiGHS

bp_params = BP_params()
bp_params.optimizer = "HiGHS"
timer = TimerOutput()
bp_params.verbose = false
max_time = 7200.0

dataset = @sprintf "%08d" 2  #  (1 - 310) instance index of synthetic kidney donor pools

filename = "00036-" * dataset
dat_file = filename * ".dat"
wmd_file = filename * ".wmd"

Downloads.download("https://www.preflib.org/static/data/kidney/" * dat_file, dat_file)
Downloads.download("https://www.preflib.org/static/data/kidney/" * wmd_file, wmd_file)

new_graph, new_weights, new_is_altruist = read_wmd_file(wmd_file)
# -

cycle_limit, chain_limit = 3, 2
@time bp_status, graph_info, subgraph_info =
    solve_with_BP(filename, cycle_limit, chain_limit, bp_params, timer, max_time)

filename = "MD-00001-00000002"
wmd_file = "MD-00001-00000002" * ".wmd"
dat_file = "MD-00001-00000002" * ".dat"
old_graph, old_weights, old_is_altruist = read_kep_file(wmd_file, dat_file)


@time bp_status, graph_info, subgraph_info =
    solve_with_BP(filename, cycle_limit, chain_limit, bp_params, timer, max_time)

new_weights ≈ old_weights
