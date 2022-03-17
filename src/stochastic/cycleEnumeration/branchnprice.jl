using Cbc
using DelimitedFiles
using GLPK
using Graphs
using JuMP
using Random
using Requires
using TimerOutputs


function ce_branch_and_price(
        instance::Instance, 
        cycles::Vector{Vector{Cycle}},
        nb_cycles::Int, 
        bp_params::ce_BP_params, 
        timer::TimerOutput, 
        time_limit::Float64
)
    graph = instance.graph

    K = instance.max_cycle_length
    start_time = time()
    nb_vertices = nv(graph)
    nb_types = size(cycles)[1]
    column_pool = Vector{Cycle}()
    if bp_params.mp2
        column_pool = cycles[1]
    end

    tree = Vector{ce_Treenode}()

    push!(tree, TreeNode(1, Inf,Vector{Pair{Int,Int}}(),Vector{Pair{Int,Int}}(),Vector{Pair{Int,Int}}(),Vector{Pair{Int,Int}}(), Vector{Int}(), Vector{Int}(), nb_vertices, 0))

    mastermodel = ce_node_master(instance, column_pool, nb_cycles, nb_types, bp_params, time_limit)

    bp_info = BP_info(-Inf,Inf,0)  # LB=-Inf, UB=Inf, nb_col_root=0
    bp_status = BP_status(bp_info, "ON_GOING", -Inf, Inf, Vector{Vector{Int}}(), Vector{Vector{Int}}(), 1, 0.0, 0)

    while length(tree) >= 1
        current_node = pop!(tree)
        if current_node.ub <= bp_status.bp_info.LB + bp_params.ϵ
            continue
        end

        activate_branching_constraints(mastermodel,  current_node, bp_params)

    end
end

function activate_branching_constraints()

end

function deactivate_branching_constraints()

end

