function ce_node_master(
        instance::Instance,
        column_pool::Vector{Cycle},
        nb_cycles::Int,
        nb_types::Int,
        bp_params::BP_params = ce_BP_params(),
        time_limit::Float64 = 10000.0
)
    graph = instance.graph

    master = create_model(time_limit, bp_params.optimizer, false, false)

    @variable(master, y[c in 1:length(column_pool)] >= 0)

    

    master[:branch_one] = Dict{Pair{Int64, Int64}, ConstraintRef}()
    master[:branch_zero] = Dict{Pair{Int64, Int64}, ConstraintRef}()
    master[:branch_one_vertex] = Dict{Int, ConstraintRef}()
    master[:branch_zero_vertex] = Dict{Int, ConstraintRef}()

end

function master_add_cycle()

end