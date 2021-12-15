mutable struct TreeNode
  index::Int  # index of the node
  ub::Float64  # node upper bound
  setzero::Vector{Pair{Int,Int}}  # list of branching arcs set to zero in the column cover, stored as Pair (vertex n => vertex m)
  setone::Vector{Pair{Int,Int}}  # list of of branching arcs set to one in the column cover, stored as Pair (vertex n => vertex m)
  setzero_pief::Vector{Pair{Int,Int}}  # list of branching arcs set to zero in the pief chain cover, stored as Pair (vertex n => vertex m)
  setone_pief::Vector{Pair{Int,Int}}  # list of branching arcs set to one in the pief chain cover, stored as Pair (vertex n => vertex m)

  function TreeNode(node::TreeNode)
    return new(node.index, node.ub, copy(node.setzero), copy(node.setone), copy(node.setzero_pief), copy(node.setone_pief))
  end

  function TreeNode(_index::Int, _ub::Float64, _setzero::Vector{Pair{Int,Int}}, _setone::Vector{Pair{Int,Int}}, _setzero_pief::Vector{Pair{Int,Int}}, _setone_pief::Vector{Pair{Int,Int}})
    return new(_index, _ub, copy(_setzero), copy(_setone), copy(_setzero_pief), copy(_setone_pief))
  end
end

mutable struct Column
  weight:: Float64  # length of the cycle
  vertices::Vector{Int}  # vertices of the cycle in order
  arcs:: Vector{Pair{Int,Int}}  # arcs of the cycle in order
  is_cycle::Bool  # true if the column is a cycle, false if it is a chain

  # use this constructor to add an artificial column
  function Column(path::Vector{Int})
    arcs = Vector{Pair{Int,Int}}(undef, 0)
    for k = 1:(length(path) - 1)
        push!(arcs, (path[k] => path[k+1]))
    end
    return new(0, path, arcs, false)
  end


  function Column(path::Vector{Int}, edge_weight::Matrix{Float64}, _is_cycle::Bool = true)
    # get the arcs of the path
    arcs = Vector{Pair{Int,Int}}(undef, 0)
    for k = 1:(length(path) - 1)
        push!(arcs, (path[k] => path[k+1]))
    end
    if _is_cycle
        push!(arcs, (path[end] => path[1]))
    end
    # compute the weight of the column
    _weight = sum(edge_weight[a[1], a[2]] for a in arcs)

    return new(_weight, path, arcs, _is_cycle)
  end

  function Column(path::Vector{Int}, vertex_weight::Array{Float64}, _is_cycle::Bool = true)
    # get the arcs of the path
    arcs = Vector{Pair{Int,Int}}(undef, 0)
    for k = 1:(length(path) - 1)
        push!(arcs, (path[k] => path[k+1]))
    end
    if _is_cycle
        push!(arcs, (path[end] => path[1]))
    end
    # compute the weight of the column
    _weight = sum(vertex_weight[v] for v in path)

    return new(_weight, path, arcs, _is_cycle)
  end
end

mutable struct BP_params
  optimizer::String  # LP and IP solver that will be used to solve the master
  verbose::Bool  # true if messages are printed during the solution
  is_pief::Bool  # true if the chains are considered in the master model using a position-indexed extended edge formulation
  reduce_vertices::Bool  # true if we try deleting useless vertices in graph copies
  reduce_arcs::Bool  # true if we try deleting useless arcs in graph copies
  is_column_disjoint::Bool  # true if we require column disjoint columns at each column generation iteration
  max_intersecting_columns::Int  # true maximum number of generated columns covering each vertex
  is_tabu_list::Bool  # true if we stop solving a subproblem as soon as it does not produce any positive cost column
  solve_master_IP::Bool  # true if we solve the master IP to find feasible solutions
  time_limit_master_IP::Float64  # time limit at each solution of the master IP
  freq_solve_master_IP::Int  # number of new columns that must be added in the master IP between two solutions of this IP

  function BP_params(_optimizer::String = "GLPK-Cbc", _verbose::Bool = true, _is_pief = false, _reduce_vertices = true, _reduce_arcs = false, _is_column_disjoint = true, _max_intersecting_columns = 6, _is_tabu_list = true, _solve_master_IP = true, _time_limit_IP = 10.0,  _freq_solve_master_IP = 1)
    return new(_optimizer, _verbose, _is_pief, _reduce_vertices, _reduce_arcs, _is_column_disjoint, _max_intersecting_columns, _is_tabu_list, _solve_master_IP, _time_limit_IP, _freq_solve_master_IP)
  end
  function BP_params(_is_pief::Bool, _verbose::Bool = false)
    return new("GLPK-Cbc", _verbose, _is_pief, true, false, true, 6, true, true, 10.0, 1)
  end
end

mutable struct BP_info
  LB::Float64  # best primal value
  UB::Float64  # best dual value
  nb_col_root::Int  # number of columns generated at root node
end

mutable struct BP_status
  bp_info:: BP_info
  status::String # status of ON_GOING, OPTIMAL or TIME_LIMIT
  objective_value::Float64  # objective value of the best solution found
  relative_gap::Float64  # final relative optimality gap
  best_cycles::Vector{Vector{Int}}  # selected cycles in the best primal value
  best_chains::Vector{Vector{Int}}  # selected chains in the best primal value
  node_count::Int  # total number of branch-and-bound nodes explored during the solution proces
  solve_time::Float64  # time spent in the solution process: this might be different from the specified time limit even if status is TIME_LIMIT, since parsing and preprocessing is counted in cpu time
end
