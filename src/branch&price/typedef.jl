mutable struct TreeNode
  index::Int  # index of the node
  ub::Float64  # node upper bound
  setzero::Vector{Pair{Int,Int}}  # list of branching arcs set to zero in the column cover, stored as Pair (vertex n => vertex m)
  setone::Vector{Pair{Int,Int}}  # list of of branching arcs set to one in the column cover, stored as Pair (vertex n => vertex m)
  setzero_pief::Vector{Pair{Int,Int}}  # list of branching arcs set to zero in the pief chain cover, stored as Pair (vertex n => vertex m)
  setone_pief::Vector{Pair{Int,Int}}  # list of branching arcs set to one in the pief chain cover, stored as Pair (vertex n => vertex m)
  setzero_vertex::Vector{Int}  # list of branching vertices set to zero in the column cover
  setone_vertex::Vector{Int}  # list of of branching vertices set to one in the column cover
  nb_cycles_max::Vector{Int}  # maximum number of cycles for each feasible length (k=2 to K)
  nb_cycles_min::Vector{Int}  # minimum number of cycles for each feasible length (k=2 to K)
  nb_chains_max::Vector{Int}  # maximum number of chains for each feasible length (l=1 to L)
  nb_chains_min::Vector{Int}  # minimum number of chains for each feasible length (l=1 to L)

  function TreeNode(node::TreeNode)
    return new(node.index, node.ub, copy(node.setzero), copy(node.setone), copy(node.setzero_pief), copy(node.setone_pief), copy(node.setzero_vertex), copy(node.setone_vertex), copy(node.nb_cycles_max), copy(node.nb_cycles_min), copy(node.nb_chains_max), copy(node.nb_chains_min))
  end

  function TreeNode(_index::Int, _ub::Float64, _setzero::Vector{Pair{Int,Int}}, _setone::Vector{Pair{Int,Int}}, _setzero_pief::Vector{Pair{Int,Int}}, _setone_pief::Vector{Pair{Int,Int}},  _setzero_vertex::Vector{Int}, _setone_vertex::Vector{Int}, _nb_cycles_max::Vector{Int}, _nb_cycles_min::Vector{Int}, _nb_chains_max::Vector{Int}, _nb_chains_min::Vector{Int})
    return new(_index, _ub, copy(_setzero), copy(_setone), copy(_setzero_pief), copy(_setone_pief), copy(_setzero_vertex), copy(_setone_vertex), copy(_nb_cycles_max), copy(_nb_cycles_min), copy(_nb_chains_max), copy(_nb_chains_min))
  end
end

mutable struct Column
  weight:: Float64  # length of the cycle
  vertices::Vector{Int}  # vertices of the cycle in order
  arcs:: Vector{Pair{Int,Int}}  # arcs of the cycle in order
  is_cycle::Bool  # true if the column is a cycle, false if it is a chain
  length::Int64  # number of vertices covered by the column

  # use this constructor to add an artificial column
  function Column(path::Vector{Int})
    arcs = Vector{Pair{Int,Int}}(undef, 0)
    for k = 1:(length(path) - 1)
        push!(arcs, (path[k] => path[k+1]))
    end
    return new(0, path, arcs, false, length(path))
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
    _length = _is_cycle ? length(path) : length(path) - 1

    return new(_weight, path, arcs, _is_cycle, _length)
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
    _length = _is_cycle ? length(path) : length(path) - 1


    return new(_weight, path, arcs, _is_cycle, _length)
  end
end


"""
  BP_params

  Mutable structure where the options of the branch-and-price solver are stored

  # Fields
  * `optimizer::String`: LP and IP solver that will be used to solve the master (default is Cbc for IPs and GLPK for LPs)
  * `verbose::Bool`: true if messages are printed during the solution (default = true)
  * `is_pief::Bool`: true if the chains are considered in the master model using a position-indexed extended edge formulation (default = false)
  * `fvs::Bool`: true if a feedback vertex set is used to reduce the number of graph copies (default = true)
  * `reduce_vertices::Bool`: true if we try deleting useless vertices in graph copies (default = true)
  * `is_column_disjoint::Bool`: true if we require column disjoint columns at each column generation iteration (default = true)
  * `max_intersecting_columns::Int`: true maximum number of generated columns covering each vertex (default = 6)
  * `is_tabu_list::Bool`: true if we stop solving a subproblem as soon as it does not produce any positive cost column; this subproblem will be considered again when proving optimality (default = true)
  * `solve_master_IP::Bool`: true if we solve the master IP to find feasible solutions (default = true)
  * `time_limit_master_IP::Float64`: time limit (seconds) at each solution of the master IP (default = 10.0)
  * `freq_solve_master_IP::Int`: number of new columns that must be added in the master IP between two solutions of this IP (default = 1)
  * `restart_for_IP::Bool`: true if the root node can be solved twice to generate more columns when the IP master could not prove optimality of the relaxation value (default = true)
"""
mutable struct BP_params
  optimizer::String
  verbose::Bool
  is_pief::Bool
  fvs::Bool
  reduce_vertices::Bool
  is_column_disjoint::Bool
  max_intersecting_columns::Int
  is_tabu_list::Bool
  solve_master_IP::Bool
  time_limit_master_IP::Float64
  freq_solve_master_IP::Int
  restart_for_IP::Bool
  branch_on_vertex::Bool

  function BP_params(_optimizer::String = "GLPK-Cbc", _verbose::Bool = true, _is_pief = false, _fvs = true,  _reduce_vertices = true, _is_column_disjoint = true, _max_intersecting_columns = 6, _is_tabu_list = true, _solve_master_IP = true, _time_limit_IP = 30.0,  _freq_solve_master_IP = 2, _restart_for_IP = true, _branch_on_vertex = false)
    return new(_optimizer, _verbose, _is_pief, _fvs, _reduce_vertices, _is_column_disjoint, _max_intersecting_columns, _is_tabu_list, _solve_master_IP, _time_limit_IP, _freq_solve_master_IP, _restart_for_IP, _branch_on_vertex)
  end
  function BP_params(_is_pief::Bool, _verbose::Bool = false)
    return new("GLPK-Cbc", _verbose, _is_pief, true, true, true, 6, true, true, 30.0, 2, true, false)
  end
end

"""
  BP_info

  Mutable structure where extra information about the branch-and-price execution is stored

  # Fields
  * `LB::Float64`: best primal value
  * `UB::Float64`: best dual value
  * `nb_col_root::Int`: number of columns generated at root node
"""
mutable struct BP_info
  LB::Float64
  UB::Float64
  nb_col_root::Int
end

"""
  BP_status

  Mutable structure where the results of the branch-and-price are stored

  # Fields
  * `bp_info::BP_info`: extra info about the branch-and-price execution
  * `status::String`: status of solution: ON_GOING, OPTIMAL or TIME_LIMIT
  * `objective_value::Float64`: objective value of the best solution found
  * `relative_gap::Float64`: final relative optimality gap
  * `best_cycles::Vector{Vector{Int}}`: selected cycles in the best primal value
  * `best_chains::Vector{Vector{Int}}`: selected chains in the best primal value
  * `node_count::Int`: total number of branch-and-bound nodes explored during the solution proces
  * `solve_time::Float64`: time spent in the solution process; this might be different from the specified time limit even if status is TIME_LIMIT, since parsing and preprocessing is counted in cpu time
  * `nb_cols_last_ip::Int`: number of columns in the master IP at last solution
"""
mutable struct BP_status
  bp_info:: BP_info
  status::String
  objective_value::Float64
  relative_gap::Float64
  best_cycles::Vector{Vector{Int}}
  best_chains::Vector{Vector{Int}}
  node_count::Int
  solve_time::Float64
  nb_cols_last_ip::Int
  node_count_last_ip::Int
  termination_status_last_ip
end
