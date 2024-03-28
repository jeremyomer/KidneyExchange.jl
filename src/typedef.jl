"""
$(TYPEDEF)
"""
mutable struct Graph_info
    nb_vertices::Int # number of vertices in the original graph
    nb_pairs::Int  # number of pair vertices
    nb_altruists::Int  # number of altruist vertices
    nb_arcs::Int  # number of arcs in the graph
    density::Float64  # density of the KEP graph

    function Graph_info(instance::Instance)
        g = instance.graph
        n_v = nv(g)
        n_e = ne(g)
        return new(
            n_v,
            instance.nb_pairs,
            instance.nb_altruists,
            n_e,
            n_e / (n_v * (n_v - 1)),
        )
    end
end

"""
$(TYPEDEF)
"""
mutable struct Subgraph_info
    nb_copies::Int  # number of copies
    nb_vertices_per_copy::Float64  # average number of vertices per copy
    nb_arcs_per_copy::Float64  # average number of arcs per copy
    density_per_copy::Float64  # average density per copy

    function Subgraph_info(subgraphs::Graph_copies)
        nb_copies = subgraphs.nb_copies
        if nb_copies > 0
            nb_vertices = [length(findall(subgraphs.is_vertex_list[c])) for c = 1:nb_copies]
            nb_vertices_per_copy = sum(nb_vertices) / nb_copies
            nb_arcs = [
                length(findall(subgraphs.is_arc_list[c])) for
                c = 1:length(subgraphs.is_arc_list)
            ]
            nb_arcs_per_copy = sum(nb_arcs) / nb_copies
            if length(nb_arcs) == nb_copies
                density_per_copy =
                    sum([
                        nb_arcs[c] / (nb_vertices[c] * (nb_vertices[c] - 1)) for
                        c = 1:nb_copies
                    ]) / nb_copies
            else
                density_per_copy = 0.0
            end
            return new(nb_copies, nb_vertices_per_copy, nb_arcs_per_copy, density_per_copy)
        else
            return new(0, 0.0, 0.0, 0.0)
        end
    end
end

"""
$(TYPEDEF)
"""
mutable struct Solution_status
    status::String  # status of OPTIMAL or TIME_LIMIT
    objective_value::Float64  # objective value of the best solution found
    relative_gap::Float64  # final relative optimality gap
    best_cycles::Vector{Vector{Int}}  # selected cycles in the best solution
    best_chains::Vector{Vector{Int}}  # selected chains in the best solution
    node_count::Int  # total number of branch-and-bound nodes explored during the solution proces
    solve_time::Float64  # time spent in the solution process: this might be different from the specified time limit even if status is TIME_LIMIT, since parsing and preprocessing is counted in cpu time

    function Solution_status()
        return new(
            "ON_GOING",
            -Inf,
            Inf,
            Vector{Vector{Int}}(),
            Vector{Vector{Int}}(),
            0,
            0.0,
        )
    end
end

"""
  Mip_model

  Enumerated structure specifying the considered formulation when solving a compact formulation of the kidney exchange problem

  # Enumerate values
  * `HPIEF`: Hybrid position-indexed extended formulation
  * `EXTENDED_EDGE`: Cycles are handled with an extended edge formulation and  chains are handled with position-indexed variables
  * `CYCLE_CUT`: Cycles are handled with an extended edge formulation and  chains are handled with cycle cuts to avoid long cycles
  * `RELAXED_ARC`: Arc formulation where the size of cycles and chains are not considered
"""
@enum Mip_model HPIEF EXTENDED_EDGE CYCLE_CUT RELAXED_ARC

"""
$(TYPEDEF)

  Mutable structure where the solving options of the compact formulation are stored

# Fields

  * `optimizer::String`: LP and IP solver that will be used to solve the master (default is HiGHS for IPs and LPs)
  * `verbose::Bool`: true if messages are printed during the solution (default = true)
  * `model_type::Mip_model`: type of MIP compact model that is to be solved (default = HPIEF)
  * `fvs::Bool`: true if a feedback vertex set is used to reduce the number of graph copies (default = true)
  * `reduce_vertices::Bool`: true if we try deleting useless arcs in graph copies (default = true)
  * `reduce_arcs::Bool`: true if we try deleting useless arcs in graph copies (default = true)
  * `symmetry_break::Bool`: true if the MIP model is modified to reduce the number of optimal solutions (default = true)
  * `nb_threads::Int`: if the LP and IP solver can be called on multiple threads, specify the maximum number of threads that will be used.
"""
mutable struct MIP_params
    optimizer::String
    verbose::Bool
    model_type::Mip_model
    fvs::Bool
    reduce_vertices::Bool
    reduce_arcs::Bool
    symmetry_break::Bool
    nb_threads::Int

    function MIP_params(
        _optimizer::String = "HiGHS",
        _verbose::Bool = true,
        _model_type = HPIEF,
        _fvs = true,
        _reduce_vertices = true,
        _reduce_arcs = true,
        _symmetry = true,
        _nb_threads = 1,
    )
        return new(
            _optimizer,
            _verbose,
            _model_type,
            _fvs,
            _reduce_vertices,
            _reduce_arcs,
            _symmetry,
            _nb_threads,
        )
    end
    function MIP_params(_model_type::Mip_model, _verbose = false)
        return new("HiGHS", _verbose, _model_type, true, true, true, true, 1)
    end
end
