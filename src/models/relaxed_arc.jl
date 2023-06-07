"""
$(SIGNATURES)

Deterministic relaxed-arc formulation for a KEP problem. No cycle length
limitation is imposed.

# Parameters
* `instance::Instance`: KEP instance to solve
* `maxtime::Real=60` : Maximum solving time in seconds
"""
function relaxed_arc(instance::Instance, params::MIP_params, maxtime::Float64 = 600)
    if params.verbose println("- build the JuMP model") end
    graph = instance.graph
    E = [(e.src=>e.dst) for e in edges(graph)]

    model = create_model(maxtime, params.optimizer, true, params.verbose)

    # Variables :
    #   x_e = 1 if the transplant e is carried out
    #   x_e = 0 else
    @variable(model, x[e in E], Bin, base_name="x")

    # Objective :
    #   Maximize the total transplant utility* `instance::Instance`: path of the input data files, this should include the name of the files, but not the .dat and .wmd extensions
    @objective(model, Max, sum(instance.edge_weight[e[1], e[2]] * x[e] for e in E))

    # Flow constraint:
    #   Each pair vertex has as many entering as exiting arcs
    @constraint(model, flow_pairs[v in instance.pairs], sum(x[u=>v] for u in inneighbors(graph, v)) == sum(x[v=>w] for w in outneighbors(graph, v)))

    # Vertex cover:
    #  Each vertex has at most one leaving edge
    @constraint(model, vertex_cover[v in vertices(graph)], sum(x[v=>w] for w in outneighbors(graph, v)) <= 1)

    # Solve the model
    if params.verbose println("- solve the model with $(params.optimizer)") end
    optimize!(model)

    # Build the MIP status structure
    mip_status = Solution_status()

    # a. first get the solution status
    if termination_status(model) == MOI.OPTIMAL
        mip_status.status = "OPTIMAL"
    elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
        mip_status.status = "TIME_LIMIT"
    else
        mip_status.status = "NOT_FEASIBLE"
        printstyled("- no solution was found\n" ; color = :red)
        return mip_status
    end

    # b. get the main characteristics of the solution
    mip_status.objective_value = floor(objective_value(model) + ϵ)
    if params.optimizer != "HiGHS"
        mip_status.relative_gap = round(10^4*relative_gap(model))/10^4
    end
    # mip_status.node_count = node_count(model)  # the function has errors in last JuMP version
    mip_status.solve_time = solve_time(model)


    # c. get the solution in a readable format, i.e. get the cycles corresponding to the solution
    flow_cycle = value.(x)
    is_covered = falses(nv(graph))
    for u in instance.pairs
        if is_covered[u] continue end
        next = u
        cycle = []
        while next != 0
            cur = next
            next = 0
            for v in outneighbors(graph, cur)
                if is_covered[v] continue end
                if flow_cycle[cur=>v] > 1-ϵ
                    next = v
                    push!(cycle, v)
                    is_covered[v] = true
                    break
                end
            end
        end
        if !isempty(cycle)
            is_covered[u] = true
            if cycle[end] != u
                error("not a cycle")
            end
            push!(mip_status.best_cycles, cycle)
        end
    end

    if params.verbose
        printstyled("\n----------------------------------------------------------\n The solution of the RELAXED ARC model is complete\n" ; color = :yellow)
        if termination_status(model) == MOI.OPTIMAL
            printstyled("- the solution is optimal\n" ; color = :yellow)
            printstyled("- best solution found: value $(mip_status.objective_value)\n" ; color = :yellow)
            printstyled("----------------------------------------------------------\n\n" ; color = :yellow)
        elseif termination_status(model) == MOI.TIME_LIMIT
            printstyled("- the time limit is exceeded\n" ; color = :yellow)
            if params.optimizer != "HiGHS"
                printstyled("- best solution found: value $(mip_status.objective_value) with gap $(mip_status.relative_gap) %\n" ; color = :yellow)
            else
                printstyled("- best solution found: value $(mip_status.objective_value)\n" ; color = :yellow)
            end
            printstyled("----------------------------------------------------------\n\n" ; color = :yellow)
        end
    end

    return mip_status
end
