"""
$(SIGNATURES)

A compact MIP formulation originally proposed by [Dickerson2016](@cite). The hybrid position indexed edge formulation (HPIEF) handles cycles and chains in the same formulation.
"""
function build_hpief_mip(instance::Instance, subgraphs::Graph_copies, params::MIP_params, maxtime::Float64 = 600)
    if params.verbose println("- build the JuMP model") end

    # extract relevant information from the instance and preprocess
    graph = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length
    nb_subgraphs = subgraphs.nb_copies
    nb_altruists = instance.nb_altruists

    # populate the position index sets for cycle variables
    E = collect(edges(graph))
    κ = Vector{Dict{Pair{Int,Int}, Any}}(undef, nb_subgraphs)
    for l in nb_altruists+1:nb_subgraphs
        κ[l]= Dict{Tuple{Int,Int}, Any}()
        for i in 1:ne(graph)
            u = E[i].src
            v = E[i].dst
            if !subgraphs.is_arc_list[l][i]
                push!(κ[l], (u=>v) => 0:-1)
            elseif u == subgraphs.sources[l]
                push!(κ[l], (u=>v) => 1:1)
            else
                first = subgraphs.d_from_vstar_list[l][u] + 1
                last = K - subgraphs.d_to_vstar_list[l][v]
                push!(κ[l], (u=>v) => first:last)
            end
        end
    end

    # initialize the model
    model = create_model(maxtime, params.optimizer, true, params.verbose)

    # create flow variables and conservation constraints for altruist subgraphs
    @variable(model, y[u in vertices(graph), v in outneighbors(graph, u), k in 1:L], Bin)
    @constraint(model, [u in instance.pairs, k in 1:L-1], sum(y[v,u,k] for v in inneighbors(graph, u)) - sum(y[u,v,k+1] for v in outneighbors(graph, u)) >= 0)
    if L >= 1
        @constraint(model,[u in instance.pairs, v in outneighbors(graph,u)], y[u,v,1] == 0)
        @constraint(model, [u in instance.altruists], sum(y[u,v,1] for v in outneighbors(graph, u)) <= 1)
    end
    @constraint(model, [u in instance.altruists, v in outneighbors(graph,u), k in 2:L], y[u,v,k] == 0)

    # create flow variables and constraints for pair subgraphs
    @variable(model, x[l in nb_altruists+1:nb_subgraphs, u in vertices(graph), v in outneighbors(graph, u), k in κ[l][u=>v]], Bin)

    @constraint(model,[u in instance.pairs, k in 1:K-1, l in nb_altruists+1:nb_subgraphs ; u != subgraphs.sources[l] && subgraphs.is_vertex_list[l][u] == true], sum(x[l,v,u,k] for v in inneighbors(graph, u) if k in κ[l][v=>u]) == sum(x[l,u,v,k+1] for v in outneighbors(graph,u) if k+1 in κ[l][u=>v]))

    @constraint(model,[l in nb_altruists+1:nb_subgraphs], sum(x[l,u,subgraphs.sources[l],k] for u in inneighbors(graph, subgraphs.sources[l]), k in κ[l][u=>subgraphs.sources[l]]) == sum(x[l,subgraphs.sources[l],v,1] for v in outneighbors(graph,subgraphs.sources[l]) if 1 in κ[l][subgraphs.sources[l]=>v]))

    # create vertex disjoint constraints for pairs
    @constraint(model, [u in instance.pairs], sum(x[l,u,v,k] for l in nb_altruists+1:nb_subgraphs, v in outneighbors(graph,u), k in κ[l][u=>v]) + sum(y[v,u,k] for v in inneighbors(graph,u), k in 1:L) <= 1)

    # objective function maximizes the weight of selected cycles
    weight = instance.edge_weight
    if params.symmetry_break
        max_weight = maximum(weight)
        @objective(model, Max, sum((1.0 + l/(max_weight*(instance.nb_pairs)^2))  * weight[u,v] * x[l,u,v,k] for l in nb_altruists+1:nb_subgraphs, u in instance.pairs, v in outneighbors(graph,u), k  in κ[l][u=>v]) + sum(weight[u,v] * y[u,v,k] for  v in instance.pairs, u in inneighbors(graph, v), k  in 1:L))
    else
        @objective(model, Max, sum(weight[u,v] * x[l,u,v,k] for l in nb_altruists+1:nb_subgraphs, u in instance.pairs, v in outneighbors(graph,u), k  in κ[l][u=>v]) + sum(weight[u,v] * y[u,v,k] for  v in instance.pairs, u in inneighbors(graph, v), k  in 1:L))
    end

    return model, κ
end

function solve_hpief_mip(model::Model, params::MIP_params, instance::Instance, κ::Vector{Dict{Pair{Int,Int}, Any}})
    if params.verbose println("- solve the model with $(params.optimizer)") end

    # Local variables
    g = instance.graph
    K = instance.max_cycle_length
    L = instance.max_chain_length

    # Solve the model
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
    mip_status.objective_value = floor(objective_value(model))
    if params.optimizer != "HiGHS"
        mip_status.relative_gap = round(10^4*relative_gap(model))/10^4
    end
    # mip_status.node_count = node_count(model)  # the function has errors in last JuMP version
    mip_status.solve_time = solve_time(model)

    # c. get the solution in a readable format
    x = value.(model[:x])
    # y = value.(model[:y])

    # get the flow corresponding to chains
    if L >= 1
        flow_chain = Dict{Pair{Int,Int}, Float64}()
        for u in vertices(g)
            for v in outneighbors(g, u)
                flow_chain[u=>v] = sum(value(model[:y][u,v,k]) for k in 1:L)
            end
        end

        # get the chains
        for u in instance.altruists
            next = u
            chain = []
            while next != 0
                cur = next
                next = 0
                for v in outneighbors(g, cur)
                    if flow_chain[cur=>v] > 1-ϵ
                        next = v
                        push!(chain, v)
                        break
                    end
                end
            end
            if !isempty(chain)
                pushfirst!(chain, u)
                push!(mip_status.best_chains, chain)
            end
        end
    end

    # get the flow corresponding to cycles
    flow_cycle = Dict{Pair{Int,Int}, Float64}()
    for u in instance.pairs
        for v in outneighbors(g,u)
            flow_cycle[u=>v] = 0
            for l in instance.nb_altruists+1:length(κ)
                if !isempty(κ[l][u=>v])
                    flow_cycle[u=>v] += sum(x[l,u,v,k] for k  in κ[l][u=>v])
                end
            end
        end
    end

    # get the cycles
    is_covered = falses(nv(g))
    for u in instance.pairs
        if is_covered[u] continue end
        next = u
        cycle = []
        while next != 0
            cur = next
            next = 0
            for v in outneighbors(g, cur)
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

    # Output a completion message if in verbose
    if params.verbose
        printstyled("\n----------------------------------------------------------\n The solution of the HPIEF model is complete\n" ; color = :yellow)
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
