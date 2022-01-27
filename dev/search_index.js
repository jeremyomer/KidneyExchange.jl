var documenterSearchIndex = {"docs":
[{"location":"#KidneyExchange.jl","page":"KidneyExchange.jl","title":"KidneyExchange.jl","text":"","category":"section"},{"location":"","page":"KidneyExchange.jl","title":"KidneyExchange.jl","text":"Documentation for KidneyExchange.jl","category":"page"},{"location":"","page":"KidneyExchange.jl","title":"KidneyExchange.jl","text":"Modules = [KidneyExchange]\nOrder   = [:type, :function]","category":"page"},{"location":"#KidneyExchange.BP_info","page":"KidneyExchange.jl","title":"KidneyExchange.BP_info","text":"BP_info\n\nMutable structure where extra information about the branch-and-price execution is stored\n\nFields\n\nLB::Float64: best primal value\nUB::Float64: best dual value\nnb_col_root::Int: number of columns generated at root node\n\n\n\n\n\n","category":"type"},{"location":"#KidneyExchange.BP_params","page":"KidneyExchange.jl","title":"KidneyExchange.BP_params","text":"BP_params\n\nMutable structure where the options of the branch-and-price solver are stored\n\nFields\n\noptimizer::String: LP and IP solver that will be used to solve the master (default is Cbc for IPs and GLPK for LPs)\nverbose::Bool: true if messages are printed during the solution (default = true)\nis_pief::Bool: true if the chains are considered in the master model using a position-indexed extended edge formulation (default = false)\nfvs::Bool: true if a feedback vertex set is used to reduce the number of graph copies (default = true)\nreduce_vertices::Bool: true if we try deleting useless vertices in graph copies (default = true)\nis_column_disjoint::Bool: true if we require column disjoint columns at each column generation iteration (default = true)\nmax_intersecting_columns::Int: true maximum number of generated columns covering each vertex (default = 6)\nis_tabu_list::Bool: true if we stop solving a subproblem as soon as it does not produce any positive cost column; this subproblem will be considered again when proving optimality (default = true)\nsolve_master_IP::Bool: true if we solve the master IP to find feasible solutions (default = true)\ntime_limit_master_IP::Float64: time limit (seconds) at each solution of the master IP (default = 10.0)\nfreq_solve_master_IP::Int: number of new columns that must be added in the master IP between two solutions of this IP (default = 1)\nrestart_for_IP::Bool: true if the root node can be solved twice to generate more columns when the IP master could not prove optimality of the relaxation value (default = true)\n\n\n\n\n\n","category":"type"},{"location":"#KidneyExchange.BP_status","page":"KidneyExchange.jl","title":"KidneyExchange.BP_status","text":"BP_params\n\nMutable structure where the results of the branch-and-price are stored\n\nFields\n\nbp_info::BP_info: extra info about the branch-and-price execution\nstatus::String: status of solution: ONGOING, OPTIMAL or TIMELIMIT\nobjective_value::Float64: objective value of the best solution found\nrelative_gap::Float64: final relative optimality gap\nbest_cycles::Vector{Vector{Int}}: selected cycles in the best primal value\nbest_chains::Vector{Vector{Int}}: selected chains in the best primal value\nnode_count::Int: total number of branch-and-bound nodes explored during the solution proces\nsolve_time::Float64: time spent in the solution process; this might be different from the specified time limit even if status is TIME_LIMIT, since parsing and preprocessing is counted in cpu time\n\n\n\n\n\n","category":"type"},{"location":"#KidneyExchange.Graph_copies","page":"KidneyExchange.jl","title":"KidneyExchange.Graph_copies","text":"Graph_copies\n\nMutable structure describing the copies of the graph: one copy per vertex or one copy per vertex of a feedback vertex set if the option is set.   It is important to note that the copies related to altruist donors always appear first in the list of copies\n\nFields\n\nsources::Vector{Int}: source of each graph copy\nis_vertex_list::Vector{BitVector}: for each copy and each vertex, true if the vertex belongs to the copy\nd_to_vstar_list::Vector{Vector{Int}}: in each copy, distance from each vertex to the source\nd_from_vstar_list::Vector{Vector{Int}}: in each copy, distance from each vertex to the source\nnb_copies::Int: number of graph copies\nis_arc_list::Vector{BitVector}: for each copy and each arc index, true if the arc appears in the copy\nchain_mip::Model: shared MIP model that will be solved every time a chain subproblem needs to be solved to optimality\n\n\n\n\n\n","category":"type"},{"location":"#KidneyExchange.Instance","page":"KidneyExchange.jl","title":"KidneyExchange.Instance","text":"Instance\n\nNon-mutable structure describing the instance that is solved\n\nFields\n\ngraph::SimpleDiGraph: graph describing compatibilities between pairs and  with altruist donors\nvertex_weight::Vector{Float64}: weight of each vertex\nedge_weight::Matrix{Float64}: weight of each arc (0 if there is no arc)\npairs::Vector{Int}: indices of the donor/patient pair vertices\naltruists::Vector{Int}: indices of the altruist donor vertices\nnb_pairs::Int: number of donor/patient pairs\nnb_altruists::Int: number of altruist donors\nmax_cycle_length::Int: maximum length of a feasible exchange cycle\nmax_chain_length::Int: maximum length of a feasible exchange chain\nis_vertex_weighted::Bool: true if all weights are actually on the vertices\n\n\n\n\n\n","category":"type"},{"location":"#KidneyExchange.MIP_params","page":"KidneyExchange.jl","title":"KidneyExchange.MIP_params","text":"MIP_params\n\nMutable structure where the solving options of the compact formulation are stored\n\nFields\n\n*optimizer::String: LP and IP solver that will be used to solve the master (default is Cbc for IPs and GLPK for LPs)\n\nverbose::Bool: true if messages are printed during the solution (default = true)\nmodel_type::Mip_model: type of MIP compact model that is to be solved (default = HPIEF)\nfvs::Bool: true if a feedback vertex set is used to reduce the number of graph copies (default = true)\nreduce_vertices::Bool: true if we try deleting useless arcs in graph copies (default = true)\nreduce_arcs::Bool: true if we try deleting useless arcs in graph copies (default = true)\nsymmetry_break::Bool: true if the MIP model is modified to reduce the number of optimal solutions (default = true)\n\n\n\n\n\n","category":"type"},{"location":"#KidneyExchange.Mip_model","page":"KidneyExchange.jl","title":"KidneyExchange.Mip_model","text":"Mip_model\n\nEnumerated structure specifying the considered formulation when solving a compact formulation of the kidney exchange problem\n\nEnumerate values\n\nHPIEF: Hybrid position-indexed extended formulation\nEXTENDED_EDGE: Cycles are handled with an extended edge formulation and  chains are handled with position-indexed variables\nCYCLE_CUT: Cycles are handled with an extended edge formulation and  chains are handled with cycle cuts to avoid long cycles\nRELAXED_ARC: Arc formulation where the size of cycles and chains are not considered\n\n\n\n\n\n","category":"type"},{"location":"#KidneyExchange.PoolGenerator","page":"KidneyExchange.jl","title":"KidneyExchange.PoolGenerator","text":"Compatibility graph generator based on the following paper:  Increasing the Opportunity of Live Kidney Donation by Matching for Two and Three Way Exchanges. S. L. Saidman, Alvin Roth, Tayfun Sonmez, Utku Unver, Frank Delmonico. Transplantation, Volume 81, Number 5, March 15, 2006.\n\nThis is known colloquially as the \"Saidman Generator\".\n\n\n\n\n\n","category":"type"},{"location":"#KidneyExchange.Bellman_Ford_cycle_search-Tuple{Graphs.SimpleGraphs.SimpleDiGraph, Vector{Float64}, Int64, Int64, BitVector, Vector{Vector{Int64}}, Vector{Float64}}","page":"KidneyExchange.jl","title":"KidneyExchange.Bellman_Ford_cycle_search","text":"Bellman_Ford_cycle_search\n\nBellman-Ford style search for one positive cost cycle\n\n#Input parameters\n\ngraph::SimpleDiGraph : The directed graph with cost on each arc\narc_cost::Matrix{Float64}:\nsource::Int : The local vertex index from which starts the search\nK::Int : The maximal length of cycles\npred[k][v] contains u if u is a predecessor of v in a path of length k from source\n\n#Output Parameters\n\ncycle::Vector{Int}: the positive cycle found, [] if none\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.MIP_chain_search","page":"KidneyExchange.jl","title":"KidneyExchange.MIP_chain_search","text":"MIP_chain_search\n\nIP model with subtour elimination constraints. The constraints are the  generalized cutset inequalities (GCS) and they are added dynamically in a row generation algorithm. Refer for instance to the following reference for a presentation of the GCS. Taccari, Leonardo. « Integer Programming Formulations for the Elementary Shortest Path Problem ». European Journal of Operational Research 252, nᵒ 1 (2016).\n\n#Input parameters\n\nmip::Model: The JuMP model initialized with the flow conservation constraints and the bound on the length of the chain\ngraph::SimpleDiGraph : The directed graph with cost on each arc\nsource::Int: Index of the source vertex\nis_vertex::BitVector: For each vertex, indicates if it is in the considered subgraph\narc_cost::Vector{Float64}: Matrix of reduced costs of every arc\n\n#Output Parameters\n\nis_positivie_chain::Bool: True if a positive chain was found\nchain::Vector{Int}: Positive chain that was found\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.SparseUNOSSaidmanPoolGenerator-Tuple{Int64}","page":"KidneyExchange.jl","title":"KidneyExchange.SparseUNOSSaidmanPoolGenerator","text":"A tweak to the published Saidman generator; distributions VERY ROUGHLY mimic the UNOS pool as of April 15, 2013.  Data taken from the KPD Work Group Data Analysis - CMR - June 2013 report. @author John P. Dickerson\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.activate_branching_constraints-Tuple{JuMP.Model, KidneyExchange.TreeNode, BP_params}","page":"KidneyExchange.jl","title":"KidneyExchange.activate_branching_constraints","text":"activate_branching_constraints\n\nActivate the all the branching constraints corresponding to a given node of the branch-and-price enumeration tree.\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.bfs","page":"KidneyExchange.jl","title":"KidneyExchange.bfs","text":"bfs\n\nBreadth-first search to calculate the shortest path in terms of number of arcs from source to other vertices of graph g. The search stops when it is at level K even if there is vertex to be visited. Because we are trying to find cycles of at most K length, vertices of level >=K are never be considered in cycles.\n\nInput parameters\n\ng::AbstractGraph{T} : The graph\nsource::Int : The source vertex\nK::Int: The maximum number of arcs\nvertex_in_subgraph::Array{Bool}: Table indicating for each vertex if it is in the subgraph under consideration\n\nOutput parameters\n\ndists::Vector{Float64}: shortest distance from source to vertices\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.bfs_reverse","page":"KidneyExchange.jl","title":"KidneyExchange.bfs_reverse","text":"bfs_reverse\n\nBreadth-first search in the graph where arcs are reversed\n\nInput parameters\n\ng::AbstractGraph{T} : The graph\nsource::Int : The source vertex\nK::Int: The maximum number of arcs\nvertex_in_subgraph::Array{Bool}: Table indicating for each vertex if it is in the subgraph under consideration\n\nOutput parameters\n\ndists::Vector{Float64}: shortest distance from source to vertices\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.branch_and_price-Tuple{KidneyExchange.Instance, KidneyExchange.Graph_copies, BP_params, TimerOutput, Real}","page":"KidneyExchange.jl","title":"KidneyExchange.branch_and_price","text":"branch_and_price\n\nCore function of the KEP solution with branch-and-price. It requires a parsed instance and the description of the graph copies. The column generation model is that of Riazcos-Alvarez et al (2020), but the many improvements have been added, in particular in the solution of the subproblem.\n\n#Input parameters\n\ninstance::Instance: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains\nsubgraphs::Graph_copies: Description of the graph copies of the extended edge formulation\nbp_params::BP_params: solution parameters of the branch-and-price\ntimer::TimerOutput: a timer that will provide detail of where the computational time was spent during the branch-and-price\ntime_limit::Float64: time limit of the algorithm, including parsing and prepreprocessing times\n\n#Output parametes\n\nbp_status::BP_status:  Structure containing every relevant information on the execution of the algorithm (including the optimal solution)\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.branch_on_arc","page":"KidneyExchange.jl","title":"KidneyExchange.branch_on_arc","text":"branch_on_arc\n\nUpdate the branch-and-bound tree with two new nodes by branching on the given arc with the specified branching (only on master problem or both in master and in subproblem)\n\nInput parameters\n\narc_to_branch::Pair{Int,Int} : The arc to be branched\nmaster::Model: Master model where the branching constraints are added\nis_cg_branching::Bool: True if the branching impacts the subproblem of the coluln generation, false if it impacts only the master problem\ntree::Vector{TreeNode}: Branch-and-bound tree\ncurrent_node::TreeNode: Branch-and-bound node currently treated\ncolumn_pool::Vector{Column}: Pool of all columns in current master problem\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.branch_on_vertex","page":"KidneyExchange.jl","title":"KidneyExchange.branch_on_vertex","text":"branch_on_vertex\n\nUpdate the branch-and-bound tree with two new nodes by branching on the given vertex.\n\nInput parameters\n\nvertex_to_branch::Int : The arc to be branched\nmaster::Model: Master model where the branching constraints are added\ntree::Vector{TreeNode}: Branch-and-bound tree\ncurrent_node::TreeNode: Branch-and-bound node currently treated\ncolumn_pool::Vector{Column}: Pool of all columns in current master problem\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.build_hpief_mip","page":"KidneyExchange.jl","title":"KidneyExchange.build_hpief_mip","text":"solve_hpief_mip\n\nA compact MIP formulation originally proposed by Dickerson et al. (2016). The hybrid position indexed edge formulation (HPIEF) handles cycles and chains in the same formulation.\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.build_reduced_extended_edge_mip","page":"KidneyExchange.jl","title":"KidneyExchange.build_reduced_extended_edge_mip","text":"reduced_extended_edge_formulation\n\nA compact MIP formulation originally proposed by Constantino et al. (2013). for the cycles-only variant of the problem.\nHere it is adapted to the graph copies based on FVS and chains are considered with position-indexed variables.\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.check_used_vertices-Tuple{Vector{Int64}, Vector{Vector{Int64}}}","page":"KidneyExchange.jl","title":"KidneyExchange.check_used_vertices","text":"check_used_vertices\n\nGives a set of indices of cycle in cycles who have common vertices with cycle\n\n#Input parameters\n\ncycle::Array{Int,1}: The vertices of cycle\ncycles::Array{Array{Int,1},1}: The set of cycles\n\n#Output parametes\n\nredundant_cycle_indices::Array{Int,1}: the set of indices of cycle in cycles who have common vertices with cycle\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.compute_arc_flow-Tuple{Vector{Float64}, Vector{KidneyExchange.Column}}","page":"KidneyExchange.jl","title":"KidneyExchange.compute_arc_flow","text":"compute_arc_flow\n\nCalculates the relaxed value of the decision variable x(i,j) of arc (i->j) in A. x(i,j) whose value >0 is stored in x::Dict{Pair{Int,Int}, Float64}\n\n#Input parameters\n\nmastersol::Vector{Float64}: solution value of the master problem, ie: y[c] for c in cycles\nnode_columns::Vector{Column}: The corresponding cycles of current node\n\n#Output parameters\n\nx:: Dict{Pair{Int,Int}, Float64}: value of x_(i,j) of arc (i->j)\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.create_chain_mip-Tuple{Graphs.SimpleGraphs.SimpleDiGraph, Int64, String, Float64}","page":"KidneyExchange.jl","title":"KidneyExchange.create_chain_mip","text":"create_chain_mip\n\nInitialize the MIP model with cycle constraint generation for the optimal search of positive cost chains. Only one model is created for every copy to save a great amount of initialization time and memory. The model will then need to be modified for each graph copy to keep only the vertices of the graph and select the right source vertex\n\n#Input parameters\n\ngraph::SimpleDiGraph : The directed graph with cost on each arc\nL::Int: The maximal length of chains\noptimizer::String: Name of the MIP sover\ntime_limit::Real: Time limit of the solver\ngurobi_env: Gurobi environment if Guroib is used (avoids many messages from the solver)\n\nln(\"   . solution found: \", arcs, \", objective value: \", objectivevalue(mip))                 println(\"   . cost of other vertices = \", [vertexcost[e[2]] for e in arcs])                 println(\"   . the mip search did not find any positive cost chain\")             end #Output Parameters\n\nmip::Model: Initial JuMP model for the search of a positive chain\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.deactivate_branching_constraints-Tuple{JuMP.Model, KidneyExchange.TreeNode, BP_params}","page":"KidneyExchange.jl","title":"KidneyExchange.deactivate_branching_constraints","text":"deactivate_branching_constraints\n\nDeactivate the all the branching constraints corresponding to a given node of the branch-and-price enumeration tree.\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.drawPatientBlood_type-Tuple{KidneyExchange.PoolGenerator}","page":"KidneyExchange.jl","title":"KidneyExchange.drawPatientBlood_type","text":"Draws a random patient's blood type from the US distribution\n@return Blood_type.{O,A,B,AB}\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.generateAltruist-Tuple{KidneyExchange.PoolGenerator}","page":"KidneyExchange.jl","title":"KidneyExchange.generateAltruist","text":"Random rolls an altruistic donor (donor with no attached patient) @param ID unique identifier for the vertex @return altruistic donor vertex KPDVertexAltruist\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.generatePair-Tuple{KidneyExchange.PoolGenerator}","page":"KidneyExchange.jl","title":"KidneyExchange.generatePair","text":"Randomly rolls a patient-donor pair (possibly compatible or incompatible) @param ID unique identifier for the vertex @return a patient-donor pair KPDVertexPair\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.generatePraIncompatibility-Tuple{KidneyExchange.PoolGenerator, Bool}","page":"KidneyExchange.jl","title":"KidneyExchange.generatePraIncompatibility","text":"Randomly generates CPRA (Calculated Panel Reactive Antibody) for a patient-donor pair, using the Saidman method.  If the patient is the donor's wife, then CPRA is increased. @param isWifePatient is the patent the wife of the donor? @return scaled CPRA double value between 0 and 1.0\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.generate_abraham_benchmark-Tuple{}","page":"KidneyExchange.jl","title":"KidneyExchange.generate_abraham_benchmark","text":"generate_abraham_benchmark\n\nGenerate a benchmark using the same method as that used in Abraham, David J., Avrim Blum, et Tuomas Sandholm. « Clearing Algorithms for Barter Exchange Markets: Enabling Nationwide Kidney Exchanges ». In Proceedings of the 8th ACM Conference on Electronic Commerce, 295‑304. EC ’07. New York, NY, USA: ACM, 2007. https://doi.org/10.1145/1250910.1250954.\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.generate_complete_benchmark-Tuple{}","page":"KidneyExchange.jl","title":"KidneyExchange.generate_complete_benchmark","text":"generate_complete_benchmark\n\nGenerate all the instances that are used in addition to the PrefLib to assess the solution methods of this package in the article describing the package.\n\nadd citation later*\n\nThe instances are stored in the subdirectories heterogeneous, sparse and saidman of the data directory. Beware that if using the package with Pkg.add(), the data directory is in the directroy where the package is stored.\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.generate_heterogeneous_instance-Tuple{Int64, Int64, Int64}","page":"KidneyExchange.jl","title":"KidneyExchange.generate_heterogeneous_instance","text":"generate_heterogeneous_instance\n\nGenerate a KEP dataset and write in two text files with extensions .dat and .wmd; the graph generator is based on the following article: Kidney Exchange in Dynamic Sparse Heterogeneous Pools. Itai Ashlagi, Patrick Jaillet, Vahideh H. Manshadi. EC-2013.  (Extended abstract.)\n\nArguments\n\nnb_pairs::Int: number of pairs of incompatible donor and patient\nnb_altruists::Int: number of altruist donors\nindex::Int: index of the dataset, which will appear in the filename\n\n\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.generate_heterogeneous_kep_graph","page":"KidneyExchange.jl","title":"KidneyExchange.generate_heterogeneous_kep_graph","text":"Compatibility graph generator based on the following paper: Kidney Exchange in Dynamic Sparse Heterogeneous Pools. Itai Ashlagi, Patrick Jaillet, Vahideh H. Manshadi. EC-2013.  (Extended abstract.)  @author John P. Dickerson\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.generate_saidman_instance-Tuple{Int64, Int64, Int64}","page":"KidneyExchange.jl","title":"KidneyExchange.generate_saidman_instance","text":"generate_saidman_instance\n\nGenerate a KEP dataset and write in two text files with extensions .dat and .wmd; the graph generator is based on the following article: \"Increasing the Opportunity of Live Kidney Donation by Matching for Two and Three Way Exchanges\". S. L. Saidman, Alvin Roth, Tayfun Sonmez, Utku Unver, Frank Delmonico. Transplantation, Volume 81, Number 5, March 15, 2006.\n\nThis generator is usually refererred to as the \"Saidman Generator\".\n\nArguments\n\nnb_pairs::Int: number of pairs of incompatible donor and patient\nnb_altruists::Int: number of altruist donors\nindex::Int: index of the dataset, which will appear in the filename\n\n\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.generate_sparse_unos_instance-Tuple{Int64, Int64, Int64}","page":"KidneyExchange.jl","title":"KidneyExchange.generate_sparse_unos_instance","text":"generate_sparse_unos_instance\n\nGenerate a KEP dataset and write in two text files with extensions .dat and .wmd. This generator is obtained by applying the Saidman generator with small modifications in the probabilities of compatibility to roughtly mimic the UNOS pool as of April 15, 2013.  The corresponding data was taken from the KPD Work Group Data Analysis - CMR - June 2013 report. Seel also the following article for the motivation and description of this generator: \"Optimizing Kidney Exchange with Transplant Chains: Theory and Reality\". John P. Dickerson, Ariel D. Procaccia, Tuomas Sandholm; Proceedings of AAMAS; 2012\n\nArguments\n\nnb_pairs::Int: number of pairs of incompatible donor and patient\nnb_altruists::Int: number of altruist donors\nindex::Int: index of the dataset, which will appear in the filename\n\n\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.get_branching_arc-Tuple{Dict{Pair{Int64, Int64}, Float64}, Dict{Pair{Int64, Int64}, Float64}}","page":"KidneyExchange.jl","title":"KidneyExchange.get_branching_arc","text":"get_branching_arc\n\nFind a fractional arc to branch. The fractional arc closest to 0.5 will be selected to branch, if there is no fractional arc in the solution\n\nInput parameters\n\ncolumn_flow::Dict{Pair{Int,Int}, Float64} : The dictionary containing the nonzero flow on each arc due to the selection of columns\npief_flow::Dict{Pair{Int,Int}, Float64} : The dictionary containing the flow on each arc from the pief model for chain search (if applicable)\n\nOutput parameters\n\narc_to_branch::Pair{Int,Int} : The arc to be branched\nis_cg_branching::Bool: True if the branching impacts the subproblem of the coluln generation, false if it impacts only the master problem\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.get_branching_vertex-Tuple{Graphs.SimpleGraphs.SimpleDiGraph, Dict{Pair{Int64, Int64}, Float64}, Dict{Pair{Int64, Int64}, Float64}}","page":"KidneyExchange.jl","title":"KidneyExchange.get_branching_vertex","text":"get_branching_vertex\n\nFind a fractional vertex cover to branch. The fractional vertex closest to 0.5 will be selected to branch\n\nInput parameters\n\ncolumn_flow::Dict{Pair{Int,Int}, Float64} : The dictionary containing the nonzero flow on each arc due to the selection of columns\npief_flow::Dict{Pair{Int,Int}, Float64} : The dictionary containing the flow on each arc from the pief model for chain search (if applicable)\n\nOutput parameters\n\nvertex_to_branch::Int : The vertex to be branched on, nothing if none was found\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.get_feasible_solution-Tuple{Vector{Float64}, Vector{KidneyExchange.Column}}","page":"KidneyExchange.jl","title":"KidneyExchange.get_feasible_solution","text":"get_feasible_solution\n\nExtract a integer feasible solution from the fractional solution by conserving a set of vertex-disjoint cycles or chains\n\n#Input parameters\n\nfractional_solution::Array{Array{Float64,1},1}: The set of value of variables indicating whether or not cycles and chains are selected\nnode_columns::Array{Array{Column,1},1}: The set of cycles and chains associated with variables\n\n#Output parametes\n\nfeas_val::Real:  The objective value of the feasible solution\ncolumns::Array{Array{Int,1},1}: The set of selected cycles and chains of the feasible solution\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.initialize_column_pool","page":"KidneyExchange.jl","title":"KidneyExchange.initialize_column_pool","text":"initialize_column_pool\n\nInitialize the pool of columns for the branch-and-price by enumerating all k-cycles up to an input maximum value of k.\n\nInput parameters\n\n*graph::SimpleDiGraph: The KEP graph *column_pool::Vector{Column}: Pool of columns where initial columns will be pushed *max_cycle_length::Int: Maximum length of the cycles to be enumerated at initialization (0 if initialization is deactivated)\n\nOutput parameters\n\ncolumn_pool::Vector{Column} : set of columns initially added to the master problem\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.initialize_master_IP","page":"KidneyExchange.jl","title":"KidneyExchange.initialize_master_IP","text":"initialize_master_IP\n\nInitialization of the master problem with integer columns\n\n#Input parameters\n\ninstance::Instance: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains.\ncolumn_pool::Vector{Column}: the set of initial columns of the master\nbp_params::BP_params: parameters of the branch-and-price\ntime_limit::Float64: time limit for each solution of the master relaxation\n\n#Output Parameters\n\nmaster::Model: the model of the restricted master problem\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.isDonorSpouse-Tuple{KidneyExchange.PoolGenerator}","page":"KidneyExchange.jl","title":"KidneyExchange.isDonorSpouse","text":"Draws a random spousal relationship between donor and patient @return true if willing donor is patient's spouse, false otherwise\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.isPatientFemale-Tuple{KidneyExchange.PoolGenerator}","page":"KidneyExchange.jl","title":"KidneyExchange.isPatientFemale","text":"Draws a random gender from the US waitlist distribution  @return true if patient is female, false otherwise\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.isPositiveCrossmatch-Tuple{Float64}","page":"KidneyExchange.jl","title":"KidneyExchange.isPositiveCrossmatch","text":"Random roll to see if a patient and donor are crossmatch compatible @param pr_PraIncompatibility probability of a PRA-based incompatibility @return true is simulated positive crossmatch, false otherwise\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.node_master","page":"KidneyExchange.jl","title":"KidneyExchange.node_master","text":"node_master\n\nInitialization of the restricted master problem\n\n#Input parameters\n\ninstance::Instance: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains.\ncolumn_pool::Vector{Column}: the set of initial columns of the master\nbp_params::BP_params: parameters of the branch-and-price\ntime_limit::Float64: time limit for each solution of the master relaxation\n\n#Output Parameters\n\nmaster::Model: the model of the restricted master problem\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.preprocess_graph_copies","page":"KidneyExchange.jl","title":"KidneyExchange.preprocess_graph_copies","text":"preprocess_graph_copies\n\nGet the data of each graph copy used in the subproblems of the column generation. The graph copies are charaterized by a source vertex and the set of vertices (optionally set of arcs) kept in the copy. Each altruist donor is the source of one graph copy. Preprocessing the vertices that can be reached with a given maximum number of arcs does not bring any improvement to the subproblem solution, so we do nothing at this stage.\n\nFor the other copies, either each donor/patient pair is the source of one copy, or we compute a feedback vertex set (FVS) that heuristically minimizes the cardinality of the set.  Each vertex of the FVS will then be the source of one copy in the column generation subproblems. In both cases, we then compute distances from and to the source in each subgraph. These will be used in the Bellman-Ford search at each solution of the subproblems.\n\nInput parameters\n\ninstance::Instance : The structure describing the characteristics of an instance including the compatibility graph\nreduce_arcs::Bool : A boolean parameter indicating whether or not arcs that are not useful in a copy should be eliminated from it (default value is false)\nreduce_vertices::Bool : A boolean parameter indicating whether or not vertices that are not useful in copy should be eliminated from it (default value is true)\nfvs::Bool : True if the copies correspond to an FVS, false if there is one copy per donor/patient pair\n\nOutput parameters\n\nsubgraphs::Graph_copies : The structure describing the preprocessed data of the graph copies\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.print_and_check_solution","page":"KidneyExchange.jl","title":"KidneyExchange.print_and_check_solution","text":"print_and_check_solution\n\nFunction that prints the main characteristics of a solution and check that it satisfies all teh constraints\n\nInput parameters\n\ncycles::Vector{Vector{Int}}: list of cycles in the solution to display\nchains::Vector{Vector{Int}}: list of chains in the solution to display\ninstance::Instance: KEP instance whose solution it is\nverbose::Bool: true if the characteristics of the solution are printed, false if only the verification need be done\n\nOutput parameters: None\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.process_node-Tuple{KidneyExchange.TreeNode, KidneyExchange.Instance, JuMP.Model, KidneyExchange.Graph_copies, KidneyExchange.BP_status, Vector{KidneyExchange.Column}, BP_params, JuMP.Model, TimerOutput, Float64}","page":"KidneyExchange.jl","title":"KidneyExchange.process_node","text":"process_node\n\nColumns generation of the node\n\n#Input parameter\n\ntree_node::TreeNode: Branch-and-bound node to solve with column generation\n\n#Output parameter\n\ncolumn_flow::Dict{Pair{Int,Int}, Float64}: value of x_(i,j) of arc (i->j)`\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.read_kep_file-Tuple{String, String}","page":"KidneyExchange.jl","title":"KidneyExchange.read_kep_file","text":"read_kep_file\n\nContruct a KEP graph from a .wmd and a .dat input files\n\nParameters\n\nwmd_file::String : Absolute path of the .wmd file.\ndat_file::String : Absolute path of the .dat file.\n\n\n\n\n\n","category":"method"},{"location":"#KidneyExchange.relaxed_arc","page":"KidneyExchange.jl","title":"KidneyExchange.relaxed_arc","text":"relaxed_arc_deterministic\n\nDeterministic relaxed-arc formulation for a KEP problem. No cycle length limitation is imposed.\n\nParameters\n\ninstance::Instance: KEP instance to solve\nmaxtime::Real=60 : Maximum solving time in seconds\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.solve_with_BP","page":"KidneyExchange.jl","title":"KidneyExchange.solve_with_BP","text":"solve_with_BP\n\nThis is the main function to call for an execution of the branch-and-price algorithm on input file with given bounds on the length of covering cycles and chains and given options\n\n#Input parameters\n\nfilename::String: path of the input data files, this should include the name of the files, but not the .dat and .wmd extensions\nK::Int: maximum length of cycles\nL::Int: maximum length of chains\nbp_params::BP_params: solution parameters of the branch-and-price\ntimer::TimerOutput: a timer that will provide detail of where the computational time was spent during the branch-and-price\ntime_limit::Float64: time limit of the algorithm, including parsing and prepreprocessing times\n\n#Output parametes\n\ninstance::Instance: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains.\nsubgraphs::Graph_copies: Description of the graph copies of the extended edge formulation\nbp_status::BP_status:  Structure containing every relevant information on the execution of the algorithm (including the optimal solution)\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.solve_with_mip","page":"KidneyExchange.jl","title":"KidneyExchange.solve_with_mip","text":"solve_with_mip\n\nThis is the main function to call to solve the input instance with given bounds on the length of covering cycles and chains and given options\n\n#Input parameters\n\nfilename::String: path of the input data files, this should include the name of the files, but not the .dat and .wmd extensions\nK::Int: maximum length of cycles\nL::Int: Maximum length of chains\nparams::MIP_params: parameters of the MIP model and solver\ntimer::TimerOutput: a timer that will provide detail of where the computational time was spent during the branch-and-price\ntime_limit::Float64: time limit of the algorithm, including parsing and prepreprocessing times\n\n#Output parametes\n\ninstance::Instance: The parsed instance that is to be solved, it contains the KEP graph and the bounds on the length of covering cycles and chains.\nsubgraphs::Graph_copies: Description of the graph copies of the extended edge formulation\nbp_status::BP_status:  Structure containing every relevant information on the execution of the algorithm (including the optimal solution)\n\n\n\n\n\n","category":"function"},{"location":"#KidneyExchange.write_kep_file-Tuple{Graphs.SimpleGraphs.SimpleDiGraph, Matrix{Float64}, Vector{KidneyExchange.Blood_type}, Vector{KidneyExchange.Blood_type}, BitArray, Vector{Float64}, BitArray, String}","page":"KidneyExchange.jl","title":"KidneyExchange.write_kep_file","text":"write_kep_file\n\nWrite a .wmd and a .dat files to store the input KEP graph in the same form as in Preflib. Store the files in ./data directory.\n\nParameters\n\nkep_graph::SimpleDiGraph: KEP graph to write\nedge_weight::Matrix{Float64}: weight of each eadge of the KEP graph\ndonorBT::Vector{Blood_type}: blood type of the donor of each pair\npatientBT::Vector{Blood_type}: blood type of the patient of each pair\nwifeP::BitArray: for each pair true iff the donor and patient are married\npatientPRA::Vector{Float64}: PRA of the patient of each pair\nis_altruist::BitArray: for each pair, true if the donor is an altruist donor\nfile_name::String: Name of the files (before the extensions)\n\n\n\n\n\n","category":"method"}]
}
