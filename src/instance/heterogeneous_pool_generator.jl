"""
$(SIGNATURES)

Compatibility graph generator based on the following paper:
Kidney Exchange in Dynamic Sparse Heterogeneous Pools. Itai Ashlagi, Patrick Jaillet, Vahideh H. Manshadi. EC-2013.  (Extended abstract.)
 @author John P. Dickerson
 """
function generate_heterogeneous_kep_graph(nb_pairs::Int, nb_altruists::Int, pct_easy_to_match::Float64 = 0.5)
	# initialize the vertices data
	nb_vertices = nb_pairs + nb_altruists
	donorBT = Vector{Blood_type}(undef, nb_vertices)
	patientBT = Vector{Blood_type}(undef, nb_vertices)
	wifep = falses(nb_vertices)
	patientPRA = Vector{Float64}(undef, nb_vertices)
	is_altruist = falses(nb_vertices)

	# initialize edges data
	ne = 0
	in_list = Vector{Vector{Int}}(undef, nb_vertices)
	out_list = Vector{Vector{Int}}(undef, nb_vertices)
	for u in 1:nb_vertices
		in_list[u] = Vector{Int}()
		out_list[u] = Vector{Int}()
	end
	weights = zeros(nb_vertices, nb_vertices)

	# initialize the global characteristics of the heterogeneous pool
	num_easy_to_match = round(Int, pct_easy_to_match * nb_pairs);
	EASY_CPRA = 0.5;
	HARD_CPRA = 1.0 - 1.0/nb_pairs # Ashlagi's model uses constant/|V| for highly-sensitized probability

	# Make n1 easy-to-match vertices with low CPRA and n2 hard-to-match with high CPRA
	for id in 1:nb_pairs
		if id < num_easy_to_match
			patientPRA[id] = EASY_CPRA
		else
			patientPRA[id] = HARD_CPRA
		end
		donorBT[id] = O
		patientBT[id] = O
		wifep[id] = false
		is_altruist[id] = false
	end

	# Connect vertices randomly according to CPRA
	for donor in 1:nb_pairs
		for patient in 1:nb_pairs
			# No self-loops (assume pairs aren't compatible)
			if donor == patient continue end

			# Forms an incoming edge with probability CPRA (either high or low)
			if rand() >= patientPRA[patient]
				ne += 1
				push!(out_list[donor], patient)
				push!(in_list[patient], donor)
				weights[donor, patient] = 1.0
			end
		end
	end

	 # Add in altruists, with high probability of edges going to easy-to-match patients and low probability otherwise
	for id in nb_pairs+1:nb_pairs+nb_altruists
		donorBT[id] = O
		patientBT[id] = O
		wifep[id] = false
		patientPRA[id] = 0.0
		is_altruist[id] = true
		for patient in 1:nb_pairs
			if rand() > patientPRA[patient]
				ne += 1
				push!(out_list[id], patient)
				push!(in_list[patient], id)
				weights[id, patient] = 1.0
			end
		end
	end
	return SimpleDiGraph(ne, out_list, in_list), weights, donorBT, patientBT, wifep, patientPRA, is_altruist
end
