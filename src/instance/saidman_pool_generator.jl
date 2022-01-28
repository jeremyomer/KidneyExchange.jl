"""

The generator code in this file is a translation of the corresponding Java functions in the code publicly shared by John Dickerson at https://github.com/JohnDickerson/KidneyExchange.
"""

"""
    PoolGenerator

Compatibility graph generator based on the following paper:
Increasing the Opportunity of Live Kidney Donation by Matching for Two and Three Way Exchanges. S. L. Saidman, Alvin Roth, Tayfun Sonmez, Utku Unver, Frank Delmonico. Transplantation, Volume 81, Number 5, March 15, 2006.

This is known colloquially as the "Saidman Generator".
"""
mutable struct PoolGenerator
	Pr_FEMALE::Float64

	Pr_SPOUSAL_DONOR::Float64

	Pr_LOW_PRA::Float64
	Pr_MED_PRA::Float64

	Pr_LOW_PRA_INCOMPATIBILITY::Float64
	Pr_MED_PRA_INCOMPATIBILITY::Float64
	Pr_HIGH_PRA_INCOMPATIBILITY::Float64

	Pr_SPOUSAL_PRA_COMPATIBILITY::Float64

	Pr_PATIENT_TYPE_O::Float64
	Pr_PATIENT_TYPE_A::Float64
	Pr_PATIENT_TYPE_B::Float64

	Pr_DONOR_TYPE_O::Float64
	Pr_DONOR_TYPE_A::Float64
	Pr_DONOR_TYPE_B::Float64

	# Current unused vertex ID for optimization graphs
	currentVertexID::Int
end

# Numbers taken from Saidman et al.'s 2006 paper "Increasing the Opportunity of Live Kidney Donation..."
function SaidmanPoolGenerator(id::Int)
	return PoolGenerator(0.4090, 0.4897,  0.7019, 0.2, 0.05, 0.45, 0.90, 0.75, 0.4814, 0.3373, 0.1428, 0.4814, 0.3373, 0.1428, id)
end


"""
    SparseUNOSSaidmanPoolGenerator(id)

A tweak to the published Saidman generator; distributions VERY ROUGHLY
mimic the UNOS pool as of April 15, 2013.  Data taken from the KPD Work
Group Data Analysis - CMR - June 2013 report.
author John P. Dickerson
"""
function SparseUNOSSaidmanPoolGenerator(id::Int)
	return PoolGenerator(0.4090, 0.4897, 0.216, 0.16, 0.50, 0.80, 0.98, 0.75, 0.651, 0.200, 0.124, 0.345, 0.459, 0.197, id)
end

"""
    drawPatientBlood_type(pool_gen)

* Draws a random patient's blood type from the US distribution
* return Blood_type.{O,A,B,AB}
"""
function drawPatientBlood_type(pool_gen::PoolGenerator)
	r = rand()

	if r <= pool_gen.Pr_PATIENT_TYPE_O
		return O
	elseif r <= pool_gen.Pr_PATIENT_TYPE_O + pool_gen.Pr_PATIENT_TYPE_A
		return A
	elseif r <= pool_gen.Pr_PATIENT_TYPE_O + pool_gen.Pr_PATIENT_TYPE_A + pool_gen.Pr_PATIENT_TYPE_B
		return B
	else
		return AB
	end
end

"""
    drawDonorBlood_type(pool_gen)

- Draws a random donor's blood type from the US distribution
- return Blood_type.{O,A,B,AB}
"""
function drawDonorBlood_type(pool_gen::PoolGenerator)
	r = rand()

	if r <= pool_gen.Pr_DONOR_TYPE_O
		return O
	elseif r <= pool_gen.Pr_DONOR_TYPE_O + pool_gen.Pr_DONOR_TYPE_A
		return A
	elseif r <= pool_gen.Pr_DONOR_TYPE_O + pool_gen.Pr_DONOR_TYPE_A + pool_gen.Pr_DONOR_TYPE_B
		return B
	else
		return AB
	end
end

"""
    isPatientFemale(pool_gen)

- Draws a random gender from the US waitlist distribution
- return true if patient is female, false otherwise

"""
function isPatientFemale(pool_gen::PoolGenerator)
	return rand() <= pool_gen.Pr_FEMALE
end

"""
    isDonorSpouse(pool_gen)

- Draws a random spousal relationship between donor and patient
- return true if willing donor is patient's spouse, false otherwise
"""
function isDonorSpouse(pool_gen::PoolGenerator)
	return rand() <= pool_gen.Pr_SPOUSAL_DONOR
end


"""
    isPositiveCrossmatch(pr_PraIncompatibility::Float64)

Random roll to see if a patient and donor are crossmatch compatible
- `pr_PraIncompatibility`: probability of a PRA-based incompatibility
- return true is simulated positive crossmatch, false otherwise
"""
function isPositiveCrossmatch(pr_PraIncompatibility::Float64)
	return rand() <= pr_PraIncompatibility
end

"""
    generatePraIncompatibility(pool_gen, isWifePatient)

Randomly generates CPRA (Calculated Panel Reactive Antibody) for a
patient-donor pair, using the Saidman method.  If the patient is the
donor's wife, then CPRA is increased.
- `isWifePatient` is the patent the wife of the donor?
- return scaled CPRA double value between 0 and 1.0
"""
function generatePraIncompatibility(pool_gen::PoolGenerator, isWifePatient::Bool)
	pr_PraIncompatibility = 0.0

	r = rand()
	if r <= pool_gen.Pr_LOW_PRA
		pr_PraIncompatibility = pool_gen.Pr_LOW_PRA_INCOMPATIBILITY;
	elseif r <= pool_gen.Pr_LOW_PRA + pool_gen.Pr_MED_PRA
		pr_PraIncompatibility = pool_gen.Pr_MED_PRA_INCOMPATIBILITY;
	else
		pr_PraIncompatibility = pool_gen.Pr_HIGH_PRA_INCOMPATIBILITY;
	end

	if !isWifePatient
		return pr_PraIncompatibility;
	else
		return 1.0 - pool_gen.Pr_SPOUSAL_PRA_COMPATIBILITY * (1.0 - pr_PraIncompatibility);
	end
end

function isCompatible(blood_type_donor::Blood_type, blood_type_patient::Blood_type, patient_cpra::Float64)
	# Donor must be blood type compatible with patient and crossmatch must be negative
	is_compatible = canGiveTo(blood_type_donor, blood_type_patient) && !isPositiveCrossmatch(patient_cpra)
	return is_compatible
end

"""
    generatePair(pool_gen)

Randomly rolls a patient-donor pair (possibly compatible or incompatible)
- `ID` unique identifier for the vertex
- return a patient-donor pair KPDVertexPair
"""
function generatePair(pool_gen::PoolGenerator)

	# draw blood types for patient and donor, along with spousal details and probability of PositiveXM
	bloodTypePatient::Blood_type = drawPatientBlood_type(pool_gen)
	bloodTypeDonor::Blood_type = drawDonorBlood_type(pool_gen)
	isWifePatient = isPatientFemale(pool_gen) && isDonorSpouse(pool_gen)
	patientCPRA = generatePraIncompatibility(pool_gen, isWifePatient);

	# can this donor donate to his or her patient? Donor must be blood type compatible with patient and crossmatch must be negative
	is_compatible = isCompatible(bloodTypeDonor, bloodTypePatient, patientCPRA)

	return bloodTypePatient, bloodTypeDonor, patientCPRA, isWifePatient, is_compatible
end

"""
    generateAltruist(pool_gen)

Random rolls an altruistic donor (donor with no attached patient)
- `ID` unique identifier for the vertex
- return altruistic donor vertex KPDVertexAltruist
"""
function generateAltruist(pool_gen::PoolGenerator)
	# Draw blood type for the altruist
	bloodTypeAltruist::Blood_type = drawDonorBlood_type(pool_gen)

	return bloodTypeAltruist
end

function generate_kep_graph(pool_gen::PoolGenerator, nb_pairs::Int, nb_altruists::Int)
	nb_vertices = nb_pairs + nb_altruists
	donorBT = Vector{Blood_type}(undef, nb_vertices)
	patientBT = Vector{Blood_type}(undef, nb_vertices)
	wifep = falses(nb_vertices)
	patientPRA = Vector{Float64}(undef, nb_vertices)
	is_altruist = falses(nb_vertices)

	# Generate enough incompatible patient-donor pair vertices
	id = 1
	while id <= nb_pairs
		bloodTypePatient, bloodTypeDonor, patientCPRA, isWifePatient, is_compatible = generatePair(pool_gen)
		if !is_compatible
			donorBT[id] = bloodTypeDonor
			patientBT[id] = bloodTypePatient
			wifep[id] = isWifePatient
			patientPRA[id] = patientCPRA
			id += 1
		end
	end

	# Generate altruistic donor vertices
	while id <= nb_vertices
		bloodTypeDonor = generateAltruist(pool_gen)
		donorBT[id] = bloodTypeDonor
		patientBT[id] = O
		wifep[id] = false
		patientPRA[id] = 0.0
		is_altruist[id] = true
		id += 1
	end

	patients_O = [v for v in 1:nb_pairs if patientBT[v] == O]
	patients_A = [v for v in 1:nb_pairs if patientBT[v] == A]
	patients_B = [v for v in 1:nb_pairs if patientBT[v] == B]
	patients_AB = [v for v in 1:nb_pairs if patientBT[v] == AB]

	# Add edges between compatible donors and other patients
	ne = 0
	in_list = Vector{Vector{Int}}(undef, nb_vertices)
	out_list = Vector{Vector{Int}}(undef, nb_vertices)
	for u in 1:nb_vertices
		in_list[u] = Vector{Int}()
		out_list[u] = Vector{Int}()
	end
	weights = zeros(nb_vertices, nb_vertices)
	for donor in 1:nb_vertices
		patients = []
		if donorBT[donor] == O
			patients = 1:nb_pairs
		elseif donorBT[donor] == A
			patients = [patients_A; patients_AB]
		elseif donorBT[donor] == B
			patients = [patients_B; patients_AB]
		else
			patients = patients_AB
		end
		for patient in patients
			if donor == patient	continue end

			if !isPositiveCrossmatch(patientPRA[patient])
				ne += 1
				push!(out_list[donor], patient)
				push!(in_list[patient], donor)
				weights[donor, patient] = 1.0
			end
		end
	end
	for u in 1:nb_vertices
		sort!(out_list[u])
		sort!(in_list[u])
	end

	return SimpleDiGraph(ne, out_list, in_list), weights, donorBT, patientBT, wifep, patientPRA, is_altruist
end

function generate_saidman_kep_graph(nb_pairs::Int, nb_altruists::Int)
	pool_gen = SaidmanPoolGenerator(1)
	return generate_kep_graph(pool_gen, nb_pairs, nb_altruists)
end

function generate_sparse_unos_kep_graph(nb_pairs::Int, nb_altruists::Int)
	pool_gen = SparseUNOSSaidmanPoolGenerator(1)
	return generate_kep_graph(pool_gen, nb_pairs, nb_altruists)
end
