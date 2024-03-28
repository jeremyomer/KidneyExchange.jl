"""
$(SIGNATURES)

Can a person with blood type  donorBT give a kidney to a patient with blood type patientBT?
"""
function canGiveTo(donorBT::Blood_type, patientBT::Blood_type)
    if (donorBT == O) || (patientBT == AB)
        # O can donate to {O,A,B,AB}, AB can receive from {O,A,B,AB}
        return true
    elseif (donorBT == A) && (patientBT == A)
        #  A gives to {A,AB}
        return true
    elseif (donorBT == B) && (patientBT == B)
        # B gives to {B,AB}
        return true
    else
        # O cannot receive from {A,B,AB}, A cannot receive from {B,AB}, B cannot receive from {A,AB}
        return false
    end
end

"""
$(SIGNATURES)

Can a person with blood type  patientBT receive a kidney of donorBT?
"""
function canGetFrom(patientBT::Blood_type, donorBT::Blood_type)
    return canGiveTo(donorBT, patientBT)
end

abstract type Vertex end

"""
$(TYPEDEF)
"""
struct VertexAltruist <: Vertex
    # Int identifier of the pair
    id::Int

    # Blood types for the patient and donor in the pair
    bloodTypeDonor::Blood_type

    function VertexAltruist(_id::Int, _bloodTypeDonor::Blood_type)
        return new(_id, _bloodTypeDonor)
    end
end

"""
$(TYPEDEF)
"""
struct VertexPair <: Vertex
    # Int identifier of the pair
    id::Int

    # Blood types for the patient and donor in the pair
    bloodTypePatient::Blood_type
    bloodTypeDonor::Blood_type

    # Patient's calculated probability of positive crossmatch, scaled [0,1]
    patientCPRA::Float64

    # Patient is wife of the donor (affects HLA)
    isWifePatient::Bool

    # Is the donor compatible with the patient
    isCompatible::Bool

    function VertexPair(
        _id::Int,
        _bloodTypePatient::Blood_type,
        _bloodTypeDonor::Blood_type,
        _patientCPRA::Float64,
        _isWifePatient::Bool,
        _isCompatible::Bool,
    )
        return new(
            _id,
            _bloodTypePatient,
            _bloodTypeDonor,
            _patientCPRA,
            _isWifePatient,
            _isCompatible,
        )
    end
end

# @Override
function isAltruist(v::VertexAltruist)
    return true
end
function isAltruist(v::VertexPair)
    return false
end
