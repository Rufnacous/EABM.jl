
abstract type JointType end

mutable struct Articulation{P}
    body_number::Integer
    state_indices::UnitRange{<: Integer}

    joint_type::JointType

    parent::Union{Articulation, Nothing}
    children::Array{Articulation}

    length::Real
    mass::Real
    I::Matrix{<: Real}
    pregeometry::Vector{<: Real}
    child_anchor::Vector{<: Real}
    mass_center::Vector{<: Real}

    geometry_transform::Matrix{<: Real}
    pregeometry_transform::Matrix{<: Real}

    properties::P
end

abstract type AbstractArticulatedBody end
mutable struct ArticulatedBody <: AbstractArticulatedBody
    articulation_zero::Articulation
    dofs::Integer
end
function ArticulatedBody(a1::Articulation)
    a0 = articulation_zero(); push!(a0.children, a1); a1.parent = a0;
    body = ArticulatedBody(a0, -1);
    dof(body);
    return body;
end

abstract type AbstractExternalForce end;
struct ExternalForce <: AbstractExternalForce
    force::Function
end
struct CompoundExternalForces <: AbstractExternalForce
    a::AbstractExternalForce
    b::AbstractExternalForce
end

abstract type AbstractInternalTorque end;
struct InternalTorque <: AbstractInternalTorque
    torque::Function
end
struct CompoundInternalTorques <: AbstractInternalTorque
    a::AbstractInternalTorque
    b::AbstractInternalTorque
end

mutable struct ArticulationHarness
    q::Vector{<: Real}
    q_dt::Vector{<: Real}
    q_d²t::Vector{<: Real}

    Xλ::Matrix{<: Real}
    λX::Matrix{<: Real}
    λX⁻ᵀ::Matrix{<: Real}

    IA::Matrix{<: Real}
    pA::Vector{<: Real}

    U::Matrix{<: Real}
    D::Matrix{<: Real}
    D⁻¹::Matrix{<: Real}
    u::Vector{<: Real}

    S::Matrix{<: Real}
    
    v::Vector{<: Real}
    c::Vector{<: Real}
    a::Vector{<: Real}

    p::Vector{<: Real}
    V::Vector{<: Real}
    X0::Matrix{<: Real}
end
struct StateHarness
    articulation_states::Vector{ArticulationHarness}
end

struct ArticulatedBodyAlgorithm
    passes::Vector{Tuple{Function, Symbol}}
end