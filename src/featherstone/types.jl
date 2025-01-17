
# Type definitions for the src/featherstone files.

abstract type JointType end
# A struct implementing JointType describes an articulation joint to be used in 
# Featherstone's algorithm. See Chapter 4.4 of Rigid Body Dynamics Algorithms by
# R. Featherstone (2008). A JointType requires no attributes but must have
# dof(<: JointType)::Integer
#     returning the number of degrees of freedom,
# and
# get_transformations(j<: JointType, q::Vector{<: Real}, qdt::Vector{<: Real})
#     which returns S, Ṡ, E, p. See the book or examples in joints.jl.
abstract type AbstractArticulatedBody end
# An AbstractArticulatedBody is a structure defined for use in Featherstone's
# algorithm. This abstract type is mostly future proofing. For most purposes so
# far seen, it is perfectly fine to go straight ahead and use ArticulatedBody
# (which <: AbstractArticulatedBody), for all simulation subjects.
# Must have articulation_zero::Articulation, and dofs::Integer. articulation_zero
# should probably always be created with articulation_zero() and have only one
# child.
abstract type AbstractExternalForce end;
abstract type AbstractInternalTorque end;
# AbstractExternalForce and AbstractInternalTorque facilitate the composition of
# multiple forces or multiple torques for one simulation. This probably only needs
# touching if you want to define forces that aren't applied to the articulation's
# mass centre, but elsewhere.

mutable struct Articulation{P}
# The discretizing segments in the EABM, and the 'bodies' in Featherstone's multi
# - body problem.
    body_number::Integer                  # unique index
    state_indices::UnitRange{<: Integer}  # maps vectors of body state to articulation states

    joint_type::JointType      # joint connecting this articulation to its parent

    parent::Union{Articulation, Nothing} # tree networking
    children::Array{Articulation} # that is to say, no cycles, please

    length::Real
    mass::Real
    I::Matrix{<: Real} # inertia. 6x6 spatial matrix.
    pregeometry::Vector{<: Real} # for articulations which do not sprout from their parent's child_anchor
    child_anchor::Vector{<: Real} # from where children sprout
    mass_center::Vector{<: Real}

    geometry_transform::Matrix{<: Real} # 3x3 rotation prior to joint rotation
    pregeometry_transform::Matrix{<: Real} # 3x3 rotation prior to pregeometry offset

    properties::P # generic properties struct, for use by forces to store e.g. drag coefficient
end

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

# Forces and torques are callable structs
# Compounds result from adding forces or torques
struct ExternalForce <: AbstractExternalForce
    force::Function
end
struct CompoundExternalForces <: AbstractExternalForce
    a::AbstractExternalForce
    b::AbstractExternalForce
end

struct InternalTorque <: AbstractInternalTorque
    torque::Function
end
struct CompoundInternalTorques <: AbstractInternalTorque
    a::AbstractInternalTorque
    b::AbstractInternalTorque
end

# Articulation and state harnesses store the kinematic state variables used in the ABA.
mutable struct ArticulationHarness
    q::Vector{<: Real} # generalized coordinate
    q_dt::Vector{<: Real} # rate of the aforementioned
    q_d²t::Vector{<: Real} # and the acceleration

    Xλ::Matrix{<: Real} # parent transform
    λX::Matrix{<: Real} # its inverse
    λX⁻ᵀ::Matrix{<: Real} # (λX)⁻ᵀ (See Featherstone 2008, Chapter 2.9)

    IA::Matrix{<: Real} # Articulated inertia
    pA::Vector{<: Real} # Articulated bias

    # Inertia projections / quantities
    U::Matrix{<: Real}
    D::Matrix{<: Real}
    D⁻¹::Matrix{<: Real}
    # Torque - Load 
    u::Vector{<: Real}

    # Joint motion subspace
    S::Matrix{<: Real}
    
    v::Vector{<: Real} # velocity (local)
    c::Vector{<: Real} # coriolis terms (local)
    a::Vector{<: Real} # acceleration (local)

    p::Vector{<: Real} # position (global)
    V::Vector{<: Real} # velocity (global)
    X0::Matrix{<: Real} # global transform
end
struct StateHarness
    articulation_states::Vector{ArticulationHarness}
end

# An "ArticulatedBodyAlgorithm" in this abstract definition is an algorithm
# which is applied to an articulated body over multiple recursive passes.
# Featherstones Articulated Body Algorithm is one such!
struct ArticulatedBodyAlgorithm
    passes::Vector{Tuple{Function, Symbol}}
end